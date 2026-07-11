#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import re
import statistics
import subprocess
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, NamedTuple, TextIO


DEFAULT_C2_TAXID = "9000008512"
DEFAULT_THREADS = 128


class PafHit(NamedTuple):
    qname: str
    tname: str
    nmatch: int
    alen: int
    mapq: int


def open_maybe_gzip(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: Iterable[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def label_taxid(label: str) -> str:
    return label.split(":", 1)[0]


def short_read_id(read_id: str) -> str:
    return read_id.split(None, 1)[0]


def extract_c2_read_ids(classify_tsv: Path, c2_taxid: str = DEFAULT_C2_TAXID) -> set[str]:
    read_ids: set[str] = set()
    with classify_tsv.open() as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t", 2)
            if len(parts) < 2:
                continue
            if label_taxid(parts[1]) == c2_taxid:
                read_ids.add(parts[0])
    return read_ids


def write_selected_fastq(fastq: Path, selected_ids: set[str], out_fastq: Path) -> dict[str, int]:
    selected_short = {short_read_id(read_id): read_id for read_id in selected_ids}
    if len(selected_short) != len(selected_ids):
        raise ValueError("selected read ids are not unique after FASTQ-header normalization")

    seen: set[str] = set()
    scanned = 0
    out_fastq.parent.mkdir(parents=True, exist_ok=True)
    with open_maybe_gzip(fastq) as src, out_fastq.open("w") as out:
        while True:
            header = src.readline()
            if not header:
                break
            seq = src.readline()
            plus = src.readline()
            qual = src.readline()
            if not qual:
                raise ValueError(f"truncated FASTQ record in {fastq}")
            scanned += 1
            read_id = header[1:].rstrip("\n")
            sid = short_read_id(read_id)
            if sid in selected_short:
                seen.add(selected_short[sid])
                out.write(header)
                out.write(seq)
                out.write(plus)
                out.write(qual)

    return {"selected": len(seen), "missing": len(selected_ids - seen), "scanned": scanned}


def parse_fasta(path: Path) -> Iterable[tuple[str, str]]:
    name: str | None = None
    chunks: list[str] = []
    with open_maybe_gzip(path) as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name = line[1:].strip()
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def safe_token(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.:-]+", "_", value).strip("_") or "na"


def build_audit_reference(reference_audit_tsv: Path, out_fasta: Path) -> list[dict[str, str]]:
    rows = read_tsv(reference_audit_tsv)
    out_fasta.parent.mkdir(parents=True, exist_ok=True)
    records: list[dict[str, str]] = []
    with out_fasta.open("w") as out:
        for row_idx, row in enumerate(rows, start=1):
            source = Path(row["path"])
            clade = row["paper_clade_projection"]
            taxid = row["taxid"]
            ref_name = row["reference_name"]
            for seq_idx, (orig_name, seq) in enumerate(parse_fasta(source), start=1):
                record_id = (
                    f"ref{row_idx:05d}_{seq_idx:05d}"
                    f"|clade={safe_token(clade)}"
                    f"|taxid={safe_token(taxid)}"
                    f"|ref={safe_token(ref_name)}"
                    f"|orig={safe_token(orig_name)}"
                )
                out.write(f">{record_id}\n")
                for i in range(0, len(seq), 80):
                    out.write(seq[i : i + 80] + "\n")
                records.append(
                    {
                        "target_id": record_id,
                        "clade": clade,
                        "taxid": taxid,
                        "reference_name": ref_name,
                        "source_path": str(source),
                    }
                )
    return records


def run_minimap2(reference_fasta: Path, reads_fastq: Path, paf_out: Path, threads: int) -> None:
    paf_out.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        "minimap2",
        "-x",
        "sr",
        "-c",
        "--secondary=yes",
        "-N",
        "100",
        "-p",
        "0.5",
        "-t",
        str(threads),
        str(reference_fasta),
        str(reads_fastq),
    ]
    with paf_out.open("w") as out:
        subprocess.run(cmd, check=True, stdout=out)


def parse_paf(path: Path) -> Iterable[PafHit]:
    with path.open() as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 12:
                continue
            yield PafHit(
                qname=parts[0],
                tname=parts[5],
                nmatch=int(parts[9]),
                alen=int(parts[10]),
                mapq=int(parts[11]),
            )


def target_clade(target_name: str) -> str:
    for part in target_name.split("|"):
        if part.startswith("clade="):
            return part.split("=", 1)[1]
    return "unknown"


def group_paf_hits(hits: Iterable[PafHit]) -> dict[str, list[PafHit]]:
    grouped: dict[str, list[PafHit]] = defaultdict(list)
    for hit in hits:
        grouped[hit.qname].append(hit)
    return dict(grouped)


def hit_score(hit: PafHit) -> tuple[int, int]:
    return (hit.nmatch, hit.alen)


def classify_read_hits(hits: list[PafHit]) -> dict[str, object]:
    if not hits:
        return {
            "identity_class": "unmapped",
            "best_clades": "",
            "best_nmatch": "",
            "best_alen": "",
            "best_identity": "",
            "best_mapq": "",
            "best_targets": "",
        }

    best_score = max(hit_score(hit) for hit in hits)
    best_hits = [hit for hit in hits if hit_score(hit) == best_score]
    best_clades = sorted({target_clade(hit.tname) for hit in best_hits})

    if best_clades == ["Fna_C2"]:
        identity_class = "strict_c2_best"
    elif "Fna_C2" in best_clades:
        identity_class = "c2_tied_best"
    else:
        identity_class = "non_c2_best"

    nmatch, alen = best_score
    identity = nmatch / alen if alen else 0.0
    return {
        "identity_class": identity_class,
        "best_clades": ",".join(best_clades),
        "best_nmatch": nmatch,
        "best_alen": alen,
        "best_identity": f"{identity:.6f}",
        "best_mapq": max(hit.mapq for hit in best_hits),
        "best_targets": ";".join(hit.tname for hit in best_hits[:10]),
    }


def sample_expected_c2_reads(row: dict[str, str], depth: int = 3_000_000) -> float:
    return float(row["paper_fna_c2_pct"]) / 100.0 * depth


def aggregate_sample(read_rows: list[dict[str, object]], sample_row: dict[str, str], selected_stats: dict[str, int]) -> dict[str, object]:
    counts = Counter(str(row["identity_class"]) for row in read_rows)
    selected = selected_stats["selected"]
    identities = [
        float(row["best_identity"])
        for row in read_rows
        if row["best_identity"] not in ("", None)
    ]

    def rate(count: int) -> str:
        return f"{(count / selected):.6f}" if selected else "0.000000"

    strict = counts["strict_c2_best"]
    tied = counts["c2_tied_best"]
    non_c2 = counts["non_c2_best"]
    unmapped = counts["unmapped"]
    supported = strict + tied
    return {
        "sample": sample_row["sample"],
        "role": sample_row["role"],
        "cohort": sample_row["cohort"],
        "condition": sample_row["condition"],
        "paper_fna_c2_pct": sample_row["paper_fna_c2_pct"],
        "expected_c2_reads_at_depth": f"{sample_expected_c2_reads(sample_row):.6f}",
        "chimera_c2_called_reads": len(read_rows) + selected_stats["missing"],
        "selected_reads": selected,
        "missing_reads": selected_stats["missing"],
        "strict_c2_best_reads": strict,
        "c2_tied_best_reads": tied,
        "c2_supported_reads": supported,
        "non_c2_best_reads": non_c2,
        "unmapped_reads": unmapped,
        "strict_c2_best_rate": rate(strict),
        "c2_tied_best_rate": rate(tied),
        "c2_supported_rate": rate(supported),
        "non_c2_best_rate": rate(non_c2),
        "unmapped_rate": rate(unmapped),
        "mean_best_identity": f"{statistics.fmean(identities):.6f}" if identities else "",
        "median_best_identity": f"{statistics.median(identities):.6f}" if identities else "",
    }


def aggregate_by(rows: list[dict[str, object]], keys: list[str]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, ...], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[tuple(str(row[key]) for key in keys)].append(row)

    out: list[dict[str, object]] = []
    for values, group in sorted(grouped.items()):
        selected = sum(int(row["selected_reads"]) for row in group)
        strict = sum(int(row["strict_c2_best_reads"]) for row in group)
        tied = sum(int(row["c2_tied_best_reads"]) for row in group)
        supported = strict + tied
        non_c2 = sum(int(row["non_c2_best_reads"]) for row in group)
        unmapped = sum(int(row["unmapped_reads"]) for row in group)

        def rate(count: int) -> str:
            return f"{(count / selected):.6f}" if selected else "0.000000"

        record: dict[str, object] = dict(zip(keys, values))
        record.update(
            {
                    "samples": len(group),
                    "selected_reads": selected,
                    "strict_c2_best_reads": strict,
                    "c2_tied_best_reads": tied,
                    "c2_supported_reads": supported,
                    "non_c2_best_reads": non_c2,
                    "unmapped_reads": unmapped,
                    "strict_c2_best_rate": rate(strict),
                    "c2_tied_best_rate": rate(tied),
                    "c2_supported_rate": rate(supported),
                    "non_c2_best_rate": rate(non_c2),
                    "unmapped_rate": rate(unmapped),
            }
        )
        out.append(record)
    return out


def run_sample(
    sample_row: dict[str, str],
    *,
    chimera_run_dir: Path,
    reference_fasta: Path,
    out_dir: Path,
    c2_taxid: str,
    threads: int,
    keep_intermediates: bool,
) -> tuple[list[dict[str, object]], dict[str, object]]:
    sample = sample_row["sample"]
    classify_tsv = chimera_run_dir / sample / "ChimeraClassify.tsv"
    if not classify_tsv.exists():
        raise FileNotFoundError(f"missing ChimeraClassify.tsv for {sample}: {classify_tsv}")

    selected_ids = extract_c2_read_ids(classify_tsv, c2_taxid)
    fastq = Path(sample_row["fastq"])
    sample_work = out_dir / "work" / sample
    selected_fastq = sample_work / f"{sample}.chimera_c2.fastq"
    paf_path = sample_work / f"{sample}.paf"

    selected_stats = write_selected_fastq(fastq, selected_ids, selected_fastq)
    if selected_stats["selected"]:
        run_minimap2(reference_fasta, selected_fastq, paf_path, threads)
        hits_by_read = group_paf_hits(parse_paf(paf_path))
    else:
        hits_by_read = {}

    read_rows: list[dict[str, object]] = []
    for full_id in sorted(selected_ids, key=short_read_id):
        sid = short_read_id(full_id)
        identity = classify_read_hits(hits_by_read.get(sid, []))
        read_rows.append(
            {
                "sample": sample,
                "role": sample_row["role"],
                "cohort": sample_row["cohort"],
                "condition": sample_row["condition"],
                "full_read_id": full_id,
                "short_read_id": sid,
                **identity,
            }
        )

    if not keep_intermediates:
        if selected_fastq.exists():
            selected_fastq.unlink()
        if paf_path.exists():
            paf_path.unlink()

    return read_rows, aggregate_sample(read_rows, sample_row, selected_stats)


READ_FIELDS = [
    "sample",
    "role",
    "cohort",
    "condition",
    "full_read_id",
    "short_read_id",
    "identity_class",
    "best_clades",
    "best_nmatch",
    "best_alen",
    "best_identity",
    "best_mapq",
    "best_targets",
]

SAMPLE_FIELDS = [
    "sample",
    "role",
    "cohort",
    "condition",
    "paper_fna_c2_pct",
    "expected_c2_reads_at_depth",
    "chimera_c2_called_reads",
    "selected_reads",
    "missing_reads",
    "strict_c2_best_reads",
    "c2_tied_best_reads",
    "c2_supported_reads",
    "non_c2_best_reads",
    "unmapped_reads",
    "strict_c2_best_rate",
    "c2_tied_best_rate",
    "c2_supported_rate",
    "non_c2_best_rate",
    "unmapped_rate",
    "mean_best_identity",
    "median_best_identity",
]

GROUP_FIELDS = [
    "cohort",
    "samples",
    "selected_reads",
    "strict_c2_best_reads",
    "c2_tied_best_reads",
    "c2_supported_reads",
    "non_c2_best_reads",
    "unmapped_reads",
    "strict_c2_best_rate",
    "c2_tied_best_rate",
    "c2_supported_rate",
    "non_c2_best_rate",
    "unmapped_rate",
]

ROLE_FIELDS = ["role", *GROUP_FIELDS[1:]]
COHORT_ROLE_FIELDS = ["cohort", "role", *GROUP_FIELDS[1:]]


def main() -> None:
    parser = argparse.ArgumentParser(description="Audit whether Chimera Fna C2 read calls align best to paper-defined C2 references.")
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--chimera-run-dir", required=True, type=Path)
    parser.add_argument("--reference-audit", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--sample", action="append", default=[])
    parser.add_argument("--c2-taxid", default=DEFAULT_C2_TAXID)
    parser.add_argument("--threads", type=int, default=DEFAULT_THREADS)
    parser.add_argument("--keep-intermediates", action="store_true")
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    reference_fasta = args.out_dir / "reference" / "paper_clade_audit_reference.fna"
    reference_records = build_audit_reference(args.reference_audit, reference_fasta)
    write_tsv(
        args.out_dir / "tables" / "reference_records.tsv",
        reference_records,
        ["target_id", "clade", "taxid", "reference_name", "source_path"],
    )

    manifest_rows = read_tsv(args.manifest)
    if args.sample:
        requested = set(args.sample)
        manifest_rows = [row for row in manifest_rows if row["sample"] in requested]
        missing = requested - {row["sample"] for row in manifest_rows}
        if missing:
            raise ValueError(f"requested samples not found in manifest: {sorted(missing)}")

    all_read_rows: list[dict[str, object]] = []
    sample_rows: list[dict[str, object]] = []
    for row in manifest_rows:
        read_rows, sample_summary = run_sample(
            row,
            chimera_run_dir=args.chimera_run_dir,
            reference_fasta=reference_fasta,
            out_dir=args.out_dir,
            c2_taxid=args.c2_taxid,
            threads=args.threads,
            keep_intermediates=args.keep_intermediates,
        )
        all_read_rows.extend(read_rows)
        sample_rows.append(sample_summary)
        print(f"done {row['sample']} selected={sample_summary['selected_reads']} supported_rate={sample_summary['c2_supported_rate']}", flush=True)

    write_tsv(args.out_dir / "tables" / "read_identity_audit.tsv", all_read_rows, READ_FIELDS)
    write_tsv(args.out_dir / "tables" / "sample_identity_summary.tsv", sample_rows, SAMPLE_FIELDS)
    write_tsv(args.out_dir / "tables" / "cohort_identity_summary.tsv", aggregate_by(sample_rows, ["cohort"]), GROUP_FIELDS)
    write_tsv(args.out_dir / "tables" / "role_identity_summary.tsv", aggregate_by(sample_rows, ["role"]), ROLE_FIELDS)
    write_tsv(args.out_dir / "tables" / "cohort_role_identity_summary.tsv", aggregate_by(sample_rows, ["cohort", "role"]), COHORT_ROLE_FIELDS)


if __name__ == "__main__":
    main()
