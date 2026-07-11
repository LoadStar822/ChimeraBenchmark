#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import re
import subprocess
import time
from collections import Counter
from pathlib import Path
from typing import Iterable, NamedTuple, TextIO


DEFAULT_THREADS = 128
TOOL_SPECS = {
    "chimera": ("Chimera", "9000008512"),
    "centrifuger": ("Centrifuger", "9000008512"),
    "kraken2": ("Kraken2_LF01", "1900008512"),
}


class BestHit(NamedTuple):
    nmatch: int
    alen: int
    mapq: int
    targets: tuple[str, ...]

    @property
    def score(self) -> tuple[int, int]:
        return self.nmatch, self.alen


def open_maybe_gzip(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open()


def short_read_id(value: str) -> str:
    return value.split(None, 1)[0]


def _update_best(hits: dict[str, BestHit], parts: list[str]) -> None:
    if len(parts) < 12:
        return
    query = short_read_id(parts[0])
    target = parts[5]
    candidate = BestHit(int(parts[9]), int(parts[10]), int(parts[11]), (target,))
    current = hits.get(query)
    if current is None or candidate.score > current.score:
        hits[query] = candidate
    elif candidate.score == current.score:
        hits[query] = BestHit(
            current.nmatch,
            current.alen,
            max(current.mapq, candidate.mapq),
            tuple(sorted(set(current.targets) | {target})),
        )


def parse_best_paf_lines(lines: Iterable[str]) -> dict[str, BestHit]:
    hits: dict[str, BestHit] = {}
    for line in lines:
        _update_best(hits, line.rstrip("\n").split("\t"))
    return hits


def parse_best_paf(path: Path) -> dict[str, BestHit]:
    with path.open() as handle:
        return parse_best_paf_lines(handle)


def classify_clade_hits(clade_hits: dict[str, BestHit]) -> dict[str, object]:
    if not clade_hits:
        return {
            "identity_class": "unmapped",
            "best_clades": "",
            "best_nmatch": "",
            "best_alen": "",
            "best_identity": "",
            "best_mapq": "",
        }

    best_score = max(hit.score for hit in clade_hits.values())
    best = {clade: hit for clade, hit in clade_hits.items() if hit.score == best_score}
    best_clades = sorted(best)
    if best_clades == ["Fna_C2"]:
        identity_class = "strict_c2_best"
    elif "Fna_C2" in best_clades:
        identity_class = "c2_tied_best"
    else:
        identity_class = "non_c2_best"

    nmatch, alen = best_score
    return {
        "identity_class": identity_class,
        "best_clades": ",".join(best_clades),
        "best_nmatch": nmatch,
        "best_alen": alen,
        "best_identity": nmatch / alen if alen else 0.0,
        "best_mapq": max(hit.mapq for hit in best.values()),
    }


def merge_clade_hits(hits_by_clade: dict[str, dict[str, BestHit]]) -> dict[str, dict[str, object]]:
    by_read: dict[str, dict[str, BestHit]] = {}
    for clade, hits in hits_by_clade.items():
        for read_id, hit in hits.items():
            by_read.setdefault(read_id, {})[clade] = hit
    return {read_id: classify_clade_hits(hits) for read_id, hits in by_read.items()}


def _identity_counts(rows: Iterable[dict[str, object]]) -> Counter[str]:
    return Counter(str(row["identity_class"]) for row in rows)


def _metrics(counts: Counter[str], denominator: int, input_reads: int) -> dict[str, int | float]:
    strict = counts["strict_c2_best"]
    tied = counts["c2_tied_best"]
    supported = strict + tied
    non_c2 = counts["non_c2_best"]
    unmapped = counts["unmapped"]

    def rate(value: int) -> float:
        return value / denominator if denominator else 0.0

    def per_million(value: int) -> float:
        return value / input_reads * 1_000_000.0 if input_reads else 0.0

    return {
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
        "strict_c2_best_reads_per_million": per_million(strict),
        "c2_supported_reads_per_million": per_million(supported),
    }


def summarize_all_reads(
    classifications: dict[str, dict[str, object]], input_reads: int
) -> dict[str, int | float]:
    if len(classifications) > input_reads:
        raise ValueError("aligned query count exceeds FASTQ record count")
    counts = _identity_counts(classifications.values())
    counts["unmapped"] += input_reads - len(classifications)
    return {"input_reads": input_reads, **_metrics(counts, input_reads, input_reads)}


def summarize_selected_reads(
    *,
    selected_ids: set[str],
    found_ids: set[str],
    classifications: dict[str, dict[str, object]],
    input_reads: int,
) -> dict[str, int | float]:
    selected_found = selected_ids & found_ids
    rows = [classifications.get(read_id, {"identity_class": "unmapped"}) for read_id in selected_found]
    counts = _identity_counts(rows)
    return {
        "tool_c2_called_reads": len(selected_ids),
        "selected_reads": len(selected_found),
        "missing_reads": len(selected_ids - found_ids),
        **_metrics(counts, len(selected_found), input_reads),
    }


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def label_taxid(label: str) -> str:
    return label.split(":", 1)[0]


def extract_tool_c2_read_ids(tool: str, path: Path, c2_taxid: str) -> set[str]:
    read_ids: set[str] = set()
    if tool == "chimera":
        with path.open() as handle:
            for line in handle:
                parts = line.rstrip("\n").split("\t", 2)
                if len(parts) >= 2 and label_taxid(parts[1]) == c2_taxid:
                    read_ids.add(parts[0])
    elif tool == "centrifuger":
        with path.open(newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if row.get("taxID") == c2_taxid:
                    read_ids.add(row["readID"])
    elif tool == "kraken2":
        with path.open() as handle:
            for line in handle:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 3 and parts[2] == c2_taxid:
                    read_ids.add(parts[1])
    else:
        raise ValueError(f"unsupported tool: {tool}")
    normalized = {short_read_id(read_id) for read_id in read_ids}
    if len(normalized) != len(read_ids):
        raise ValueError(f"{tool} read ids are not unique after FASTQ normalization")
    return normalized


def tool_output_path(tool: str, run_dir: Path, sample: str) -> Path:
    if tool == "chimera":
        return run_dir / sample / "ChimeraClassify.tsv"
    if tool == "centrifuger":
        return run_dir / f"{sample}.centrifuger.tsv"
    if tool == "kraken2":
        return run_dir / f"{sample}.kraken.tsv"
    raise ValueError(f"unsupported tool: {tool}")


def scan_fastq(path: Path, selected_ids: set[str]) -> tuple[int, set[str]]:
    records = 0
    found: set[str] = set()
    with open_maybe_gzip(path) as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            seq = handle.readline()
            plus = handle.readline()
            qual = handle.readline()
            if not qual:
                raise ValueError(f"truncated FASTQ record in {path}")
            if not header.startswith("@") or not plus.startswith("+"):
                raise ValueError(f"invalid FASTQ record in {path} at record {records + 1}")
            records += 1
            read_id = short_read_id(header[1:].rstrip("\n"))
            if read_id in selected_ids:
                found.add(read_id)
    return records, found


def run_minimap2_best(
    reference: Path,
    reads: Path,
    *,
    threads: int,
    stderr_path: Path,
) -> tuple[dict[str, BestHit], float]:
    stderr_path.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        "minimap2",
        "-x",
        "sr",
        "-c",
        "--secondary=no",
        "-t",
        str(threads),
        str(reference),
        str(reads),
    ]
    started = time.monotonic()
    with stderr_path.open("w") as stderr:
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=stderr, text=True)
        assert process.stdout is not None
        hits = parse_best_paf_lines(process.stdout)
        return_code = process.wait()
    if return_code:
        raise subprocess.CalledProcessError(return_code, cmd)
    return hits, time.monotonic() - started


def parse_references(values: list[str]) -> dict[str, Path]:
    references: dict[str, Path] = {}
    for value in values:
        if "=" not in value:
            raise ValueError(f"reference must be CLADE=PATH: {value}")
        clade, raw_path = value.split("=", 1)
        if clade in references:
            raise ValueError(f"duplicate reference clade: {clade}")
        path = Path(raw_path)
        if not path.exists():
            raise FileNotFoundError(path)
        references[clade] = path
    if "Fna_C2" not in references:
        raise ValueError("Fna_C2 reference is required")
    return references


def sample_input_reads_from_name(row: dict[str, str], fallback: int) -> int:
    status = row.get("download_status", "")
    prefix = "validated_actual_source_reads_"
    if status.startswith(prefix):
        return int(status[len(prefix) :])
    for key in ("sample", "panel"):
        match = re.search(r"head(\d+)k", row.get(key, ""))
        if match:
            return int(match.group(1)) * 1000
    return fallback


BASE_FIELDS = [
    "sample",
    "role",
    "cohort",
    "condition",
    "paper_fna_c2_pct",
    "input_reads",
]
METRIC_FIELDS = [
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
    "strict_c2_best_reads_per_million",
    "c2_supported_reads_per_million",
]


def completed_checkpoint_samples(
    direct_rows: list[dict[str, object]],
    tool_rows: dict[str, list[dict[str, object]]],
) -> set[str]:
    expected = {str(row["sample"]) for row in direct_rows}
    for tool, rows in tool_rows.items():
        observed = {str(row["sample"]) for row in rows}
        if observed != expected:
            raise ValueError(
                f"checkpoint sample mismatch for {tool}: "
                f"direct={len(expected)} tool={len(observed)}"
            )
    return expected


def load_checkpoint(
    out_dir: Path,
) -> tuple[
    list[dict[str, object]],
    dict[str, list[dict[str, object]]],
    list[dict[str, object]],
]:
    direct_path = out_dir / "direct_alignment_sample_summary.tsv"
    direct_rows: list[dict[str, object]] = read_tsv(direct_path) if direct_path.exists() else []
    tool_rows: dict[str, list[dict[str, object]]] = {}
    for tool in TOOL_SPECS:
        path = out_dir / f"{tool}_candidate_identity_summary.tsv"
        tool_rows[tool] = read_tsv(path) if path.exists() else []
    runtime_path = out_dir / "runtime.tsv"
    runtime_rows: list[dict[str, object]] = read_tsv(runtime_path) if runtime_path.exists() else []
    completed_checkpoint_samples(direct_rows, tool_rows)
    return direct_rows, tool_rows, runtime_rows


def write_checkpoint(
    out_dir: Path,
    direct_rows: list[dict[str, object]],
    tool_rows: dict[str, list[dict[str, object]]],
    runtime_rows: list[dict[str, object]],
) -> None:
    write_tsv(
        out_dir / "direct_alignment_sample_summary.tsv",
        direct_rows,
        BASE_FIELDS + METRIC_FIELDS,
    )
    for tool, rows in tool_rows.items():
        write_tsv(
            out_dir / f"{tool}_candidate_identity_summary.tsv",
            rows,
            BASE_FIELDS
            + ["tool_c2_called_reads", "selected_reads", "missing_reads"]
            + METRIC_FIELDS,
        )
    write_tsv(
        out_dir / "runtime.tsv",
        runtime_rows,
        ["sample", "stage", "seconds", "aligned_queries"],
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Audit Fna C2 calls and all-read direct alignment with one best hit per clade."
    )
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--chimera-run-dir", required=True, type=Path)
    parser.add_argument("--centrifuger-run-dir", required=True, type=Path)
    parser.add_argument("--kraken2-run-dir", required=True, type=Path)
    parser.add_argument("--reference", action="append", required=True)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--sample", action="append", default=[])
    parser.add_argument("--sample-list", type=Path)
    parser.add_argument("--threads", type=int, default=DEFAULT_THREADS)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()

    references = parse_references(args.reference)
    requested = set(args.sample)
    if args.sample_list:
        requested.update(
            line.strip() for line in args.sample_list.read_text().splitlines() if line.strip()
        )

    manifest_rows = read_tsv(args.manifest)
    if requested:
        manifest_rows = [row for row in manifest_rows if row["sample"] in requested]
        missing = requested - {row["sample"] for row in manifest_rows}
        if missing:
            raise ValueError(f"requested samples not found in manifest: {sorted(missing)}")

    run_dirs = {
        "chimera": args.chimera_run_dir,
        "centrifuger": args.centrifuger_run_dir,
        "kraken2": args.kraken2_run_dir,
    }
    if args.resume:
        direct_rows, tool_rows, runtime_rows = load_checkpoint(args.out_dir)
        completed = completed_checkpoint_samples(direct_rows, tool_rows)
        manifest_rows = [row for row in manifest_rows if row["sample"] not in completed]
        print(f"resume completed={len(completed)} remaining={len(manifest_rows)}", flush=True)
    else:
        direct_rows = []
        tool_rows = {tool: [] for tool in TOOL_SPECS}
        runtime_rows = []

    for sample_row in manifest_rows:
        sample = sample_row["sample"]
        fastq = Path(sample_row["fastq"])
        selected_by_tool: dict[str, set[str]] = {}
        selected_union: set[str] = set()
        for tool, (_, taxid) in TOOL_SPECS.items():
            output = tool_output_path(tool, run_dirs[tool], sample)
            if not output.exists():
                raise FileNotFoundError(output)
            selected = extract_tool_c2_read_ids(tool, output, taxid)
            selected_by_tool[tool] = selected
            selected_union.update(selected)

        scan_started = time.monotonic()
        input_reads, found_union = scan_fastq(fastq, selected_union)
        scan_seconds = time.monotonic() - scan_started
        expected = sample_input_reads_from_name(sample_row, input_reads)
        if input_reads != expected:
            raise ValueError(
                f"FASTQ depth mismatch for {sample}: scanned={input_reads}, expected={expected}"
            )

        hits_by_clade: dict[str, dict[str, BestHit]] = {}
        for clade, reference in references.items():
            hits, seconds = run_minimap2_best(
                reference,
                fastq,
                threads=args.threads,
                stderr_path=args.out_dir / "logs" / sample / f"{clade}.minimap2.log",
            )
            hits_by_clade[clade] = hits
            runtime_rows.append(
                {
                    "sample": sample,
                    "stage": f"align_{clade}",
                    "seconds": f"{seconds:.6f}",
                    "aligned_queries": len(hits),
                }
            )

        classifications = merge_clade_hits(hits_by_clade)
        base = {
            "sample": sample,
            "role": sample_row["role"],
            "cohort": sample_row["cohort"],
            "condition": sample_row["condition"],
            "paper_fna_c2_pct": sample_row["paper_fna_c2_pct"],
            "input_reads": input_reads,
        }
        direct_rows.append({**base, **summarize_all_reads(classifications, input_reads)})
        runtime_rows.append(
            {
                "sample": sample,
                "stage": "scan_fastq",
                "seconds": f"{scan_seconds:.6f}",
                "aligned_queries": "",
            }
        )

        for tool, selected in selected_by_tool.items():
            tool_rows[tool].append(
                {
                    **base,
                    **summarize_selected_reads(
                        selected_ids=selected,
                        found_ids=found_union,
                        classifications=classifications,
                        input_reads=input_reads,
                    ),
                }
            )
        write_checkpoint(args.out_dir, direct_rows, tool_rows, runtime_rows)
        counts = " ".join(f"{tool}={len(ids)}" for tool, ids in selected_by_tool.items())
        print(
            f"done {sample} reads={input_reads} aligned={len(classifications)} {counts}",
            flush=True,
        )

    write_checkpoint(args.out_dir, direct_rows, tool_rows, runtime_rows)


if __name__ == "__main__":
    main()
