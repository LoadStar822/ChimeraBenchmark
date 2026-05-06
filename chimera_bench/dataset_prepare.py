from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, Iterable, Iterator, Tuple


def _iter_fastq_records(path: Path) -> Iterator[Tuple[str, str, str, str]]:
    with path.open("r", encoding="utf-8", errors="strict") as fh:
        while True:
            header = fh.readline()
            if not header:
                return
            seq = fh.readline()
            plus = fh.readline()
            qual = fh.readline()
            if not seq or not plus or not qual:
                raise ValueError(f"FASTQ record truncated: {path}")
            yield header, seq, plus, qual


def _read_name(header: str, prefix: str) -> str:
    if not header.startswith(prefix):
        raise ValueError(f"Invalid FASTQ header: {header.rstrip()}")
    return header[1:].strip().split()[0]


def _pair_base(name: str) -> tuple[str, str]:
    if name.endswith("/1"):
        return name[:-2], "/1"
    if name.endswith("/2"):
        return name[:-2], "/2"
    raise ValueError(f"FASTQ read does not end with /1 or /2: {name}")


def _manifest_is_current(manifest_path: Path, outputs: Iterable[Path], inputs: Iterable[Path]) -> bool:
    if not manifest_path.exists():
        return False
    if not all(path.exists() and path.stat().st_size > 0 for path in outputs):
        return False
    try:
        manifest = json.loads(manifest_path.read_text())
    except json.JSONDecodeError:
        return False
    return manifest.get("inputs") == [str(path) for path in inputs]


def _write_manifest(
    manifest_path: Path,
    *,
    inputs: Iterable[Path],
    outputs: Iterable[Path],
    records: int,
) -> None:
    payload = {
        "inputs": [str(path) for path in inputs],
        "outputs": [str(path) for path in outputs],
        "records": records,
    }
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(payload, indent=2))


def resolve_strain_madness_reads(root: Path, sample_ids: Iterable[int | str]) -> list[Path]:
    reads: list[Path] = []
    for sid in sample_ids:
        sample_root = root / f"strmgCAMI2_sample_{sid}"
        candidates = sorted(sample_root.glob("**/anonymous_reads.fq"))
        if not candidates:
            raise FileNotFoundError(f"anonymous_reads.fq not found under {sample_root}")
        reads.append(candidates[0])
    return reads


def prepare_dataset_inputs(dataset: Dict, *, force: bool = False) -> Dict:
    source = dataset.get("source") or {}
    if source.get("kind") == "strain_madness":
        sample_ids = [int(sid) for sid in (dataset.get("sample_ids") or source.get("sample_ids") or [])]
        if not sample_ids:
            raise ValueError("strain_madness dataset requires sample_ids")
        dataset["sample_ids"] = sample_ids
        source_root = Path(source["root"])
        source_reads = resolve_strain_madness_reads(source_root, sample_ids)

        read_type = source.get("read_type")
        if read_type == "long":
            dataset["reads"] = [str(path) for path in source_reads]
            return dataset
        if read_type == "short":
            prepare = dataset.setdefault("prepare", {})
            prepare.setdefault("kind", "strain_madness_short_paired")
            prepare["interleaved_reads"] = [str(path) for path in source_reads]
        else:
            raise ValueError(f"unsupported strain_madness read_type: {read_type}")

    prepare = dataset.get("prepare")
    if not prepare:
        return dataset

    kind = prepare.get("kind")
    if kind != "strain_madness_short_paired":
        raise ValueError(f"unsupported dataset prepare kind: {kind}")

    interleaved = [Path(path) for path in prepare.get("interleaved_reads", [])]
    paired = [Path(path) for path in dataset.get("paired", [])]
    if not interleaved:
        raise ValueError("strain_madness_short_paired requires prepare.interleaved_reads")
    if len(paired) != 2:
        raise ValueError("strain_madness_short_paired requires dataset.paired with two output files")

    manifest_path = Path(prepare.get("manifest") or (str(paired[0]) + ".manifest.json"))
    for out in paired:
        out.parent.mkdir(parents=True, exist_ok=True)

    if not force and _manifest_is_current(manifest_path, paired, interleaved):
        return dataset

    tmp_outputs = [Path(str(out) + ".tmp") for out in paired]
    for tmp in tmp_outputs:
        if tmp.exists():
            tmp.unlink()

    records = 0
    try:
        with (
            tmp_outputs[0].open("w", encoding="utf-8") as r1_out,
            tmp_outputs[1].open("w", encoding="utf-8") as r2_out,
        ):
            for source_path in interleaved:
                record_iter = _iter_fastq_records(source_path)
                while True:
                    try:
                        rec1 = next(record_iter)
                    except StopIteration:
                        break
                    try:
                        rec2 = next(record_iter)
                    except StopIteration as exc:
                        raise ValueError(f"paired FASTQ mate missing in {source_path}") from exc

                    name1 = _read_name(rec1[0], "@")
                    name2 = _read_name(rec2[0], "@")
                    base1, suffix1 = _pair_base(name1)
                    base2, suffix2 = _pair_base(name2)
                    if suffix1 != "/1" or suffix2 != "/2" or base1 != base2:
                        raise ValueError(
                            f"invalid interleaved pair order in {source_path}: {name1} then {name2}"
                        )
                    if not rec1[2].startswith("+") or not rec2[2].startswith("+"):
                        raise ValueError(f"invalid FASTQ plus line in {source_path}")

                    r1_out.write("".join(rec1))
                    r2_out.write("".join(rec2))
                    records += 2

        for tmp, final in zip(tmp_outputs, paired):
            tmp.replace(final)
    finally:
        for tmp in tmp_outputs:
            if tmp.exists():
                tmp.unlink()

    _write_manifest(manifest_path, inputs=interleaved, outputs=paired, records=records)
    return dataset
