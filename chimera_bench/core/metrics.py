from __future__ import annotations

import gzip
import re
from pathlib import Path
from typing import Dict, Iterable, Tuple

RANKS_DEFAULT = ("species", "genus")


def _open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="ignore")
    return path.open("r", encoding="utf-8", errors="ignore")


def load_taxonomy(path: Path) -> Dict[int, Tuple[int, str]]:
    taxonomy: Dict[int, Tuple[int, str]] = {}
    with _open_text(path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            try:
                taxid = int(parts[0])
                parent = int(parts[1])
            except ValueError:
                continue
            rank = parts[2]
            taxonomy[taxid] = (parent, rank)
    return taxonomy


def taxid_to_rank(taxid: int | None, rank: str, taxonomy: Dict[int, Tuple[int, str]]):
    if taxid is None:
        return None
    seen = set()
    current = taxid
    while current not in seen:
        seen.add(current)
        info = taxonomy.get(current)
        if info is None:
            return None
        parent, cur_rank = info
        if cur_rank == rank:
            return current
        if parent == current:
            return None
        current = parent
    return None


def parse_classify_tsv(path: Path) -> Dict[str, int | None]:
    preds: Dict[str, int | None] = {}
    with _open_text(path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            read_id = parts[0]
            tokens = [t for t in parts[1:] if t]
            if not tokens or tokens[0] == "unclassified":
                preds[read_id] = None
                continue
            first = tokens[0]
            taxid_str = first.split(":", 1)[0]
            try:
                preds[read_id] = int(taxid_str)
            except ValueError:
                preds[read_id] = None
    return preds


def parse_ganon_one(path: Path) -> Dict[str, int]:
    preds: Dict[str, int] = {}
    with _open_text(path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            read_id = None
            tokens: list[str]
            if parts[0].startswith("H") and len(parts) >= 3:
                read_id = parts[1]
                tokens = parts[2:]
            else:
                read_id = parts[0]
                tokens = parts[1:]
            taxid = None
            for tok in tokens:
                tok = tok.split(":", 1)[0]
                if tok.isdigit():
                    taxid = int(tok)
                    break
            if read_id and taxid is not None:
                preds[read_id] = taxid
    return preds


def load_cami_mapping(paths: Iterable[Path]):
    truth_map: Dict[str, int] = {}
    abundance: Dict[int, int] = {}
    for path in paths:
        with _open_text(path) as fh:
            for raw in fh:
                line = raw.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 5:
                    continue
                contig_id = parts[0]
                try:
                    taxid = int(parts[2])
                except ValueError:
                    continue
                try:
                    reads = int(parts[4])
                except ValueError:
                    reads = 1
                truth_map[contig_id] = taxid
                abundance[taxid] = abundance.get(taxid, 0) + max(1, reads)
    return truth_map, abundance


def parse_tre_counts(path: Path, ranks: Iterable[str]) -> Dict[str, Dict[int, int]]:
    out: Dict[str, Dict[int, int]] = {r: {} for r in ranks}
    ranks_set = set(ranks)
    with _open_text(path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            rank = parts[0]
            if rank not in ranks_set:
                continue
            taxid_str = parts[1]
            try:
                taxid = int(taxid_str)
            except ValueError:
                continue
            try:
                count = int(float(parts[-2]))
            except (ValueError, IndexError):
                continue
            if count <= 0:
                continue
            bucket = out.setdefault(rank, {})
            bucket[taxid] = bucket.get(taxid, 0) + count
    return out


def _safe_prf(tp: int, fp: int, fn: int):
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1 = (2 * precision * recall / (precision + recall)) if (precision + recall) > 0 else 0.0
    return precision, recall, f1


def compute_per_read_metrics(
    truth: Dict[str, int],
    preds: Dict[str, int | None],
    taxonomy: Dict[int, Tuple[int, str]],
    ranks: Iterable[str],
):
    metrics: Dict[str, float] = {}
    total = len(truth)
    classified = 0
    for read_id, true_taxid in truth.items():
        pred_taxid = preds.get(read_id)
        if pred_taxid is not None:
            classified += 1
    if total > 0:
        metrics["per_read_classified_rate"] = classified / total
        metrics["per_read_unclassified_rate"] = (total - classified) / total

    for rank in ranks:
        tp = fp = fn = 0
        for read_id, true_taxid in truth.items():
            true_rank = taxid_to_rank(true_taxid, rank, taxonomy)
            pred_taxid = preds.get(read_id)
            if pred_taxid is None:
                fn += 1
                continue
            pred_rank = taxid_to_rank(pred_taxid, rank, taxonomy)
            if pred_rank is not None and true_rank is not None and pred_rank == true_rank:
                tp += 1
            else:
                fp += 1
                fn += 1
        precision, recall, f1 = _safe_prf(tp, fp, fn)
        metrics[f"per_read_precision_{rank}"] = precision
        metrics[f"per_read_recall_{rank}"] = recall
        metrics[f"per_read_f1_{rank}"] = f1
    return metrics


def compute_abundance_metrics(
    truth: Dict[str, Dict[int, int]],
    preds: Dict[str, Dict[int, int]],
    ranks: Iterable[str],
    presence_tau: float = 0.0,
):
    metrics: Dict[str, float] = {}
    for rank in ranks:
        truth_rank = truth.get(rank)
        pred_rank = preds.get(rank)
        if not truth_rank and not pred_rank:
            continue
        truth_vec = truth_rank or {}
        pred_vec = pred_rank or {}
        truth_total = sum(truth_vec.values())
        pred_total = sum(pred_vec.values())
        if truth_total > 0:
            truth_norm = {k: v / truth_total for k, v in truth_vec.items()}
        else:
            truth_norm = {}
        if pred_total > 0:
            pred_norm = {k: v / pred_total for k, v in pred_vec.items()}
        else:
            pred_norm = {}
        keys = set(truth_norm) | set(pred_norm)
        l1 = sum(abs(pred_norm.get(k, 0.0) - truth_norm.get(k, 0.0)) for k in keys)
        metrics[f"abundance_l1_{rank}"] = l1
        metrics[f"abundance_tv_{rank}"] = 0.5 * l1
        metrics[f"abundance_bc_{rank}"] = 0.5 * l1

        tp = fp = fn = 0
        for k in keys:
            t_val = truth_norm.get(k, 0.0)
            p_val = pred_norm.get(k, 0.0)
            t_present = t_val > presence_tau
            p_present = p_val > presence_tau
            if t_present and p_present:
                tp += 1
            elif (not t_present) and p_present:
                fp += 1
            elif t_present and (not p_present):
                fn += 1
        precision, recall, f1 = _safe_prf(tp, fp, fn)
        metrics[f"presence_precision_{rank}"] = precision
        metrics[f"presence_recall_{rank}"] = recall
        metrics[f"presence_f1_{rank}"] = f1
    return metrics


def _resolve_taxonomy(exp: dict) -> Path | None:
    tax_path = exp.get("taxonomy") or exp.get("tax_path")
    if tax_path:
        path = Path(tax_path)
        if path.exists():
            return path
    db_prefix = exp.get("db") or exp.get("db_prefix")
    if db_prefix:
        path = Path(f"{db_prefix}.tax")
        if path.exists():
            return path
    return None


def _resolve_mapping_paths(dataset: dict) -> list[Path]:
    direct = dataset.get("truth_map") or dataset.get("truth_mapping") or dataset.get("truth_maps")
    if direct:
        if isinstance(direct, (str, Path)):
            return [Path(direct)]
        return [Path(p) for p in direct]

    truth_dir = dataset.get("truth_dir")
    reads = dataset.get("reads") or []
    if truth_dir:
        truth_dir = Path(truth_dir)
        sample_id = None
        for read_path in reads:
            match = re.search(r"sample_(\d+)", str(read_path))
            if match:
                sample_id = match.group(1)
                break
        if sample_id:
            candidates = list(truth_dir.rglob(f"*sample_{sample_id}*gsa_mapping.tsv"))
            if not candidates:
                candidates = list(truth_dir.rglob(f"*sample_{sample_id}*gsa_mapping.tsv.gz"))
            if candidates:
                return candidates

    if reads:
        read_path = Path(reads[0])
        if "fasta" in read_path.parts:
            mapping_path = Path(str(read_path).replace("/fasta/", "/mapping/"))
            mapping_path = Path(str(mapping_path).replace("_anonymous_gsa.fasta", "_gsa_mapping.tsv"))
            if mapping_path.exists():
                return [mapping_path]
            gz = Path(str(mapping_path) + ".gz")
            if gz.exists():
                return [gz]
    return []


def evaluate_with_truth(exp: dict, dataset: dict, outputs: dict) -> Dict[str, float]:
    ranks = tuple(exp.get("ranks", RANKS_DEFAULT))
    mapping_paths = _resolve_mapping_paths(dataset)
    if not mapping_paths:
        return {}
    taxonomy_path = _resolve_taxonomy(exp)
    if taxonomy_path is None:
        return {}
    taxonomy = load_taxonomy(taxonomy_path)
    truth_reads, truth_abundance = load_cami_mapping(mapping_paths)

    metrics: Dict[str, float] = {}

    classify_path = outputs.get("classify_tsv")
    preds = None
    if classify_path:
        preds = parse_classify_tsv(Path(classify_path))
    else:
        classify_one = outputs.get("classify_one")
        if classify_one:
            preds = parse_ganon_one(Path(classify_one))
    if preds is not None:
        metrics.update(compute_per_read_metrics(truth_reads, preds, taxonomy, ranks))

    truth_by_rank: Dict[str, Dict[int, int]] = {r: {} for r in ranks}
    for taxid, count in truth_abundance.items():
        for rank in ranks:
            mapped = taxid_to_rank(taxid, rank, taxonomy)
            if mapped is None:
                continue
            bucket = truth_by_rank.setdefault(rank, {})
            bucket[mapped] = bucket.get(mapped, 0) + count

    pred_source = None
    if outputs.get("report_abundance_tre"):
        pred_source = Path(outputs["report_abundance_tre"])
    elif outputs.get("report_reads_tre"):
        pred_source = Path(outputs["report_reads_tre"])
    if pred_source and pred_source.exists():
        pred_by_rank = parse_tre_counts(pred_source, ranks)
        metrics.update(compute_abundance_metrics(truth_by_rank, pred_by_rank, ranks))

    return metrics
