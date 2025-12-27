from __future__ import annotations

import gzip
import re
from pathlib import Path
from typing import Dict, Iterable, Tuple

RANKS_DEFAULT = ("species", "genus")
UNK_TAXID = -1


def _open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="ignore")
    return path.open("r", encoding="utf-8", errors="ignore")


def load_taxonomy(
    path: Path,
) -> tuple[Dict[int, Tuple[int, str]], Dict[str, int], Dict[tuple[str, str], int]]:
    taxonomy: Dict[int, Tuple[int, str]] = {}
    file_to_taxid: Dict[str, int] = {}
    name_to_taxid: Dict[tuple[str, str], int] = {}
    with _open_text(path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            rank = parts[2]
            if rank == "file":
                if parts[0] and not parts[0].isdigit():
                    try:
                        file_to_taxid[parts[0]] = int(parts[1])
                    except ValueError:
                        continue
                continue
            name = parts[3] if len(parts) > 3 else ""
            try:
                taxid = int(parts[0])
                parent = int(parts[1])
            except ValueError:
                continue
            taxonomy[taxid] = (parent, rank)
            if name:
                key = (rank, name)
                if key not in name_to_taxid:
                    name_to_taxid[key] = taxid
    return taxonomy, file_to_taxid, name_to_taxid


def load_nodes_taxonomy(path: Path) -> Dict[int, Tuple[int, str]]:
    taxonomy: Dict[int, Tuple[int, str]] = {}
    with path.open("r", encoding="utf-8", errors="ignore") as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            parts = [p.strip() for p in line.split("|")]
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


def build_coverage_sets(
    target_tsv: Path,
    taxonomy: Dict[int, Tuple[int, str]],
    ranks: Iterable[str],
) -> Dict[str, set[int]]:
    covered: Dict[str, set[int]] = {r: set() for r in ranks}
    with _open_text(target_tsv) as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            try:
                taxid = int(parts[1])
            except ValueError:
                continue
            for rank in ranks:
                mapped = taxid_to_rank(taxid, rank, taxonomy)
                if mapped is None:
                    continue
                covered[rank].add(mapped)
    return covered


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


def is_descendant(taxid: int | None, ancestor: int, taxonomy: Dict[int, Tuple[int, str]]) -> bool:
    if taxid is None:
        return False
    current = taxid
    seen = set()
    while current not in seen:
        if current == ancestor:
            return True
        seen.add(current)
        info = taxonomy.get(current)
        if info is None:
            return False
        parent, _rank = info
        if parent == current:
            return False
        current = parent
    return False


def collapse_pred_to_truth(
    pred: Dict[int, float],
    truth: Dict[int, float],
    taxonomy: Dict[int, Tuple[int, str]],
) -> Dict[int, float]:
    if not truth:
        return dict(pred)
    truth_set = set(truth.keys())
    out: Dict[int, float] = {}
    for taxid, value in pred.items():
        if taxid == UNK_TAXID:
            out[taxid] = out.get(taxid, 0.0) + value
            continue
        current = taxid
        seen = set()
        target = None
        while current not in seen:
            if current in truth_set:
                target = current
                break
            seen.add(current)
            info = taxonomy.get(current)
            if info is None:
                break
            parent, _rank = info
            if parent == current:
                break
            current = parent
        key = target if target is not None else taxid
        out[key] = out.get(key, 0.0) + value
    return out


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


def parse_ganon_one(path: Path, file_to_taxid: Dict[str, int] | None = None) -> Dict[str, int | None]:
    preds: Dict[str, int] = {}
    with _open_text(path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            read_id = parts[0]
            ref_id = parts[1] if len(parts) > 1 else None
            if parts[0].startswith("H") and len(parts) >= 3:
                read_id = parts[1]
                ref_id = parts[2]
            taxid = None
            if ref_id and file_to_taxid:
                taxid = file_to_taxid.get(ref_id)
            if taxid is None:
                tokens = parts[2:] if len(parts) > 2 else parts[1:]
                for tok in tokens:
                    tok = tok.split(":", 1)[0]
                    if tok.isdigit():
                        taxid = int(tok)
                        break
            if read_id and taxid is not None:
                preds[read_id] = taxid
    return preds


def parse_truth_profile(path: Path) -> Dict[str, float]:
    entries: Dict[str, float] = {}
    header = None
    with _open_text(path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if header is None:
                if any(token.lower() in {"species", "taxon", "name"} for token in parts):
                    header = [p.lower() for p in parts]
                    continue
                header = []
            if len(parts) < 2:
                continue
            if header:
                try:
                    name_idx = header.index("species")
                except ValueError:
                    name_idx = 0
                value_idx = 1 if name_idx == 0 else 0
                name = parts[name_idx].strip()
                value_text = parts[value_idx].strip()
            else:
                name = parts[0].strip()
                value_text = parts[1].strip()
            if not name:
                continue
            try:
                value = float(value_text)
            except ValueError:
                continue
            entries[name] = entries.get(name, 0.0) + value
    if not entries:
        return {}
    total = sum(entries.values())
    if total > 1.5:
        return {k: v / 100.0 for k, v in entries.items()}
    return entries


def _candidate_genome_ids(raw: str) -> list[str]:
    name = raw.strip()
    if not name:
        return []
    base = Path(name).name
    candidates = [name, base]
    for ext in (
        ".fna.gz",
        ".fa.gz",
        ".fasta.gz",
        ".fna",
        ".fa",
        ".fasta",
    ):
        if base.endswith(ext):
            candidates.append(base[: -len(ext)])
            break
    match = re.search(r"(GCF|GCA)_\d+\.\d+", base)
    if match:
        candidates.append(match.group(0))
    return candidates


def parse_sylph_profile(
    path: Path, file_to_taxid: Dict[str, int]
) -> tuple[Dict[int, float], float, int]:
    entries: Dict[int, float] = {}
    unmapped_mass = 0.0
    unmapped_count = 0
    header = None
    genome_idx = None
    value_idx = None
    with _open_text(path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            parts = line.split("\t")
            if header is None:
                header = [p.strip() for p in parts]
                lower = [p.lower() for p in header]
                if "genome_file" in lower:
                    genome_idx = lower.index("genome_file")
                if "taxonomic_abundance" in lower:
                    value_idx = lower.index("taxonomic_abundance")
                elif "sequence_abundance" in lower:
                    value_idx = lower.index("sequence_abundance")
                continue
            if genome_idx is None or value_idx is None:
                continue
            if len(parts) <= max(genome_idx, value_idx):
                continue
            genome = parts[genome_idx].strip()
            value_text = parts[value_idx].strip()
            try:
                value = float(value_text)
            except ValueError:
                continue
            taxid = None
            for candidate in _candidate_genome_ids(genome):
                taxid = file_to_taxid.get(candidate)
                if taxid is not None:
                    break
            if taxid is None:
                unmapped_mass += value
                unmapped_count += 1
                continue
            entries[taxid] = entries.get(taxid, 0.0) + value
    if not entries and unmapped_count == 0:
        return {}, 0.0, 0
    total = sum(entries.values()) + unmapped_mass
    if total > 1.5:
        entries = {k: v / 100.0 for k, v in entries.items()}
        unmapped_mass = unmapped_mass / 100.0
    return entries, unmapped_mass, unmapped_count


def map_species_profile(
    profile: Dict[str, float],
    taxonomy: Dict[int, Tuple[int, str]],
    name_to_taxid: Dict[tuple[str, str], int],
    ranks: Iterable[str],
    covered_by_rank: Dict[str, set[int]] | None = None,
):
    truth_by_rank: Dict[str, Dict[int, float]] = {r: {} for r in ranks}
    rank_unmapped_mass: Dict[str, float] = {r: 0.0 for r in ranks}
    unmapped_count = 0
    unmapped_mass = 0.0

    for name, value in profile.items():
        taxid = name_to_taxid.get(("species", name))
        if taxid is None:
            unmapped_count += 1
            unmapped_mass += value
            for rank in ranks:
                rank_unmapped_mass[rank] += value
            continue
        for rank in ranks:
            mapped = taxid_to_rank(taxid, rank, taxonomy)
            if mapped is None:
                rank_unmapped_mass[rank] += value
                continue
            if covered_by_rank is not None and mapped not in covered_by_rank.get(rank, set()):
                rank_unmapped_mass[rank] += value
                continue
            bucket = truth_by_rank.setdefault(rank, {})
            bucket[mapped] = bucket.get(mapped, 0.0) + value
    return truth_by_rank, unmapped_count, unmapped_mass, rank_unmapped_mass


def map_taxid_profile_to_rank(
    profile: Dict[int, float],
    taxonomy: Dict[int, Tuple[int, str]],
    ranks: Iterable[str],
    covered_by_rank: Dict[str, set[int]] | None = None,
) -> tuple[Dict[str, Dict[int, float]], Dict[str, float]]:
    truth_by_rank: Dict[str, Dict[int, float]] = {r: {} for r in ranks}
    rank_unmapped_mass: Dict[str, float] = {r: 0.0 for r in ranks}
    for taxid, value in profile.items():
        for rank in ranks:
            mapped = taxid_to_rank(taxid, rank, taxonomy)
            if mapped is None:
                rank_unmapped_mass[rank] += value
                continue
            if covered_by_rank is not None and mapped not in covered_by_rank.get(rank, set()):
                rank_unmapped_mass[rank] += value
                continue
            bucket = truth_by_rank.setdefault(rank, {})
            bucket[mapped] = bucket.get(mapped, 0.0) + value
    return truth_by_rank, rank_unmapped_mass


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
    covered_by_rank: Dict[str, set[int]] | None = None,
):
    metrics: Dict[str, float] = {}
    total = len(truth)
    classified = 0
    for read_id in truth:
        if preds.get(read_id) is not None:
            classified += 1
    if total > 0:
        metrics["per_read_classified_rate"] = classified / total
        metrics["per_read_unclassified_rate"] = (total - classified) / total

    for rank in ranks:
        tp = fp = fn = 0
        truth_mapped = 0
        pred_mapped = 0
        covered = covered_by_rank.get(rank) if covered_by_rank is not None else None
        for read_id, true_taxid in truth.items():
            true_rank = taxid_to_rank(true_taxid, rank, taxonomy)
            true_is_mapped = true_rank is not None and (covered is None or true_rank in covered)
            if not true_is_mapped:
                continue
            truth_mapped += 1
            pred_taxid = preds.get(read_id)
            if pred_taxid is not None:
                pred_mapped += 1
            if pred_taxid is None:
                fn += 1
                continue
            if is_descendant(pred_taxid, true_rank, taxonomy):
                tp += 1
            else:
                fp += 1
                fn += 1
        if total > 0:
            metrics[f"per_read_truth_mapped_rate_{rank}"] = truth_mapped / total
            metrics[f"per_read_pred_mapped_rate_{rank}"] = pred_mapped / total
        if truth_mapped > 0:
            precision, recall, f1 = _safe_prf(tp, fp, fn)
            metrics[f"per_read_precision_{rank}"] = precision
            metrics[f"per_read_recall_{rank}"] = recall
            metrics[f"per_read_f1_{rank}"] = f1
    return metrics


def compute_per_read_metrics_exact(
    truth: Dict[str, int],
    preds: Dict[str, int | None],
    taxonomy: Dict[int, Tuple[int, str]],
    ranks: Iterable[str],
    covered_by_rank: Dict[str, set[int]] | None = None,
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
        truth_mapped = 0
        pred_mapped = 0
        covered = covered_by_rank.get(rank) if covered_by_rank is not None else None
        for read_id, true_taxid in truth.items():
            true_rank = taxid_to_rank(true_taxid, rank, taxonomy)
            pred_taxid = preds.get(read_id)
            pred_rank = taxid_to_rank(pred_taxid, rank, taxonomy) if pred_taxid is not None else None
            pred_is_mapped = pred_rank is not None and (covered is None or pred_rank in covered)
            if pred_is_mapped:
                pred_mapped += 1
            true_is_mapped = true_rank is not None and (covered is None or true_rank in covered)
            if not true_is_mapped:
                continue
            truth_mapped += 1
            if not pred_is_mapped:
                fn += 1
                continue
            if pred_rank is not None and true_rank is not None and pred_rank == true_rank:
                tp += 1
            else:
                fp += 1
                fn += 1
        if total > 0:
            metrics[f"per_read_truth_mapped_rate_{rank}"] = truth_mapped / total
            metrics[f"per_read_pred_mapped_rate_{rank}"] = pred_mapped / total
        if truth_mapped > 0:
            precision, recall, f1 = _safe_prf(tp, fp, fn)
            metrics[f"per_read_precision_{rank}"] = precision
            metrics[f"per_read_recall_{rank}"] = recall
            metrics[f"per_read_f1_{rank}"] = f1
    return metrics


def compute_per_read_metrics_unk(
    truth: Dict[str, int],
    preds: Dict[str, int | None],
    taxonomy: Dict[int, Tuple[int, str]],
    ranks: Iterable[str],
    covered_by_rank: Dict[str, set[int]] | None = None,
):
    metrics: Dict[str, float] = {}
    for rank in ranks:
        tp = fp = fn = 0
        covered = covered_by_rank.get(rank) if covered_by_rank is not None else None
        for read_id, true_taxid in truth.items():
            true_rank = taxid_to_rank(true_taxid, rank, taxonomy)
            if true_rank is None or (covered is not None and true_rank not in covered):
                true_rank = UNK_TAXID
            pred_taxid = preds.get(read_id)
            if pred_taxid is None:
                pred_rank = UNK_TAXID
            else:
                if true_rank != UNK_TAXID and is_descendant(pred_taxid, true_rank, taxonomy):
                    pred_rank = true_rank
                else:
                    pred_rank = taxid_to_rank(pred_taxid, rank, taxonomy)
                    if pred_rank is None or (covered is not None and pred_rank not in covered):
                        pred_rank = UNK_TAXID
            if pred_rank == true_rank:
                tp += 1
            else:
                fp += 1
                fn += 1
        precision, recall, f1 = _safe_prf(tp, fp, fn)
        metrics[f"per_read_precision_{rank}_unk"] = precision
        metrics[f"per_read_recall_{rank}_unk"] = recall
        metrics[f"per_read_f1_{rank}_unk"] = f1
    return metrics


def compute_per_read_metrics_unk_exact(
    truth: Dict[str, int],
    preds: Dict[str, int | None],
    taxonomy: Dict[int, Tuple[int, str]],
    ranks: Iterable[str],
    covered_by_rank: Dict[str, set[int]] | None = None,
):
    metrics: Dict[str, float] = {}
    for rank in ranks:
        tp = fp = fn = 0
        covered = covered_by_rank.get(rank) if covered_by_rank is not None else None
        for read_id, true_taxid in truth.items():
            true_rank = taxid_to_rank(true_taxid, rank, taxonomy)
            pred_taxid = preds.get(read_id)
            pred_rank = taxid_to_rank(pred_taxid, rank, taxonomy) if pred_taxid is not None else None
            if true_rank is None or (covered is not None and true_rank not in covered):
                true_rank = UNK_TAXID
            if pred_rank is None or (covered is not None and pred_rank not in covered):
                pred_rank = UNK_TAXID
            if pred_rank == true_rank:
                tp += 1
            else:
                fp += 1
                fn += 1
        precision, recall, f1 = _safe_prf(tp, fp, fn)
        metrics[f"per_read_precision_{rank}_unk"] = precision
        metrics[f"per_read_recall_{rank}_unk"] = recall
        metrics[f"per_read_f1_{rank}_unk"] = f1
    return metrics


def compute_abundance_metrics(
    truth: Dict[str, Dict[int, float]],
    preds: Dict[str, Dict[int, float]],
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


def collapse_pred_by_rank(
    pred_by_rank: Dict[str, Dict[int, float]],
    truth_by_rank: Dict[str, Dict[int, float]],
    taxonomy: Dict[int, Tuple[int, str]],
    ranks: Iterable[str],
) -> Dict[str, Dict[int, float]]:
    out: Dict[str, Dict[int, float]] = {}
    for rank in ranks:
        pred_rank = pred_by_rank.get(rank, {})
        truth_rank = truth_by_rank.get(rank, {})
        out[rank] = collapse_pred_to_truth(pred_rank, truth_rank, taxonomy)
    return out


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


def _resolve_nodes_path(exp: dict) -> Path | None:
    for key in (
        "coverage_nodes_dmp",
        "taxonomy_nodes_dmp",
        "nodes_dmp",
        "taxonomy_nodes",
    ):
        value = exp.get(key)
        if value:
            path = Path(value)
            if path.exists():
                return path
    tax_dir = exp.get("coverage_taxonomy_dir")
    if tax_dir:
        path = Path(tax_dir) / "nodes.dmp"
        if path.exists():
            return path
    return None


def _resolve_coverage_target(exp: dict) -> Path | None:
    for key in ("coverage_target_tsv", "coverage_target", "coverage_tsv"):
        value = exp.get(key)
        if value:
            path = Path(value)
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
    taxonomy_path = _resolve_taxonomy(exp)
    if taxonomy_path is None:
        return {}
    taxonomy_db, file_to_taxid, name_to_taxid = load_taxonomy(taxonomy_path)
    nodes_path = _resolve_nodes_path(exp)
    taxonomy = taxonomy_db
    if nodes_path is not None:
        taxonomy = load_nodes_taxonomy(nodes_path)

    covered_by_rank = None
    target_tsv = _resolve_coverage_target(exp)
    if target_tsv is not None:
        covered_by_rank = build_coverage_sets(target_tsv, taxonomy, ranks)

    metrics: Dict[str, float] = {}

    mapping_paths = _resolve_mapping_paths(dataset)
    truth_reads: Dict[str, int] = {}
    truth_abundance: Dict[int, int] = {}
    if mapping_paths:
        truth_reads, truth_abundance = load_cami_mapping(mapping_paths)

        classify_path = outputs.get("classify_tsv")
        preds = None
        if classify_path:
            preds = parse_classify_tsv(Path(classify_path))
        else:
            classify_one = outputs.get("classify_one")
            if classify_one:
                preds = parse_ganon_one(Path(classify_one), file_to_taxid)
        if preds is not None:
            desc_metrics = compute_per_read_metrics(truth_reads, preds, taxonomy, ranks, covered_by_rank)
            exact_metrics = compute_per_read_metrics_exact(truth_reads, preds, taxonomy, ranks, covered_by_rank)
            metrics.update(desc_metrics)
            for key, value in exact_metrics.items():
                metrics[f"exact_{key}"] = value

            desc_unk = compute_per_read_metrics_unk(truth_reads, preds, taxonomy, ranks, covered_by_rank)
            exact_unk = compute_per_read_metrics_unk_exact(truth_reads, preds, taxonomy, ranks, covered_by_rank)
            metrics.update(desc_unk)
            for key, value in exact_unk.items():
                metrics[f"exact_{key}"] = value

    truth_by_rank: Dict[str, Dict[int, float]] = {r: {} for r in ranks}
    truth_unmapped_mass: Dict[str, float] = {r: 0.0 for r in ranks}
    truth_profile_path = dataset.get("truth_profile") or dataset.get("truth_profile_path")
    if not truth_profile_path:
        truth_candidate = dataset.get("truth")
        if truth_candidate and str(truth_candidate).endswith(".tsv"):
            truth_profile_path = truth_candidate

    unmapped = None
    if truth_profile_path:
        profile = parse_truth_profile(Path(truth_profile_path))
        if profile:
            (
                truth_by_rank,
                unmapped,
                unmapped_mass,
                truth_unmapped_mass,
            ) = map_species_profile(profile, taxonomy, name_to_taxid, ranks, covered_by_rank)
            total_species = len(profile)
            metrics["truth_profile_species_total"] = total_species
            metrics["truth_profile_species_unmapped"] = unmapped
            metrics["truth_profile_species_mapped"] = total_species - unmapped
            if total_species > 0:
                metrics["truth_profile_species_unmapped_rate"] = unmapped / total_species
            total_mass = sum(profile.values())
            metrics["truth_profile_mass_total"] = total_mass
            metrics["truth_profile_mass_unmapped"] = unmapped_mass
            metrics["truth_profile_mass_mapped"] = total_mass - unmapped_mass
            if total_mass > 0:
                metrics["truth_profile_mass_mapped_rate"] = (total_mass - unmapped_mass) / total_mass
    elif truth_abundance:
        truth_by_rank, truth_unmapped_mass = map_taxid_profile_to_rank(
            {k: float(v) for k, v in truth_abundance.items()},
            taxonomy,
            ranks,
            covered_by_rank,
        )

    pred_source = None
    if outputs.get("report_abundance_tre"):
        pred_source = Path(outputs["report_abundance_tre"])
    elif outputs.get("report_reads_tre"):
        pred_source = Path(outputs["report_reads_tre"])
    pred_by_rank: Dict[str, Dict[int, float]] | None = None
    pred_unmapped_mass: Dict[str, float] = {r: 0.0 for r in ranks}
    if pred_source and pred_source.exists():
        pred_by_rank = {r: {} for r in ranks}
        raw_pred = parse_tre_counts(pred_source, ranks)
        for rank, bucket in raw_pred.items():
            for taxid, value in bucket.items():
                mapped = taxid_to_rank(taxid, rank, taxonomy)
                if mapped is None:
                    pred_unmapped_mass[rank] += float(value)
                    continue
                if covered_by_rank is not None and mapped not in covered_by_rank.get(rank, set()):
                    pred_unmapped_mass[rank] += float(value)
                    continue
                pred_by_rank[rank][mapped] = pred_by_rank[rank].get(mapped, 0.0) + float(value)
    elif outputs.get("profile_tsv"):
        profile_path = Path(outputs["profile_tsv"])
        if profile_path.exists():
            pred_species, unmapped_mass, _unmapped_count = parse_sylph_profile(
                profile_path, file_to_taxid
            )
            if pred_species or unmapped_mass > 0:
                pred_by_rank = {r: {} for r in ranks}
                for rank in ranks:
                    pred_unmapped_mass[rank] += unmapped_mass
                for taxid, value in pred_species.items():
                    for rank in ranks:
                        mapped = taxid_to_rank(taxid, rank, taxonomy)
                        if mapped is None:
                            pred_unmapped_mass[rank] += value
                            continue
                        if covered_by_rank is not None and mapped not in covered_by_rank.get(rank, set()):
                            pred_unmapped_mass[rank] += value
                            continue
                        bucket = pred_by_rank.setdefault(rank, {})
                        bucket[mapped] = bucket.get(mapped, 0.0) + value

    if truth_by_rank or any(v > 0 for v in truth_unmapped_mass.values()):
        for rank in ranks:
            mapped_mass = sum(truth_by_rank.get(rank, {}).values())
            unmapped_mass = truth_unmapped_mass.get(rank, 0.0)
            total_mass = mapped_mass + unmapped_mass
            metrics[f"truth_mapped_mass_{rank}"] = mapped_mass
            metrics[f"truth_unmapped_mass_{rank}"] = unmapped_mass
            if total_mass > 0:
                metrics[f"truth_mapped_mass_rate_{rank}"] = mapped_mass / total_mass

    if pred_by_rank is not None:
        for rank in ranks:
            mapped_mass = sum(pred_by_rank.get(rank, {}).values())
            unmapped_mass = pred_unmapped_mass.get(rank, 0.0)
            total_mass = mapped_mass + unmapped_mass
            metrics[f"pred_mapped_mass_{rank}"] = mapped_mass
            metrics[f"pred_unmapped_mass_{rank}"] = unmapped_mass
            if total_mass > 0:
                metrics[f"pred_mapped_mass_rate_{rank}"] = mapped_mass / total_mass

    if pred_by_rank is not None:
        has_truth_mapped = any(truth_by_rank.get(rank) for rank in ranks)
        desc_pred_by_rank = collapse_pred_by_rank(pred_by_rank, truth_by_rank, taxonomy, ranks)
        if has_truth_mapped:
            desc_metrics = compute_abundance_metrics(truth_by_rank, desc_pred_by_rank, ranks)
            exact_metrics = compute_abundance_metrics(truth_by_rank, pred_by_rank, ranks)
            metrics.update(desc_metrics)
            for key, value in exact_metrics.items():
                metrics[f"exact_{key}"] = value
        if has_truth_mapped or any(v > 0 for v in truth_unmapped_mass.values()):
            truth_unk = {r: dict(truth_by_rank.get(r, {})) for r in ranks}
            pred_unk = {r: dict(pred_by_rank.get(r, {})) for r in ranks}
            for rank in ranks:
                unk_mass = truth_unmapped_mass.get(rank, 0.0)
                if unk_mass > 0:
                    truth_unk[rank][UNK_TAXID] = truth_unk[rank].get(UNK_TAXID, 0.0) + unk_mass
                pred_unk_mass = pred_unmapped_mass.get(rank, 0.0)
                if pred_unk_mass > 0:
                    pred_unk[rank][UNK_TAXID] = pred_unk[rank].get(UNK_TAXID, 0.0) + pred_unk_mass
            desc_pred_unk = collapse_pred_by_rank(pred_unk, truth_unk, taxonomy, ranks)
            desc_unk_metrics = compute_abundance_metrics(truth_unk, desc_pred_unk, ranks)
            exact_unk_metrics = compute_abundance_metrics(truth_unk, pred_unk, ranks)
            for key, value in desc_unk_metrics.items():
                metrics[f"{key}_unk"] = value
            for key, value in exact_unk_metrics.items():
                metrics[f"exact_{key}_unk"] = value

    if metrics:
        metrics.setdefault("metric_version", "descendant-aware-v1")
    return metrics
