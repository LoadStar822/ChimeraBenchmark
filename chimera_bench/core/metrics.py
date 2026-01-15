from __future__ import annotations

import gzip
import re
from pathlib import Path
from typing import Dict, Iterable, Tuple, Union

RANKS_DEFAULT = ("species", "genus")

TaxKey = Union[int, str]

NAME_CLASSES = {
    "synonym",
    "authority",
    "equivalent name",
    "genbank synonym",
    "genbank anamorph",
    "genbank common name",
    "common name",
    "includes",
}


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


def build_name_maps(names_path: Path):
    """
    Build taxid->scientific name and synonym->scientific name mappings from NCBI `names.dmp`.

    Mirrors `bench/scripts/eval_abundance.py` logic to mitigate naming drift.
    """
    sci_for_taxid: Dict[int, str] = {}
    raw_syn: Dict[str, int] = {}
    with names_path.open("r", encoding="utf-8", errors="ignore") as fh:
        for raw in fh:
            parts = [p.strip() for p in raw.split("|")]
            if len(parts) < 4:
                continue
            try:
                taxid = int(parts[0])
            except ValueError:
                continue
            name = parts[1]
            cls = parts[3]
            if not name:
                continue
            if cls == "scientific name":
                sci_for_taxid[taxid] = name
                continue
            if cls in NAME_CLASSES:
                raw_syn.setdefault(name, taxid)
                alt = re.sub(r"\s+\(.*?\)$", "", name)
                if alt and alt != name:
                    raw_syn.setdefault(alt, taxid)
                if "(" in name or ")" in name:
                    words = name.split()
                    if len(words) >= 2:
                        alt2 = " ".join(words[:2])
                        raw_syn.setdefault(alt2, taxid)
    syn_to_sci = {syn: sci_for_taxid[taxid] for syn, taxid in raw_syn.items() if taxid in sci_for_taxid}
    sci_names = set(sci_for_taxid.values())
    return sci_for_taxid, syn_to_sci, sci_names


def build_name_taxid_maps(sci_for_taxid: Dict[int, str]) -> Dict[str, int]:
    return {name: taxid for taxid, name in sci_for_taxid.items() if name}


def normalize_name(name: str, syn_to_sci: Dict[str, str], sci_names: set[str]) -> str:
    key = name.strip()
    if key in sci_names:
        return key
    return syn_to_sci.get(key, key)


def taxid_for_name(
    name: str,
    name_to_taxid: Dict[str, int],
    syn_to_sci: Dict[str, str],
    sci_names: set[str],
) -> int | None:
    taxid = name_to_taxid.get(name)
    if taxid is not None:
        return taxid
    sci = normalize_name(name, syn_to_sci, sci_names)
    return name_to_taxid.get(sci)


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
    pred: Dict[TaxKey, float],
    truth: Dict[TaxKey, float],
    taxonomy: Dict[int, Tuple[int, str]],
) -> Dict[TaxKey, float]:
    if not truth:
        return dict(pred)
    truth_set = set(truth.keys())
    out: Dict[TaxKey, float] = {}
    for taxid, value in pred.items():
        if not isinstance(taxid, int):
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


def parse_cami_profile(path: Path) -> Dict[str, Dict[int, float]]:
    buckets: Dict[str, Dict[int, float]] = {}
    taxid_idx = None
    rank_idx = None
    perc_idx = None
    with _open_text(path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            if line.startswith("@@"):
                tokens = line.lstrip("@").split("\t") if "\t" in line else line.lstrip("@").split()
                tokens = [t.strip().lstrip("@") for t in tokens if t.strip()]
                upper = [t.upper() for t in tokens]
                taxid_idx = upper.index("TAXID") if "TAXID" in upper else None
                rank_idx = upper.index("RANK") if "RANK" in upper else None
                perc_idx = upper.index("PERCENTAGE") if "PERCENTAGE" in upper else None
                continue
            if line.startswith("@"):
                continue
            parts = line.split("\t")
            if taxid_idx is None or rank_idx is None or perc_idx is None:
                if len(parts) < 5:
                    continue
                taxid_text = parts[0].strip()
                rank_text = parts[1].strip()
                perc_text = parts[-1].strip()
            else:
                if len(parts) <= max(taxid_idx, rank_idx, perc_idx):
                    continue
                taxid_text = parts[taxid_idx].strip()
                rank_text = parts[rank_idx].strip()
                perc_text = parts[perc_idx].strip()
            taxid = taxid_text.split(":", 1)[0].split("|", 1)[0].strip()
            if not taxid or taxid in {"0", "-", "NA"}:
                continue
            if "." in taxid:
                taxid = taxid.split(".", 1)[0]
            if not taxid.isdigit():
                continue
            try:
                taxid_int = int(taxid)
            except ValueError:
                continue
            if taxid_int <= 0:
                continue
            rank = rank_text.lower()
            if not rank or rank == "unclassified":
                continue
            try:
                value = float(perc_text)
            except ValueError:
                continue
            bucket = buckets.setdefault(rank, {})
            bucket[taxid_int] = bucket.get(taxid_int, 0.0) + value
    return buckets


def map_species_profile(
    profile: Dict[str, float],
    taxonomy: Dict[int, Tuple[int, str]],
    name_to_taxid: Dict[str, int],
    syn_to_sci: Dict[str, str],
    sci_names: set[str],
    ranks: Iterable[str],
    covered_by_rank: Dict[str, set[int]] | None = None,
):
    mapped_by_rank: Dict[str, Dict[TaxKey, float]] = {r: {} for r in ranks}
    full_by_rank: Dict[str, Dict[TaxKey, float]] = {r: {} for r in ranks}
    unmapped_count = 0
    unmapped_mass = 0.0

    for name, value in profile.items():
        taxid = taxid_for_name(name, name_to_taxid, syn_to_sci, sci_names)
        if taxid is None:
            unmapped_count += 1
            unmapped_mass += value
            key = normalize_name(name, syn_to_sci, sci_names)
            for rank in ranks:
                bucket = full_by_rank.setdefault(rank, {})
                bucket[key] = bucket.get(key, 0.0) + value
            continue
        for rank in ranks:
            mapped = taxid_to_rank(taxid, rank, taxonomy)
            if mapped is None or (
                covered_by_rank is not None and mapped not in covered_by_rank.get(rank, set())
            ):
                bucket = full_by_rank.setdefault(rank, {})
                bucket[taxid] = bucket.get(taxid, 0.0) + value
                continue
            bucket_full = full_by_rank.setdefault(rank, {})
            bucket_full[mapped] = bucket_full.get(mapped, 0.0) + value
            bucket_mapped = mapped_by_rank.setdefault(rank, {})
            bucket_mapped[mapped] = bucket_mapped.get(mapped, 0.0) + value
    return mapped_by_rank, full_by_rank, unmapped_count, unmapped_mass


def map_taxid_profile_to_rank(
    profile: Dict[int, float],
    taxonomy: Dict[int, Tuple[int, str]],
    ranks: Iterable[str],
    covered_by_rank: Dict[str, set[int]] | None = None,
) -> tuple[Dict[str, Dict[TaxKey, float]], Dict[str, Dict[TaxKey, float]]]:
    mapped_by_rank: Dict[str, Dict[TaxKey, float]] = {r: {} for r in ranks}
    full_by_rank: Dict[str, Dict[TaxKey, float]] = {r: {} for r in ranks}
    for taxid, value in profile.items():
        for rank in ranks:
            mapped = taxid_to_rank(taxid, rank, taxonomy)
            if mapped is None or (
                covered_by_rank is not None and mapped not in covered_by_rank.get(rank, set())
            ):
                bucket = full_by_rank.setdefault(rank, {})
                bucket[taxid] = bucket.get(taxid, 0.0) + value
                continue
            bucket_full = full_by_rank.setdefault(rank, {})
            bucket_full[mapped] = bucket_full.get(mapped, 0.0) + value
            bucket_mapped = mapped_by_rank.setdefault(rank, {})
            bucket_mapped[mapped] = bucket_mapped.get(mapped, 0.0) + value
    return mapped_by_rank, full_by_rank


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
            if pred_taxid is None:
                fn += 1
                continue
            if not pred_is_mapped:
                fp += 1
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


def compute_abundance_metrics(
    truth: Dict[str, Dict[TaxKey, float]],
    preds: Dict[str, Dict[TaxKey, float]],
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
            truth_pct = {k: 100.0 * v / truth_total for k, v in truth_vec.items()}
        else:
            truth_pct = {}
        if pred_total > 0:
            pred_pct = {k: 100.0 * v / pred_total for k, v in pred_vec.items()}
        else:
            pred_pct = {}
        keys = set(truth_pct) | set(pred_pct)
        # Match Chimera/bench convention:
        # - L1 is reported in "percentage points" (0..200) after renormalization
        # - Bray-Curtis = L1 / 200, which equals total variation distance under this normalization
        l1_pct = sum(abs(pred_pct.get(k, 0.0) - truth_pct.get(k, 0.0)) for k in keys)
        metrics[f"abundance_l1_{rank}"] = l1_pct
        metrics[f"abundance_tv_{rank}"] = l1_pct / 200.0
        metrics[f"abundance_bc_{rank}"] = l1_pct / 200.0

        tp = fp = fn = 0
        for k in keys:
            t_val = truth_pct.get(k, 0.0)
            p_val = pred_pct.get(k, 0.0)
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
    pred_by_rank: Dict[str, Dict[TaxKey, float]],
    truth_by_rank: Dict[str, Dict[TaxKey, float]],
    taxonomy: Dict[int, Tuple[int, str]],
    ranks: Iterable[str],
) -> Dict[str, Dict[TaxKey, float]]:
    out: Dict[str, Dict[TaxKey, float]] = {}
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
        "taxonomy_nodes_dmp",
        "nodes_dmp",
        "taxonomy_nodes",
        "coverage_nodes_dmp",
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


def _resolve_names_path(exp: dict, nodes_path: Path | None = None) -> Path | None:
    for key in (
        "taxonomy_names_dmp",
        "names_dmp",
        "taxonomy_names",
        "coverage_names_dmp",
    ):
        value = exp.get(key)
        if value:
            path = Path(value)
            if path.exists():
                return path
    tax_dir = exp.get("coverage_taxonomy_dir")
    if tax_dir:
        path = Path(tax_dir) / "names.dmp"
        if path.exists():
            return path
    if nodes_path is not None:
        path = nodes_path.parent / "names.dmp"
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
    taxonomy_db, file_to_taxid, _name_to_taxid = load_taxonomy(taxonomy_path)
    nodes_path = _resolve_nodes_path(exp)
    taxonomy = taxonomy_db
    if nodes_path is not None:
        taxonomy = load_nodes_taxonomy(nodes_path)

    covered_by_rank = None
    use_coverage_filter = bool(exp.get("use_coverage_filter") or exp.get("coverage_filter"))
    if use_coverage_filter:
        target_tsv = _resolve_coverage_target(exp)
        if target_tsv is not None:
            covered_by_rank = build_coverage_sets(target_tsv, taxonomy, ranks)

    names_path = _resolve_names_path(exp, nodes_path)
    name_to_taxid: Dict[str, int] = {}
    syn_to_sci: Dict[str, str] = {}
    sci_names: set[str] = set()
    if names_path is not None:
        sci_for_taxid, syn_to_sci, sci_names = build_name_maps(names_path)
        name_to_taxid = build_name_taxid_maps(sci_for_taxid)
    else:
        # Fallback: use DB `.tax` names directly (no synonym collapse).
        name_to_taxid = {name: taxid for (rank, name), taxid in _name_to_taxid.items() if rank == "species"}

    metrics: Dict[str, float] = {}

    mapping_paths = _resolve_mapping_paths(dataset)
    truth_reads: Dict[str, int] = {}
    truth_abundance: Dict[int, int] = {}
    if mapping_paths:
        truth_reads, truth_abundance = load_cami_mapping(mapping_paths)

        classify_path = outputs.get("classify_tsv")
        preds = None
        if classify_path:
            classify_path = Path(classify_path)
            if classify_path.exists():
                preds = parse_classify_tsv(classify_path)
            else:
                metrics["classify_tsv_missing"] = 1
        else:
            classify_one = outputs.get("classify_one")
            if classify_one:
                classify_one = Path(classify_one)
                if classify_one.exists():
                    preds = parse_ganon_one(classify_one, file_to_taxid)
                else:
                    metrics["classify_one_missing"] = 1
        if preds is not None:
            desc_metrics = compute_per_read_metrics(truth_reads, preds, taxonomy, ranks, covered_by_rank)
            exact_metrics = compute_per_read_metrics_exact(truth_reads, preds, taxonomy, ranks, covered_by_rank)
            metrics.update(desc_metrics)
            for key, value in exact_metrics.items():
                metrics[f"exact_{key}"] = value

    truth_by_rank_mapped: Dict[str, Dict[TaxKey, float]] = {r: {} for r in ranks}
    truth_by_rank_full: Dict[str, Dict[TaxKey, float]] = {r: {} for r in ranks}
    truth_profile_path = dataset.get("truth_profile") or dataset.get("truth_profile_path")
    if not truth_profile_path:
        truth_candidate = dataset.get("truth")
        if truth_candidate and str(truth_candidate).endswith(".tsv"):
            truth_profile_path = truth_candidate

    unmapped = None
    if truth_profile_path:
        truth_profile_path = Path(truth_profile_path)
        if not truth_profile_path.exists():
            metrics["truth_profile_missing"] = 1
        else:
            profile = parse_truth_profile(truth_profile_path)
            if profile:
                truth_by_rank_mapped, truth_by_rank_full, unmapped, unmapped_mass = map_species_profile(
                    profile,
                    taxonomy,
                    name_to_taxid,
                    syn_to_sci,
                    sci_names,
                    ranks,
                    covered_by_rank,
                )
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
        truth_by_rank_mapped, truth_by_rank_full = map_taxid_profile_to_rank(
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
    pred_by_rank_mapped: Dict[str, Dict[TaxKey, float]] | None = None
    pred_by_rank_full: Dict[str, Dict[TaxKey, float]] | None = None
    if pred_source and pred_source.exists():
        pred_by_rank_mapped = {r: {} for r in ranks}
        pred_by_rank_full = {r: {} for r in ranks}
        raw_pred = parse_tre_counts(pred_source, ranks)
        for rank, bucket in raw_pred.items():
            for taxid, value in bucket.items():
                mapped = taxid_to_rank(taxid, rank, taxonomy)
                if mapped is None or (
                    covered_by_rank is not None and mapped not in covered_by_rank.get(rank, set())
                ):
                    key: TaxKey = taxid
                    pred_by_rank_full[rank][key] = pred_by_rank_full[rank].get(key, 0.0) + float(value)
                    continue
                pred_by_rank_full[rank][mapped] = pred_by_rank_full[rank].get(mapped, 0.0) + float(value)
                pred_by_rank_mapped[rank][mapped] = pred_by_rank_mapped[rank].get(mapped, 0.0) + float(value)
    else:
        profile_path_str = outputs.get("profile_tsv") or outputs.get("sylph_profile_tsv")
        if profile_path_str:
            profile_path = Path(profile_path_str)
            if profile_path.exists():
                pred_species, unmapped_mass, _unmapped_count = parse_sylph_profile(
                    profile_path, file_to_taxid
                )
                if pred_species or unmapped_mass > 0:
                    pred_by_rank_mapped = {r: {} for r in ranks}
                    pred_by_rank_full = {r: {} for r in ranks}
                    if unmapped_mass > 0:
                        for rank in ranks:
                            pred_by_rank_full[rank]["unmapped"] = pred_by_rank_full[rank].get("unmapped", 0.0) + unmapped_mass
                    for taxid, value in pred_species.items():
                        for rank in ranks:
                            mapped = taxid_to_rank(taxid, rank, taxonomy)
                            if mapped is None or (
                                covered_by_rank is not None and mapped not in covered_by_rank.get(rank, set())
                            ):
                                key = taxid
                                pred_by_rank_full[rank][key] = pred_by_rank_full[rank].get(key, 0.0) + value
                                continue
                            pred_by_rank_full[rank][mapped] = pred_by_rank_full[rank].get(mapped, 0.0) + value
                            pred_by_rank_mapped[rank][mapped] = pred_by_rank_mapped[rank].get(mapped, 0.0) + value
        else:
            cami_profile_path_str = outputs.get("cami_profile_tsv") or outputs.get("taxor_profile_tsv")
            if cami_profile_path_str:
                cami_path = Path(cami_profile_path_str)
                if cami_path.exists():
                    cami_by_rank = parse_cami_profile(cami_path)
                    if cami_by_rank:
                        pred_by_rank_mapped = {r: {} for r in ranks}
                        pred_by_rank_full = {r: {} for r in ranks}
                        for rank in ranks:
                            if rank == "species":
                                candidates = ("species", "strain")
                            elif rank == "genus":
                                candidates = ("species", "strain", "genus")
                            else:
                                candidates = (rank,)
                            src_rank = next((c for c in candidates if c in cami_by_rank), None)
                            if src_rank is None:
                                continue
                            for taxid, value in cami_by_rank[src_rank].items():
                                mapped = taxid_to_rank(taxid, rank, taxonomy)
                                if mapped is None or (
                                    covered_by_rank is not None and mapped not in covered_by_rank.get(rank, set())
                                ):
                                    key: TaxKey = taxid
                                    pred_by_rank_full[rank][key] = pred_by_rank_full[rank].get(key, 0.0) + float(value)
                                    continue
                                pred_by_rank_full[rank][mapped] = pred_by_rank_full[rank].get(mapped, 0.0) + float(value)
                                pred_by_rank_mapped[rank][mapped] = pred_by_rank_mapped[rank].get(mapped, 0.0) + float(value)

    if pred_by_rank_full is not None:
        has_truth = any(truth_by_rank_full.get(rank) for rank in ranks)
        if has_truth:
            desc_pred_by_rank = collapse_pred_by_rank(pred_by_rank_full, truth_by_rank_full, taxonomy, ranks)
            desc_metrics = compute_abundance_metrics(truth_by_rank_full, desc_pred_by_rank, ranks)
            metrics.update(desc_metrics)

        if pred_by_rank_mapped is not None and any(truth_by_rank_mapped.get(rank) for rank in ranks):
            exact_metrics = compute_abundance_metrics(truth_by_rank_mapped, pred_by_rank_mapped, ranks)
            for key, value in exact_metrics.items():
                metrics[f"exact_{key}"] = value

    if metrics:
        metrics.setdefault("metric_version", "descendant-aware-v1")
    return metrics
