from __future__ import annotations

import gzip
import re
from functools import lru_cache
from pathlib import Path
from typing import Dict, Iterable, Tuple, Union

RANKS_DEFAULT = ("species", "genus")
PROFILE_RANK_PRIORITY = (
    "strain",
    "species",
    "genus",
    "family",
    "order",
    "class",
    "phylum",
    "superkingdom",
)
METRIC_VERSION = "per-read-descendant-aware-v1+profile-opal-v2"

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

_NAME_MAP_CACHE: dict[
    tuple[str, int | None],
    tuple[Dict[int, str], Dict[str, str], set[str]],
] = {}


def _open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="ignore")
    return path.open("r", encoding="utf-8", errors="ignore")


@lru_cache(maxsize=16)
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


@lru_cache(maxsize=16)
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


def _name_aliases(name: str) -> set[str]:
    clean = re.sub(r"\s+", " ", name.replace('"', "").strip())
    if not clean:
        return set()
    aliases = {clean}

    no_paren = re.sub(r"\s+\([^)]*\)", "", clean).strip()
    no_paren = re.sub(r"\s+", " ", no_paren)
    if no_paren:
        aliases.add(no_paren)

    tokens = no_paren.split()
    if len(tokens) >= 2:
        genus_raw = tokens[0]
        genus = genus_raw.strip("[]")
        epithet = tokens[1].strip("[]")
        blocked = {"sp", "sp.", "cf", "cf.", "aff", "aff.", "bacterium"}
        if (
            genus
            and epithet
            and epithet.lower() not in blocked
            and re.match(r"^[A-Za-z][A-Za-z.-]*$", genus)
            and re.match(r"^[a-z][a-z.-]*$", epithet)
        ):
            aliases.add(f"{genus} {epithet}")
            aliases.add(f"{genus_raw} {epithet}")

    return {alias for alias in aliases if alias}


def _species_taxid_for_name_alias(
    taxid: int,
    taxonomy: Dict[int, Tuple[int, str]] | None,
) -> int:
    if taxonomy is None:
        return taxid
    species_taxid = taxid_to_rank(taxid, "species", taxonomy)
    return species_taxid or taxid


def build_name_maps(
    names_path: Path,
    taxonomy: Dict[int, Tuple[int, str]] | None = None,
):
    """
    Build taxid->scientific name and synonym->scientific name mappings from NCBI `names.dmp`.

    Mirrors `bench/scripts/eval_abundance.py` logic to mitigate naming drift.
    """
    cache_key = (str(names_path.resolve()), id(taxonomy) if taxonomy is not None else None)
    cached = _NAME_MAP_CACHE.get(cache_key)
    if cached is not None:
        return cached

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
                alias_taxid = _species_taxid_for_name_alias(taxid, taxonomy)
                for alias in _name_aliases(name):
                    if alias != name:
                        raw_syn.setdefault(alias, alias_taxid)
                continue
            if cls in NAME_CLASSES:
                alias_taxid = _species_taxid_for_name_alias(taxid, taxonomy)
                for alias in _name_aliases(name):
                    raw_syn.setdefault(alias, alias_taxid)
    syn_to_sci = {syn: sci_for_taxid[taxid] for syn, taxid in raw_syn.items() if taxid in sci_for_taxid}
    sci_names = set(sci_for_taxid.values())
    result = (sci_for_taxid, syn_to_sci, sci_names)
    _NAME_MAP_CACHE[cache_key] = result
    return result


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
            read_id = normalize_read_id(parts[0])
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


def normalize_read_id(read_id: str) -> str:
    return read_id.strip().split()[0] if read_id.strip() else ""


def _paired_mate_id(read_id: str) -> str | None:
    if read_id.endswith("/1"):
        return f"{read_id[:-2]}/2"
    if read_id.endswith("/2"):
        return f"{read_id[:-2]}/1"
    return None


def _prediction_for_read(
    preds: Dict[str, int | None],
    read_id: str,
    truth: Dict[str, int] | None = None,
) -> int | None:
    normalized = normalize_read_id(read_id)
    if normalized in preds:
        return preds[normalized]
    if normalized.endswith("/1") or normalized.endswith("/2"):
        base_pred = preds.get(normalized[:-2])
        if base_pred is not None or normalized[:-2] in preds:
            return base_pred
        mate = _paired_mate_id(normalized)
        if (
            mate is not None
            and truth is not None
            and mate in preds
            and mate in truth
            and normalized in truth
            and truth[mate] == truth[normalized]
        ):
            return preds[mate]
    return None


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
    name_idx = None
    value_idx = None
    with _open_text(path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if header is None:
                lowered = [p.lower() for p in parts]
                if any(token in {"species_label", "species_name", "species", "taxon", "name"} for token in lowered):
                    header = [p.lower() for p in parts]
                    for candidate in ("species_label", "species_name", "taxon", "name", "species"):
                        if candidate in header:
                            name_idx = header.index(candidate)
                            break
                    for candidate in (
                        "truth_abundance",
                        "abundance",
                        "relative_abundance",
                        "percentage",
                        "percent",
                        "value",
                    ):
                        if candidate in header:
                            value_idx = header.index(candidate)
                            break
                    continue
                header = []
            if len(parts) < 2:
                continue
            if header:
                if name_idx is None:
                    name_idx = 0
                if value_idx is None:
                    value_idx = 1 if name_idx == 0 else 0
                if len(parts) <= max(name_idx, value_idx) or name_idx == value_idx:
                    continue
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


def load_species_label_mapping(
    paths: Iterable[Path],
    name_to_taxid: Dict[str, int],
    syn_to_sci: Dict[str, str],
    sci_names: set[str],
) -> tuple[Dict[str, int], Dict[int, int], Dict[str, int], int, int]:
    truth_map: Dict[str, int] = {}
    abundance: Dict[int, int] = {}
    read_weight: Dict[str, int] = {}
    unmapped: set[str] = set()
    unmapped_rows = 0
    for path in paths:
        with _open_text(path) as fh:
            header: list[str] | None = None
            read_idx = None
            name_idx = None
            for raw in fh:
                line = raw.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if header is None:
                    header = [p.lower() for p in parts]
                    if "read_id" not in header:
                        raise ValueError(f"species-label truth missing read_id column: {path}")
                    read_idx = header.index("read_id")
                    for candidate in ("species_label", "species_name", "taxon", "name", "species"):
                        if candidate in header:
                            name_idx = header.index(candidate)
                            break
                    if name_idx is None:
                        raise ValueError(f"species-label truth missing species label column: {path}")
                    continue
                if read_idx is None or name_idx is None or len(parts) <= max(read_idx, name_idx):
                    continue
                read_id = parts[read_idx].strip()
                species_name = parts[name_idx].strip()
                if not read_id or not species_name:
                    continue
                taxid = taxid_for_name(species_name, name_to_taxid, syn_to_sci, sci_names)
                if taxid is None:
                    unmapped.add(species_name)
                    unmapped_rows += 1
                    continue
                previous = truth_map.get(read_id)
                if previous is not None and previous != taxid:
                    raise ValueError(f"conflicting species-label truth for read {read_id}: {previous} vs {taxid}")
                if previous is None:
                    truth_map[read_id] = taxid
                    abundance[taxid] = abundance.get(taxid, 0) + 1
                    read_weight[read_id] = 1
    return truth_map, abundance, read_weight, unmapped_rows, len(unmapped)


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
    contig_weight: Dict[str, int] = {}
    for path in paths:
        with _open_text(path) as fh:
            for raw in fh:
                line = raw.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 4:
                    continue
                contig_id = parts[0]
                try:
                    taxid = int(parts[2])
                except ValueError:
                    continue
                if len(parts) >= 5:
                    try:
                        reads = int(parts[4])
                    except ValueError:
                        reads = 1
                else:
                    reads = 1
                weight = max(1, reads)
                truth_map[contig_id] = taxid
                contig_weight[contig_id] = weight
                abundance[taxid] = abundance.get(taxid, 0) + weight
    return truth_map, abundance, contig_weight


def parse_tre_counts(path: Path, ranks: Iterable[str] | None = None) -> Dict[str, Dict[int, float]]:
    out: Dict[str, Dict[int, float]] = {}
    ranks_set = set(ranks) if ranks is not None else None
    with _open_text(path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            rank = parts[0]
            if ranks_set is not None and rank not in ranks_set:
                continue
            taxid_str = parts[1]
            try:
                taxid = int(taxid_str)
            except ValueError:
                continue
            try:
                count = float(parts[-2])
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


def _normalize_profile(profile: Dict[int, float]) -> Dict[int, float]:
    cleaned = {taxid: float(value) for taxid, value in profile.items() if value > 0}
    total = sum(cleaned.values())
    if total <= 0:
        return {}
    return {taxid: value / total for taxid, value in cleaned.items()}


def _select_most_specific_profile(profile_by_rank: Dict[str, Dict[int, float]]) -> Dict[int, float]:
    for rank in PROFILE_RANK_PRIORITY:
        bucket = profile_by_rank.get(rank)
        if bucket:
            return dict(bucket)
    for _rank, bucket in profile_by_rank.items():
        if bucket:
            return dict(bucket)
    return {}


def _taxid_depth(
    taxid: int,
    taxonomy: Dict[int, Tuple[int, str]],
    cache: Dict[int, int],
) -> int | None:
    if taxid in cache:
        return cache[taxid]
    info = taxonomy.get(taxid)
    if info is None:
        return None
    parent, _rank = info
    if parent == taxid or parent not in taxonomy:
        cache[taxid] = 0
        return 0
    parent_depth = _taxid_depth(parent, taxonomy, cache)
    if parent_depth is None:
        return None
    cache[taxid] = parent_depth + 1
    return cache[taxid]


def _accumulate_tree_masses(
    profile: Dict[int, float],
    taxonomy: Dict[int, Tuple[int, str]],
) -> Dict[int, float]:
    masses: Dict[int, float] = {}
    for taxid, value in _normalize_profile(profile).items():
        current = taxid
        seen = set()
        while current not in seen:
            info = taxonomy.get(current)
            if info is None:
                break
            masses[current] = masses.get(current, 0.0) + value
            seen.add(current)
            parent, _rank = info
            if parent == current:
                break
            current = parent
    return masses


def compute_weighted_unifrac(
    truth: Dict[int, float],
    preds: Dict[int, float],
    taxonomy: Dict[int, Tuple[int, str]],
) -> float:
    truth_mass = _accumulate_tree_masses(truth, taxonomy)
    pred_mass = _accumulate_tree_masses(preds, taxonomy)
    if not truth_mass and not pred_mass:
        return 0.0

    depth_cache: Dict[int, int] = {}
    total = 0.0
    for taxid in set(truth_mass) | set(pred_mass):
        depth = _taxid_depth(taxid, taxonomy, depth_cache)
        if depth is None or depth <= 0:
            continue
        branch_length = 1.0 / depth
        total += branch_length * abs(truth_mass.get(taxid, 0.0) - pred_mass.get(taxid, 0.0))
    return total


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
        if _prediction_for_read(preds, read_id, truth) is not None:
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
            pred_taxid = _prediction_for_read(preds, read_id, truth)
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
        pred_taxid = _prediction_for_read(preds, read_id, truth)
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
            pred_taxid = _prediction_for_read(preds, read_id, truth)
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


def compute_opal_profile_metrics(
    truth: Dict[str, Dict[TaxKey, float]],
    preds: Dict[str, Dict[TaxKey, float]],
    ranks: Iterable[str],
):
    metrics: Dict[str, float] = {}
    for rank in ranks:
        truth_rank = truth.get(rank)
        pred_rank = preds.get(rank)
        if not truth_rank and not pred_rank:
            continue
        truth_vec = {k: float(v) for k, v in (truth_rank or {}).items() if isinstance(k, int) and v > 0}
        pred_vec = {k: float(v) for k, v in (pred_rank or {}).items() if isinstance(k, int) and v > 0}
        truth_rel = _normalize_profile(truth_vec)
        pred_rel = _normalize_profile(pred_vec)
        keys = set(truth_rel) | set(pred_rel)
        metrics[f"l1_norm_{rank}"] = sum(abs(pred_rel.get(k, 0.0) - truth_rel.get(k, 0.0)) for k in keys)

        tp = len(set(truth_rel) & set(pred_rel))
        fp = len(set(pred_rel) - set(truth_rel))
        fn = len(set(truth_rel) - set(pred_rel))
        purity = tp / (tp + fp) if (tp + fp) > 0 else 0.0
        completeness = tp / (tp + fn) if (tp + fn) > 0 else 0.0
        metrics[f"purity_{rank}"] = purity
        metrics[f"completeness_{rank}"] = completeness
    return metrics


def map_named_profile_to_taxids(
    profile: Dict[str, float],
    name_to_taxid: Dict[str, int],
    syn_to_sci: Dict[str, str],
    sci_names: set[str],
) -> tuple[Dict[int, float], int, float]:
    mapped: Dict[int, float] = {}
    unmapped_count = 0
    unmapped_mass = 0.0
    for name, value in profile.items():
        taxid = taxid_for_name(name, name_to_taxid, syn_to_sci, sci_names)
        if taxid is None:
            unmapped_count += 1
            unmapped_mass += value
            continue
        mapped[taxid] = mapped.get(taxid, 0.0) + value
    return mapped, unmapped_count, unmapped_mass


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
        sample_ids = [str(sid) for sid in (dataset.get("sample_ids") or [])]
        if not sample_ids:
            source = dataset.get("source") or {}
            sample_ids = [str(sid) for sid in (source.get("sample_ids") or [])]
        for read_path in reads:
            match = re.search(r"sample_(\d+)", str(read_path))
            if match:
                sample_ids.append(match.group(1))
        if sample_ids:
            # Preserve read-list order while de-duplicating.
            seen: set[str] = set()
            ordered = []
            for sid in sample_ids:
                if sid in seen:
                    continue
                seen.add(sid)
                ordered.append(sid)

            paths: list[Path] = []
            for sid in ordered:
                candidates = list(truth_dir.rglob(f"*sample_{sid}_*gsa_mapping.tsv"))
                if not candidates:
                    candidates = list(truth_dir.rglob(f"*sample_{sid}_*gsa_mapping.tsv.gz"))
                if not candidates:
                    candidates = list(truth_dir.glob(f"*sample_{sid}/**/reads_mapping.tsv"))
                if not candidates:
                    candidates = list(truth_dir.glob(f"*sample_{sid}/**/reads_mapping.tsv.gz"))
                if candidates:
                    paths.extend(sorted(candidates))
            if paths:
                return paths

    if reads:
        paths = []
        for read in reads:
            read_path = Path(read)
            if "fasta" not in read_path.parts:
                continue
            mapping_path = Path(str(read_path).replace("/fasta/", "/mapping/"))
            mapping_path = Path(str(mapping_path).replace("_anonymous_gsa.fasta", "_gsa_mapping.tsv"))
            if mapping_path.exists():
                paths.append(mapping_path)
                continue
            gz = Path(str(mapping_path) + ".gz")
            if gz.exists():
                paths.append(gz)
        if paths:
            return paths
    return []


def evaluate_with_truth(
    exp: dict,
    dataset: dict,
    outputs: dict,
    *,
    include_per_read: bool = True,
    include_profile: bool = True,
) -> Dict[str, float]:
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
        sci_for_taxid, syn_to_sci, sci_names = build_name_maps(names_path, taxonomy)
        name_to_taxid = build_name_taxid_maps(sci_for_taxid)
    else:
        # Fallback: use DB `.tax` names directly (no synonym collapse).
        name_to_taxid = {name: taxid for (rank, name), taxid in _name_to_taxid.items() if rank == "species"}

    metrics: Dict[str, float] = {}

    mapping_paths = _resolve_mapping_paths(dataset)
    truth_reads: Dict[str, int] = {}
    truth_abundance: Dict[int, int] = {}
    contig_weight: Dict[str, int] = {}
    if mapping_paths:
        truth_format = str(dataset.get("truth_map_format") or dataset.get("truth_format") or "cami").lower()
        if truth_format in {"species_label", "species-label", "prjna637878_species_label"}:
            (
                truth_reads,
                truth_abundance,
                contig_weight,
                unmapped_rows,
                unmapped_names,
            ) = load_species_label_mapping(
                mapping_paths,
                name_to_taxid,
                syn_to_sci,
                sci_names,
            )
            metrics["truth_map_species_label_mapped_rows"] = len(truth_reads)
            metrics["truth_map_species_label_unmapped_rows"] = unmapped_rows
            metrics["truth_map_species_label_unmapped_names"] = unmapped_names
        elif truth_format == "cami":
            truth_reads, truth_abundance, contig_weight = load_cami_mapping(mapping_paths)
        else:
            raise ValueError(f"unsupported truth_map_format: {truth_format}")

    classify_path = outputs.get("classify_tsv")
    classify_one = outputs.get("classify_one")
    preds: Dict[str, int | None] | None = None
    preds_loaded = False

    def _load_preds() -> Dict[str, int | None] | None:
        nonlocal preds, preds_loaded
        if preds_loaded:
            return preds
        preds_loaded = True
        if classify_path:
            path = Path(classify_path)
            if path.exists():
                preds = parse_classify_tsv(path)
            else:
                metrics["classify_tsv_missing"] = 1
            return preds
        if classify_one:
            path = Path(classify_one)
            if path.exists():
                preds = parse_ganon_one(path, file_to_taxid)
            else:
                metrics["classify_one_missing"] = 1
        return preds

    if include_per_read and truth_reads:
        preds = _load_preds()
    if include_per_read and preds is not None and truth_reads:
        desc_metrics = compute_per_read_metrics(truth_reads, preds, taxonomy, ranks, covered_by_rank)
        exact_metrics = compute_per_read_metrics_exact(truth_reads, preds, taxonomy, ranks, covered_by_rank)
        metrics.update(desc_metrics)
        for key, value in exact_metrics.items():
            metrics[f"exact_{key}"] = value

    truth_by_rank: Dict[str, Dict[TaxKey, float]] = {r: {} for r in ranks}
    truth_taxid_profile: Dict[int, float] = {}
    if include_profile:
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
                    truth_taxid_profile, unmapped, unmapped_mass = map_named_profile_to_taxids(
                        profile,
                        name_to_taxid,
                        syn_to_sci,
                        sci_names,
                    )
                    truth_by_rank, _truth_by_rank_full = map_taxid_profile_to_rank(
                        truth_taxid_profile,
                        taxonomy,
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
            truth_taxid_profile = {k: float(v) for k, v in truth_abundance.items()}
            truth_by_rank, _truth_by_rank_full = map_taxid_profile_to_rank(
                truth_taxid_profile,
                taxonomy,
                ranks,
                covered_by_rank,
            )

    pred_by_rank: Dict[str, Dict[TaxKey, float]] | None = None
    pred_taxid_profile: Dict[int, float] = {}
    pred_source = None
    explicit_profile_seen = False
    if outputs.get("report_abundance_tre"):
        pred_source = Path(outputs["report_abundance_tre"])
    elif outputs.get("report_reads_tre"):
        pred_source = Path(outputs["report_reads_tre"])
    if pred_source and pred_source.exists():
        explicit_profile_seen = True
        pred_taxid_profile = _select_most_specific_profile(parse_tre_counts(pred_source))
        if pred_taxid_profile:
            _pred_by_rank_mapped, pred_by_rank = map_taxid_profile_to_rank(
                pred_taxid_profile,
                taxonomy,
                ranks,
                covered_by_rank,
            )
    else:
        profile_path_str = outputs.get("profile_tsv") or outputs.get("sylph_profile_tsv")
        if profile_path_str:
            profile_path = Path(profile_path_str)
            if profile_path.exists():
                explicit_profile_seen = True
                pred_taxid_profile, _unmapped_mass, _unmapped_count = parse_sylph_profile(
                    profile_path, file_to_taxid
                )
                if pred_taxid_profile:
                    _pred_by_rank_mapped, pred_by_rank = map_taxid_profile_to_rank(
                        pred_taxid_profile,
                        taxonomy,
                        ranks,
                        covered_by_rank,
                    )
        else:
            cami_profile_path_str = outputs.get("cami_profile_tsv") or outputs.get("taxor_profile_tsv")
            if cami_profile_path_str:
                cami_path = Path(cami_profile_path_str)
                if cami_path.exists():
                    explicit_profile_seen = True
                    cami_by_rank = parse_cami_profile(cami_path)
                    if cami_by_rank:
                        pred_taxid_profile = _select_most_specific_profile(cami_by_rank)
                        if pred_taxid_profile:
                            _pred_by_rank_mapped, pred_by_rank = map_taxid_profile_to_rank(
                                pred_taxid_profile,
                                taxonomy,
                                ranks,
                                covered_by_rank,
                            )
            else:
                chimera_profile_path_str = outputs.get("chimera_profile_tsv") or outputs.get("chimera_profile")
                if chimera_profile_path_str:
                    profile_path = Path(chimera_profile_path_str)
                    if not profile_path.exists():
                        metrics["chimera_profile_missing"] = 1
                    else:
                        explicit_profile_seen = True
                        pred_profile = parse_truth_profile(profile_path)
                        if pred_profile:
                            pred_taxid_profile, _unmapped, _unmapped_mass = map_named_profile_to_taxids(
                                pred_profile,
                                name_to_taxid,
                                syn_to_sci,
                                sci_names,
                            )
                            if pred_taxid_profile:
                                _pred_by_rank_mapped, pred_by_rank = map_taxid_profile_to_rank(
                                pred_taxid_profile,
                                taxonomy,
                                ranks,
                                covered_by_rank,
                            )

    if include_profile and explicit_profile_seen and pred_by_rank is None:
        pred_by_rank = {rank: {} for rank in ranks}

    if include_profile and pred_by_rank is not None and any(truth_by_rank.get(rank) for rank in ranks):
        metrics.update(compute_opal_profile_metrics(truth_by_rank, pred_by_rank, ranks))
        metrics["weighted_unifrac"] = compute_weighted_unifrac(
            truth_taxid_profile,
            pred_taxid_profile,
            taxonomy,
        )
        metrics["profile_metric_version"] = METRIC_VERSION

    if metrics:
        metrics.setdefault("metric_version", METRIC_VERSION)
    return metrics
