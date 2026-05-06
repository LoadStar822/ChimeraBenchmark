from __future__ import annotations

from pathlib import Path
from typing import Any, Dict

import yaml


DATASET_COLLECTION_KEYS = {"samples", "catalog"}


def load_yaml(path: Path) -> Dict[str, Any]:
    with path.open("r", encoding="utf-8") as f:
        data = yaml.safe_load(f) or {}
    if not isinstance(data, dict):
        raise ValueError(f"Invalid YAML object in {path}")
    return data


def load_yaml_dir(dir_path: Path) -> Dict[str, Dict[str, Any]]:
    out: Dict[str, Dict[str, Any]] = {}
    for path in sorted(dir_path.glob("*.y*ml")):
        data = load_yaml(path)
        name = data.get("name") or path.stem
        out[name] = data
    return out


def expand_dataset_config(name: str, data: Dict[str, Any]) -> list[Dict[str, Any]]:
    collection_name = str(data.get("name") or name)
    samples = data.get("samples")
    if samples is None:
        item = dict(data)
        item["name"] = item.get("name", collection_name)
        return [item]
    if not isinstance(samples, list) or not samples:
        raise ValueError(f"dataset samples must be a non-empty list: {collection_name}")

    inherited = {k: v for k, v in data.items() if k not in DATASET_COLLECTION_KEYS}
    expanded: list[Dict[str, Any]] = []
    for idx, sample in enumerate(samples, start=1):
        if not isinstance(sample, dict):
            raise ValueError(f"dataset sample must be an object: {collection_name}[{idx}]")
        sample_id = sample.get("sample_id") or sample.get("id")
        if not sample_id:
            raise ValueError(f"dataset sample is missing sample_id: {collection_name}[{idx}]")
        item = dict(inherited)
        item.update(sample)
        item["sample_id"] = str(sample_id)
        item["dataset_collection"] = collection_name
        item["name"] = str(sample.get("name") or f"{collection_name}.{sample_id}")
        expanded.append(item)
    return expanded


def expand_datasets(datasets: Dict[str, Dict[str, Any]]) -> Dict[str, Dict[str, Any]]:
    expanded: Dict[str, Dict[str, Any]] = {}
    for name, data in datasets.items():
        for item in expand_dataset_config(name, data):
            item_name = str(item["name"])
            if item_name in expanded:
                raise ValueError(f"duplicate expanded dataset name: {item_name}")
            expanded[item_name] = item
    return expanded
