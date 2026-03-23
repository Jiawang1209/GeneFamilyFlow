from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml


REQUIRED_TOP_LEVEL_KEYS = {
    "project_name",
    "output_dir",
    "work_dir",
    "log_dir",
    "species",
    "reference_species",
    "pfam",
    "blast",
    "tree",
    "motif",
    "promoter",
    "tools",
}

REQUIRED_SPECIES_KEYS = {"name", "pep_fasta", "cds_fasta", "gff3"}


@dataclass
class LoadedConfig:
    config_path: Path
    data: dict[str, Any]


def load_config(config_path: str | Path) -> LoadedConfig:
    path = Path(config_path).expanduser().resolve()
    if not path.exists():
        raise FileNotFoundError(f"Config file not found: {path}")
    with path.open("r", encoding="utf-8") as fh:
        data = yaml.safe_load(fh) or {}
    validate_config(data, path)
    return LoadedConfig(config_path=path, data=data)


def validate_config(data: dict[str, Any], path: Path | None = None) -> None:
    missing = REQUIRED_TOP_LEVEL_KEYS.difference(data.keys())
    if missing:
        where = str(path) if path else "config"
        raise ValueError(f"Missing required keys in {where}: {sorted(missing)}")

    species = data.get("species", [])
    if not isinstance(species, list) or not species:
        raise ValueError("'species' must be a non-empty list")

    for idx, item in enumerate(species, start=1):
        if not isinstance(item, dict):
            raise ValueError(f"species[{idx}] must be a mapping")
        missing_species = REQUIRED_SPECIES_KEYS.difference(item.keys())
        if missing_species:
            raise ValueError(
                f"species[{idx}] is missing required keys: {sorted(missing_species)}"
            )

    ref = data.get("reference_species", {})
    if not isinstance(ref, dict):
        raise ValueError("'reference_species' must be a mapping")
    if "name" not in ref or "family_member_ids" not in ref:
        raise ValueError(
            "'reference_species' must include 'name' and 'family_member_ids'"
        )


def as_path(config: LoadedConfig, key: str) -> Path:
    raw = config.data.get(key)
    if raw is None:
        raise KeyError(f"Missing config key: {key}")
    return Path(raw).expanduser().resolve()
