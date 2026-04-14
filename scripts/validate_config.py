"""Config validator for GeneFamilyFlow.

Loads `config/schema.json` and validates a Snakemake config dict against it,
then runs cross-field semantic checks that JSON Schema cannot express cleanly
(e.g. target_species must be in species, step13 species subset of species,
step2.domains[].seed_references subset of step3.reference_species).

Designed to be called from the top of the Snakefile (after `configfile:`)
so misconfigured runs fail with a clear message *before* DAG resolution.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_SCHEMA_PATH = REPO_ROOT / "config" / "schema.json"


class ConfigValidationError(ValueError):
    """Raised when the config fails schema or semantic validation."""


@dataclass
class _Report:
    errors: list[str] = field(default_factory=list)

    def add(self, msg: str) -> None:
        self.errors.append(msg)

    def raise_if_any(self) -> None:
        if not self.errors:
            return
        bullet = "\n  - "
        raise ConfigValidationError(
            "Config validation failed:" + bullet + bullet.join(self.errors)
        )


def _load_schema(schema_path: Path) -> dict[str, Any]:
    with schema_path.open() as fh:
        return json.load(fh)


def _format_jsonschema_path(path: Any) -> str:
    parts = [str(p) for p in path]
    return ".".join(parts) if parts else "<root>"


def _run_jsonschema(config: dict[str, Any], schema: dict[str, Any], report: _Report) -> None:
    try:
        import jsonschema
    except ImportError as exc:  # pragma: no cover - env failure
        raise ConfigValidationError(
            "jsonschema is required for config validation but is not installed. "
            "Install it via `pip install jsonschema` or add it to envs/genefamily.yaml."
        ) from exc

    validator = jsonschema.Draft7Validator(schema)
    for err in sorted(validator.iter_errors(config), key=lambda e: list(e.absolute_path)):
        loc = _format_jsonschema_path(err.absolute_path)
        report.add(f"[schema] {loc}: {err.message}")


def _run_semantic_checks(config: dict[str, Any], report: _Report) -> None:
    species: list[str] = list(config.get("species") or [])
    species_set = set(species)

    target = config.get("target_species")
    if isinstance(target, str) and species and target not in species_set:
        report.add(
            f"target_species '{target}' must be one of top-level species {species}"
        )

    step2 = config.get("step2") or {}
    step3 = config.get("step3") or {}
    refs = step3.get("reference_species") or {}
    ref_keys = set(refs.keys())

    for i, domain in enumerate(step2.get("domains") or []):
        if not isinstance(domain, dict):
            continue
        seeds = domain.get("seed_references")
        if seeds is None:
            continue
        missing = [s for s in seeds if s not in ref_keys]
        if missing:
            pfam = domain.get("pfam_id", f"index {i}")
            report.add(
                f"step2.domains[{pfam}].seed_references contains species not in "
                f"step3.reference_species: {missing}"
            )

    s13 = config.get("step13_go_kegg") or {}
    if s13.get("enabled"):
        s13_species = list(s13.get("species") or [])
        if not s13_species:
            report.add("step13_go_kegg.enabled=true but step13_go_kegg.species is empty")
        unknown = [sp for sp in s13_species if sp not in species_set]
        if unknown:
            report.add(
                f"step13_go_kegg.species contains unknown species {unknown}; "
                f"must be a subset of top-level species {species}"
            )
        if s13.get("run_eggnog"):
            data_dir = s13.get("eggnog_data_dir") or ""
            if not data_dir or data_dir == "/path/to/eggnog_data":
                report.add(
                    "step13_go_kegg.run_eggnog=true but eggnog_data_dir is not set "
                    "(still the placeholder). Point it at your eggNOG data directory."
                )
        else:
            pre = s13.get("precomputed_dir") or ""
            if not pre:
                report.add(
                    "step13_go_kegg.run_eggnog=false but precomputed_dir is empty; "
                    "supply a directory containing {species}.emapper.annotations files."
                )

    s14 = config.get("step14_qrtpcr") or {}
    if s14.get("enabled"):
        if not s14.get("expression_data"):
            report.add("step14_qrtpcr.enabled=true but expression_data is empty")


def validate_config(
    config: dict[str, Any],
    schema_path: Path | str = DEFAULT_SCHEMA_PATH,
) -> None:
    """Validate a config dict in-place. Raises ConfigValidationError on failure."""
    schema = _load_schema(Path(schema_path))
    report = _Report()
    _run_jsonschema(config, schema, report)
    _run_semantic_checks(config, report)
    report.raise_if_any()


def main(argv: list[str] | None = None) -> int:
    import argparse
    import sys

    import yaml

    parser = argparse.ArgumentParser(description="Validate a GeneFamilyFlow config file.")
    parser.add_argument("config", type=Path, help="Path to config YAML")
    parser.add_argument(
        "--schema",
        type=Path,
        default=DEFAULT_SCHEMA_PATH,
        help="Path to schema JSON (default: config/schema.json)",
    )
    args = parser.parse_args(argv)

    with args.config.open() as fh:
        cfg = yaml.safe_load(fh)

    try:
        validate_config(cfg, args.schema)
    except ConfigValidationError as exc:
        print(str(exc), file=sys.stderr)
        return 1

    print(f"OK: {args.config} passes schema and semantic checks.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
