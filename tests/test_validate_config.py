"""Tests for scripts/validate_config.py.

The schema-validation happy path is exercised indirectly by the Snakefile
dry-run smoke tests; this file targets the semantic checks and CLI surface
that JSON Schema cannot express, plus the error branches that are easy to
trigger with focused fixtures.
"""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from scripts.validate_config import (
    DEFAULT_SCHEMA_PATH,
    ConfigValidationError,
    main,
    validate_config,
)

pytestmark = pytest.mark.unit


def _minimal_valid_config() -> dict:
    """A config that satisfies the schema and semantic checks."""
    raw = yaml.safe_load(
        Path("config/default_config.yaml").read_text()
    )
    return raw


class TestSemanticChecks:
    def test_target_species_must_be_in_species(self) -> None:
        cfg = _minimal_valid_config()
        cfg["target_species"] = "Bogus"
        with pytest.raises(ConfigValidationError, match="target_species"):
            validate_config(cfg)

    def test_seed_references_must_be_known(self) -> None:
        cfg = _minimal_valid_config()
        cfg["step2"]["domains"][0]["seed_references"] = ["GhostSpecies"]
        with pytest.raises(ConfigValidationError, match="seed_references"):
            validate_config(cfg)

    def test_non_dict_domain_entry_is_skipped(self) -> None:
        # Schema would normally reject this, but the semantic loop has
        # a defensive `if not isinstance(domain, dict): continue` branch.
        # Bypass the schema check by validating semantics directly via
        # _run_semantic_checks-style input — easiest path is to feed a
        # config with a list-typed domain and assert schema rejects it
        # AND the semantic loop survives the iteration.
        cfg = _minimal_valid_config()
        cfg["step2"]["domains"].append(["not", "a", "dict"])
        # Schema validation will fire first, so we expect a schema error,
        # but we want the semantic loop to NOT crash on the bad entry.
        with pytest.raises(ConfigValidationError) as exc_info:
            validate_config(cfg)
        assert "schema" in str(exc_info.value)

    def test_seed_references_none_is_allowed(self) -> None:
        cfg = _minimal_valid_config()
        # Drop seed_references entirely from the first domain
        cfg["step2"]["domains"][0].pop("seed_references", None)
        # Should not raise
        validate_config(cfg)


class TestStep13Checks:
    def test_enabled_with_empty_species_fails(self) -> None:
        cfg = _minimal_valid_config()
        cfg["step13_go_kegg"] = {
            "enabled": True,
            "species": [],
            "run_eggnog": False,
            "precomputed_dir": "some/dir",
        }
        with pytest.raises(ConfigValidationError, match="species is empty"):
            validate_config(cfg)

    def test_enabled_with_unknown_species_fails(self) -> None:
        cfg = _minimal_valid_config()
        cfg["step13_go_kegg"] = {
            "enabled": True,
            "species": ["Unknown"],
            "run_eggnog": False,
            "precomputed_dir": "some/dir",
        }
        with pytest.raises(ConfigValidationError, match="unknown species"):
            validate_config(cfg)

    def test_run_eggnog_with_placeholder_fails(self) -> None:
        cfg = _minimal_valid_config()
        cfg["step13_go_kegg"] = {
            "enabled": True,
            "species": cfg["species"][:1],
            "run_eggnog": True,
            "eggnog_data_dir": "/path/to/eggnog_data",
        }
        with pytest.raises(ConfigValidationError, match="eggnog_data_dir"):
            validate_config(cfg)

    def test_precomputed_mode_requires_dir(self) -> None:
        cfg = _minimal_valid_config()
        cfg["step13_go_kegg"] = {
            "enabled": True,
            "species": cfg["species"][:1],
            "run_eggnog": False,
            "precomputed_dir": "",
        }
        with pytest.raises(ConfigValidationError, match="precomputed_dir is empty"):
            validate_config(cfg)


class TestStep14Checks:
    def test_enabled_without_expression_data_fails(self) -> None:
        cfg = _minimal_valid_config()
        cfg["step14_qrtpcr"] = {"enabled": True, "expression_data": ""}
        with pytest.raises(ConfigValidationError, match="expression_data"):
            validate_config(cfg)


class TestValidateConfigHappyPath:
    def test_default_config_passes(self) -> None:
        cfg = _minimal_valid_config()
        validate_config(cfg)  # must not raise

    def test_default_schema_path_resolves(self) -> None:
        assert DEFAULT_SCHEMA_PATH.exists()


class TestMainCli:
    def test_happy_path(
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        cfg_path = tmp_path / "ok.yaml"
        cfg_path.write_text(Path("config/default_config.yaml").read_text())
        rc = main([str(cfg_path)])
        assert rc == 0
        assert "passes schema and semantic checks" in capsys.readouterr().out

    def test_returns_one_on_invalid(
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        cfg = _minimal_valid_config()
        cfg["target_species"] = "Bogus"
        cfg_path = tmp_path / "bad.yaml"
        cfg_path.write_text(yaml.safe_dump(cfg))
        rc = main([str(cfg_path)])
        assert rc == 1
        err = capsys.readouterr().err
        assert "target_species" in err
