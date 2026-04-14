"""Tests for scripts/validate_config.py.

Covers:
- default_config.yaml passes (the canonical reference config must stay valid).
- Structural schema failures: missing required keys, wrong types, bad enums.
- Semantic cross-field checks: target_species, step13 species subset,
  seed_references, run_eggnog/precomputed routing, step14 expression_data.
"""

from __future__ import annotations

import copy
from pathlib import Path

import pytest
import yaml

from scripts.validate_config import ConfigValidationError, validate_config

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_CONFIG = REPO_ROOT / "config" / "default_config.yaml"

pytestmark = pytest.mark.unit


@pytest.fixture
def default_cfg() -> dict:
    with DEFAULT_CONFIG.open() as fh:
        return yaml.safe_load(fh)


def _assert_error_mentions(exc: ConfigValidationError, needle: str) -> None:
    msg = str(exc)
    assert needle in msg, f"expected error to mention {needle!r}, got:\n{msg}"


def test_default_config_passes(default_cfg: dict) -> None:
    validate_config(default_cfg)


def test_missing_required_key(default_cfg: dict) -> None:
    cfg = copy.deepcopy(default_cfg)
    del cfg["species"]
    with pytest.raises(ConfigValidationError) as exc_info:
        validate_config(cfg)
    _assert_error_mentions(exc_info.value, "species")


def test_species_must_be_non_empty(default_cfg: dict) -> None:
    cfg = copy.deepcopy(default_cfg)
    cfg["species"] = []
    with pytest.raises(ConfigValidationError) as exc_info:
        validate_config(cfg)
    _assert_error_mentions(exc_info.value, "species")


def test_target_species_not_in_species(default_cfg: dict) -> None:
    cfg = copy.deepcopy(default_cfg)
    cfg["target_species"] = "NotInList"
    with pytest.raises(ConfigValidationError) as exc_info:
        validate_config(cfg)
    _assert_error_mentions(exc_info.value, "NotInList")


def test_tree_tool_invalid_enum(default_cfg: dict) -> None:
    cfg = copy.deepcopy(default_cfg)
    cfg["step6"]["tree_tool"] = "raxml"
    with pytest.raises(ConfigValidationError) as exc_info:
        validate_config(cfg)
    _assert_error_mentions(exc_info.value, "tree_tool")


def test_pfam_id_invalid_pattern(default_cfg: dict) -> None:
    cfg = copy.deepcopy(default_cfg)
    cfg["step2"]["domains"][0]["pfam_id"] = "not_a_pfam_id"
    with pytest.raises(ConfigValidationError) as exc_info:
        validate_config(cfg)
    _assert_error_mentions(exc_info.value, "pfam_id")


def test_seed_references_unknown_species(default_cfg: dict) -> None:
    cfg = copy.deepcopy(default_cfg)
    cfg["step2"]["domains"][0]["seed_references"] = ["Athaliana", "Bogus"]
    with pytest.raises(ConfigValidationError) as exc_info:
        validate_config(cfg)
    _assert_error_mentions(exc_info.value, "Bogus")


def test_step13_rejects_unknown_species(default_cfg: dict) -> None:
    cfg = copy.deepcopy(default_cfg)
    cfg["step13_go_kegg"]["enabled"] = True
    cfg["step13_go_kegg"]["species"] = ["Athaliana", "NotASpecies"]
    with pytest.raises(ConfigValidationError) as exc_info:
        validate_config(cfg)
    _assert_error_mentions(exc_info.value, "NotASpecies")


def test_step13_run_eggnog_requires_data_dir(default_cfg: dict) -> None:
    cfg = copy.deepcopy(default_cfg)
    cfg["step13_go_kegg"]["enabled"] = True
    cfg["step13_go_kegg"]["run_eggnog"] = True
    # Leaves the placeholder /path/to/eggnog_data in place
    with pytest.raises(ConfigValidationError) as exc_info:
        validate_config(cfg)
    _assert_error_mentions(exc_info.value, "eggnog_data_dir")


def test_step13_run_eggnog_accepts_real_path(default_cfg: dict) -> None:
    cfg = copy.deepcopy(default_cfg)
    cfg["step13_go_kegg"]["enabled"] = True
    cfg["step13_go_kegg"]["run_eggnog"] = True
    cfg["step13_go_kegg"]["eggnog_data_dir"] = "/data/eggnog"
    validate_config(cfg)


def test_step14_enabled_requires_expression_data(default_cfg: dict) -> None:
    cfg = copy.deepcopy(default_cfg)
    cfg["step14_qrtpcr"]["enabled"] = True
    cfg["step14_qrtpcr"]["expression_data"] = ""
    with pytest.raises(ConfigValidationError) as exc_info:
        validate_config(cfg)
    _assert_error_mentions(exc_info.value, "expression_data")


def test_wrong_type_on_max_threads(default_cfg: dict) -> None:
    cfg = copy.deepcopy(default_cfg)
    cfg["max_threads"] = "eight"
    with pytest.raises(ConfigValidationError) as exc_info:
        validate_config(cfg)
    _assert_error_mentions(exc_info.value, "max_threads")
