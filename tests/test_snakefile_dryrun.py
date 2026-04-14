"""Smoke tests: `snakemake -n` must succeed for every supported toggle combo.

These tests guard the Snakefile DAG against regressions when rules are
re-organized, renamed, or split. They do NOT execute any rule — they only
validate that Snakemake can resolve the DAG end-to-end.

Each test builds a temporary config by patching the default config in-memory
and writes it to the pytest tmp dir, then invokes `snakemake -n` from the
repository root. A non-zero exit code indicates a broken DAG (missing inputs,
wildcard mismatches, cyclic dependencies, unknown rules, etc.).
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest
import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_CONFIG = REPO_ROOT / "config" / "default_config.yaml"

pytestmark = pytest.mark.skipif(
    shutil.which("snakemake") is None,
    reason="snakemake CLI not available on PATH",
)


EMAPPER_STUB_HEADER = (
    "## emapper-2.1.12\n"
    "#query\tseed_ortholog\tevalue\tscore\teggNOG_OGs\tmax_annot_lvl\t"
    "COG_category\tDescription\tPreferred_name\tGOs\tEC\tKEGG_ko\t"
    "KEGG_Pathway\tKEGG_Module\tKEGG_Reaction\tKEGG_rclass\tBRITE\t"
    "KEGG_TC\tCAZy\tBiGG_Reaction\tPFAMs\n"
)
EMAPPER_STUB_ROW = (
    "gene1\t-\t1e-50\t100\t-\t-\t-\t-\t-\tGO:0008150\t-\t-\t"
    "ko04075\t-\t-\t-\t-\t-\t-\t-\t-\n"
)


def _load_default_config() -> dict:
    with DEFAULT_CONFIG.open() as fh:
        return yaml.safe_load(fh)


def _write_config(tmp_path: Path, patch: dict) -> Path:
    cfg = _load_default_config()

    def _deep_merge(dst: dict, src: dict) -> None:
        for k, v in src.items():
            if isinstance(v, dict) and isinstance(dst.get(k), dict):
                _deep_merge(dst[k], v)
            else:
                dst[k] = v

    _deep_merge(cfg, patch)
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump(cfg))
    return path


def _run_dryrun(config_path: Path) -> subprocess.CompletedProcess:
    return subprocess.run(
        ["snakemake", "-n", "--configfile", str(config_path)],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        timeout=120,
    )


def _assert_dag_ok(result: subprocess.CompletedProcess) -> None:
    assert result.returncode == 0, (
        f"snakemake -n failed (exit {result.returncode})\n"
        f"STDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
    )
    # Snakemake prints "total" in the job-stats table when DAG resolved.
    assert "total" in result.stdout, (
        f"snakemake -n did not emit a job stats table:\n{result.stdout}"
    )


@pytest.mark.integration
def test_dryrun_default(tmp_path: Path) -> None:
    """Baseline: default config (step13/14 disabled) resolves."""
    cfg_path = _write_config(tmp_path, {})
    _assert_dag_ok(_run_dryrun(cfg_path))


@pytest.mark.integration
def test_dryrun_step13_precomputed(tmp_path: Path) -> None:
    """Step 13 with run_eggnog=false reads precomputed emapper annotations."""
    stub_dir = tmp_path / "eggnog_stub"
    stub_dir.mkdir()
    for sp in ("Athaliana", "Osativa"):
        (stub_dir / f"{sp}.emapper.annotations").write_text(
            EMAPPER_STUB_HEADER + EMAPPER_STUB_ROW
        )
    cfg_path = _write_config(
        tmp_path,
        {
            "step13_go_kegg": {
                "enabled": True,
                "run_eggnog": False,
                "species": ["Athaliana", "Osativa"],
                "precomputed_dir": str(stub_dir),
            }
        },
    )
    _assert_dag_ok(_run_dryrun(cfg_path))


@pytest.mark.integration
def test_dryrun_step13_run_eggnog(tmp_path: Path) -> None:
    """Step 13 with run_eggnog=true wires the emapper rule into the DAG."""
    cfg_path = _write_config(
        tmp_path,
        {
            "step13_go_kegg": {
                "enabled": True,
                "run_eggnog": True,
                "species": ["Athaliana", "Osativa"],
                "eggnog_data_dir": "/tmp/fake_eggnog_data",
            }
        },
    )
    result = _run_dryrun(cfg_path)
    _assert_dag_ok(result)
    # Confirm the per-species eggnog rule is actually scheduled.
    assert "step13_eggnog_run" in result.stdout, (
        f"step13_eggnog_run missing from DAG:\n{result.stdout}"
    )


@pytest.mark.integration
def test_dryrun_step14(tmp_path: Path) -> None:
    """Step 14 (qRT-PCR) toggles cleanly into the DAG."""
    cfg_path = _write_config(tmp_path, {"step14_qrtpcr": {"enabled": True}})
    result = _run_dryrun(cfg_path)
    _assert_dag_ok(result)
    assert "step14_qrt_pcr" in result.stdout


@pytest.mark.integration
def test_dryrun_step13_rejects_unknown_species(tmp_path: Path) -> None:
    """Config validation must fire when step13.species contains an unknown sp."""
    cfg_path = _write_config(
        tmp_path,
        {
            "step13_go_kegg": {
                "enabled": True,
                "species": ["Athaliana", "NotASpecies"],
            }
        },
    )
    result = _run_dryrun(cfg_path)
    assert result.returncode != 0, (
        "Expected dry-run to fail for unknown species in step13_go_kegg.species"
    )
    combined = result.stdout + result.stderr
    assert "NotASpecies" in combined, f"Validation error missing species name:\n{combined}"
