"""
GeneFamilyFlow — Multi-species gene family analysis pipeline.

Usage:
    snakemake --configfile config/default_config.yaml -n -p           # dry-run
    snakemake --configfile config/default_config.yaml -j 8 --use-conda # full run
    snakemake --configfile config/default_config.yaml --until step05   # partial run

Per-step rules live under rules/*.smk and are included below.
"""

configfile: "config/default_config.yaml"

import itertools
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Config schema + semantic validation (fails fast before DAG resolution)
# ---------------------------------------------------------------------------
sys.path.insert(0, str(Path(workflow.basedir) / "scripts"))
from validate_config import validate_config, ConfigValidationError  # noqa: E402

try:
    validate_config(dict(config), Path(workflow.basedir) / "config" / "schema.json")
except ConfigValidationError as _err:
    sys.exit(str(_err))

# ---------------------------------------------------------------------------
# Globals from config
# ---------------------------------------------------------------------------
SPECIES       = config["species"]
TARGET        = config["target_species"]
DB_DIR        = config["database_dir"]
OUT_DIR       = config["output_dir"]
WORK_DIR      = config["work_dir"]
LOG_DIR       = config.get("log_dir", "logs")
MAX_THREADS   = config.get("max_threads", 8)

# Species-pair combinations for JCVI synteny
SPECIES_PAIRS = list(itertools.combinations(SPECIES, 2))

# Multi-domain support
DOMAINS     = config["step2"]["domains"]
DOMAIN_IDS  = [d["pfam_id"] for d in DOMAINS]


def get_domain_config(pfam_id):
    """Look up per-domain config by Pfam ID."""
    for d in DOMAINS:
        if d["pfam_id"] == pfam_id:
            return d
    raise ValueError(f"Domain {pfam_id} not found in config step2.domains")


def get_seed_references(pfam_id):
    """Return the list of reference species keys for a domain's BLAST seed.

    If the domain entry omits `seed_references`, fall back to every species
    declared in step3.reference_species.
    """
    domain_cfg = get_domain_config(pfam_id)
    all_refs = config["step3"]["reference_species"]
    refs = domain_cfg.get("seed_references")
    if refs is None:
        return list(all_refs.keys())
    missing = [r for r in refs if r not in all_refs]
    if missing:
        raise ValueError(
            f"Domain {pfam_id} references unknown species {missing}; "
            f"add them to step3.reference_species."
        )
    return refs


# Species ID prefix map (used when assigning gene IDs to species in steps 5/6/8).
# Falls back to the first two characters of the species name when not set.
_SPECIES_PREFIX_OVERRIDES = config.get("species_id_prefix", {}) or {}


def _species_prefix(sp):
    return _SPECIES_PREFIX_OVERRIDES.get(sp, sp[:2])


def build_species_map():
    return ",".join(f"{_species_prefix(sp)}={sp}" for sp in SPECIES)


# Tree tool
TREE_TOOL = config["step6"].get("tree_tool", "iqtree")

# Whether to run pfam_scan.pl for final Pfam verification in step 4.
# When false, step 4 trusts the per-domain hmmsearch results from step 2 and
# merges them according to `step4.domain_combine` directly, skipping pfam_scan.pl.
USE_PFAM_SCAN = config["step4"].get("use_pfam_scan", True)

# Step-level toggles
STEP8_PRECOMPUTED = config["step8"].get("precomputed", True)
STEP8_KAKS_FILE = (
    config["step8"]["kaks_file"] if STEP8_PRECOMPUTED
    else f"{WORK_DIR}/08_jcvi/kaks.tab.xls"
)

STEP9_PRECOMPUTED = config["step9"].get("precomputed", True)


def _step09_file(wc, ext):
    """Resolve MCScanX output path (precomputed or computed)."""
    if STEP9_PRECOMPUTED:
        return f"example/9.mcscanx/{wc.species}.{ext}"
    return f"{WORK_DIR}/09_mcscanx/{wc.species}.{ext}"


STEP10_COMPUTE_PROMOTER = config["step10"].get("compute_promoter", False)


def family_ids_for(sp):
    """Return the per-species gene family member ID file for ``sp``.

    Priority:
      1. ``step9.species_config[sp].npf_id_file`` — explicit user override
         (historical NPF fixture files live here).
      2. ``{WORK_DIR}/04_identification/{sp}.family.ID`` — derived at runtime
         by ``rule step04_split_by_species`` from the step4 final family set
         intersected with the species' longest-transcript proteome.

    Consumed by step09 circos, step10 promoter, and step13 GO/KEGG so all
    downstream rules share one canonical per-species membership list.
    """
    sp_cfg = (
        (config.get("step9", {}) or {}).get("species_config", {}) or {}
    ).get(sp, {}) or {}
    override = sp_cfg.get("npf_id_file")
    if override:
        return override
    return f"{WORK_DIR}/04_identification/{sp}.family.ID"

# Optional steps
STEP13_ENABLED = config.get("step13_go_kegg", {}).get("enabled", False)
STEP14_ENABLED = config.get("step14_qrtpcr", {}).get("enabled", False)

# Step 13 config routing (only referenced inside rules/step13_go_kegg.smk)
if STEP13_ENABLED:
    _STEP13_CFG = config["step13_go_kegg"]
    STEP13_SPECIES = _STEP13_CFG.get("species", []) or []
    STEP13_RUN_EGGNOG = _STEP13_CFG.get("run_eggnog", False)
    STEP13_PRECOMPUTED_DIR = _STEP13_CFG.get("precomputed_dir", "").rstrip("/")

    def _step13_annotations(wc):
        """Route per-species emapper annotation source based on run_eggnog."""
        if STEP13_RUN_EGGNOG:
            return f"{WORK_DIR}/13_go_kegg/eggnog/{wc.species}.emapper.annotations"
        return f"{STEP13_PRECOMPUTED_DIR}/{wc.species}.emapper.annotations"


# ---------------------------------------------------------------------------
# rule all — define final targets
# ---------------------------------------------------------------------------
ALL_TARGETS = [
    # Step 4 (final gene family set)
    f"{OUT_DIR}/04_identification/identify.ID.clean.fa",
    # Step 5
    f"{OUT_DIR}/05_genefamily_info/Gene_Information.xlsx",
    # Step 6
    f"{OUT_DIR}/06_tree/phylogenetic_tree_circular.pdf",
    f"{OUT_DIR}/06_tree/phylogenetic_tree_rectangular.pdf",
    # Step 7
    f"{OUT_DIR}/07_motif/meme_location.txt",
    f"{OUT_DIR}/07_motif/Tree_Domain_Motif_GeneStructure.pdf",
    # Step 8
    f"{OUT_DIR}/08_jcvi_kaks/08_jcvi_kaks.pdf",
    # Step 9
    *expand(f"{OUT_DIR}/09_circos/circos_{{sp}}.pdf", sp=SPECIES),
    # Step 10
    f"{OUT_DIR}/10_promoter/promoter_elements.pdf",
    # Step 11
    f"{OUT_DIR}/11_ppi/PPI_network.pdf",
]

if STEP13_ENABLED:
    ALL_TARGETS.extend(
        expand(f"{OUT_DIR}/13_go_kegg/{{sp}}_enrichment.pdf", sp=STEP13_SPECIES)
    )

if STEP14_ENABLED:
    ALL_TARGETS.append(f"{OUT_DIR}/14_qrt_pcr/qRT_PCR.pdf")


rule all:
    input:
        ALL_TARGETS,


# ---------------------------------------------------------------------------
# Per-step rules
# ---------------------------------------------------------------------------
include: "rules/step01_database.smk"
include: "rules/step02_hmmsearch.smk"
include: "rules/step03_blast.smk"
include: "rules/step04_identification.smk"
include: "rules/step05_info.smk"
include: "rules/step06_tree.smk"
include: "rules/step07_motif.smk"
include: "rules/step08_kaks.smk"
include: "rules/step09_circos.smk"
include: "rules/step10_promoter.smk"
include: "rules/step11_ppi.smk"
include: "rules/step13_go_kegg.smk"
include: "rules/step14_qrtpcr.smk"
