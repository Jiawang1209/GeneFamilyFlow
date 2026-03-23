"""Pipeline step implementations for GeneFamilyFlow."""

from . import (
    step01_database,
    step02_hmmsearch,
    step03_blast,
    step04_identification,
    step05_genefamily_info,
    step06_tree,
    step07_motif,
    step08_collinearity,
    step09_mcscanx,
    step10_promoter,
    step11_ppi,
)

STEP_REGISTRY = {
    "step01_database": step01_database.run,
    "step02_hmmsearch": step02_hmmsearch.run,
    "step03_blast": step03_blast.run,
    "step04_identification": step04_identification.run,
    "step05_genefamily_info": step05_genefamily_info.run,
    "step06_tree": step06_tree.run,
    "step07_motif": step07_motif.run,
    "step08_collinearity": step08_collinearity.run,
    "step09_mcscanx": step09_mcscanx.run,
    "step10_promoter": step10_promoter.run,
    "step11_ppi": step11_ppi.run,
}
