from __future__ import annotations

import logging
import shutil
from pathlib import Path

from ..utils import ensure_dir


def run(config: dict, logger: logging.Logger) -> None:
    """Prepare per-species working database files.

    This initial implementation standardizes inputs into:
    work/01_database/{species}/
    """
    base_dir = ensure_dir(Path(config["work_dir"]) / "01_database")
    for sp in config["species"]:
        name = sp["name"]
        species_dir = ensure_dir(base_dir / name)
        logger.info("Preparing species database: %s", name)
        for key in ("pep_fasta", "cds_fasta", "gff3"):
            src = Path(sp[key]).expanduser()
            if not src.exists():
                raise FileNotFoundError(f"Missing {key} for {name}: {src}")
            dst = species_dir / src.name
            shutil.copyfile(src, dst)
        genome = sp.get("genome_fasta")
        if genome:
            src = Path(genome).expanduser()
            if src.exists():
                shutil.copyfile(src, species_dir / src.name)
            else:
                logger.warning("genome_fasta not found for %s: %s", name, src)
