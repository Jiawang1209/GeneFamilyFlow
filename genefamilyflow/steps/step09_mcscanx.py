from __future__ import annotations

import logging
from pathlib import Path

from ..utils import ensure_dir, run_command


def run(config: dict, logger: logging.Logger) -> None:
    step_dir = ensure_dir(Path(config["work_dir"]) / "09_mcscanx")
    logger.info("Step09 scaffolds MCScanX output parsing and visualization.")
    run_command(
        [
            config["tools"]["Rscript"],
            "R/09_circos.R",
            "--outdir",
            str(step_dir),
        ],
        logger,
    )
