from __future__ import annotations

import logging
from pathlib import Path

from ..utils import ensure_dir, run_command


def run(config: dict, logger: logging.Logger) -> None:
    step_dir = ensure_dir(Path(config["work_dir"]) / "11_ppi")
    logger.info("Step11 PPI network wrapper.")
    run_command(
        [
            config["tools"]["Rscript"],
            "R/11_ppi.R",
            "--outdir",
            str(step_dir),
        ],
        logger,
    )
