from __future__ import annotations

import logging
from pathlib import Path

from ..utils import ensure_dir, run_command


def run(config: dict, logger: logging.Logger) -> None:
    step_dir = ensure_dir(Path(config["work_dir"]) / "08_collinearity")
    logger.info("Step08 currently expects JCVI input preparation by user data config.")
    run_command(
        [
            config["tools"]["Rscript"],
            "R/08_jcvi_kaks.R",
            "--outdir",
            str(step_dir),
        ],
        logger,
    )
