from __future__ import annotations

import logging
from pathlib import Path

from ..utils import ensure_dir, run_command


def run(config: dict, logger: logging.Logger) -> None:
    step_dir = ensure_dir(Path(config["work_dir"]) / "10_promoter")
    logger.info("Step10 promoter extraction wrapper.")
    run_command(
        [
            config["tools"]["Rscript"],
            "R/10_promoter.R",
            "--outdir",
            str(step_dir),
        ],
        logger,
    )
