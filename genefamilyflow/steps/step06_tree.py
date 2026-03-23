from __future__ import annotations

import logging
from pathlib import Path

from ..utils import ensure_dir, run_command


def run(config: dict, logger: logging.Logger) -> None:
    step_dir = ensure_dir(Path(config["work_dir"]) / "06_tree")
    fasta = Path(config["work_dir"]) / "05_genefamily_info" / "identify.ID.fa"
    aln = step_dir / "identify.ID.aln"
    run_command([config["tools"]["muscle"], "-in", str(fasta), "-out", str(aln)], logger)
    run_command(
        [
            config["tools"]["iqtree"],
            "-s",
            str(aln),
            "-m",
            str(config["tree"]["iqtree_model"]),
            "-bb",
            str(config["tree"]["bootstrap"]),
            "-nt",
            "AUTO",
            "-redo",
        ],
        logger,
    )
    run_command(
        [
            config["tools"]["Rscript"],
            "R/06_tree.R",
            "--treefile",
            str(aln) + ".treefile",
            "--outdir",
            str(step_dir),
        ],
        logger,
    )
