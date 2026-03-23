from __future__ import annotations

import logging
import re
from pathlib import Path

from ..utils import ensure_dir, run_command


def run(config: dict, logger: logging.Logger) -> None:
    step_dir = ensure_dir(Path(config["work_dir"]) / "07_motif")
    fasta = Path(config["work_dir"]) / "05_genefamily_info" / "identify.ID.fa"
    meme_out = step_dir / "meme_out"
    run_command(
        [
            config["tools"]["meme"],
            str(fasta),
            "-oc",
            str(meme_out),
            "-mod",
            "anr",
            "-protein",
            "-nmotifs",
            str(config["motif"]["nmotifs"]),
            "-minw",
            str(config["motif"]["minw"]),
            "-maxw",
            str(config["motif"]["maxw"]),
        ],
        logger,
    )
    _parse_meme_txt(meme_out / "meme.txt", step_dir)
    run_command(
        [
            config["tools"]["Rscript"],
            "R/07_domain_motif_structure.R",
            "--outdir",
            str(step_dir),
        ],
        logger,
    )


def _parse_meme_txt(meme_txt: Path, out_dir: Path) -> None:
    motif_name = ""
    with (out_dir / "meme_info.txt").open("w", encoding="utf-8") as fw1, (
        out_dir / "meme_location.txt"
    ).open("w", encoding="utf-8") as fw2:
        with meme_txt.open("r", encoding="utf-8") as fr:
            for line in fr:
                line_tmp = line.strip()
                if re.match(r"Motif.*Description", line_tmp):
                    fields = line_tmp.split()
                    motif_name = fields[2]
                    fw1.write(f"{fields[2]}\t{fields[1]}\t{len(fields[1])}\n")
                if re.match(r"^\S+\s+\d+\s+\d+.*", line_tmp) and motif_name:
                    fields = line_tmp.split()
                    start = int(fields[1])
                    end = start + len(fields[4]) - 1
                    fw2.write(f"{fields[0]}\t{start}\t{end}\t{motif_name}\n")
