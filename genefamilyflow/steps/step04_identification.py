from __future__ import annotations

import logging
from pathlib import Path

from ..utils import ensure_dir, run_command, write_lines


def run(config: dict, logger: logging.Logger) -> None:
    step_dir = ensure_dir(Path(config["work_dir"]) / "04_identification")
    hmm_ids = _read_id_set(Path(config["work_dir"]) / "02_hmmsearch" / "2st_id")
    blast_ids = _read_id_set(Path(config["work_dir"]) / "03_blast" / "species.blast.id")
    inter_ids = sorted(hmm_ids.intersection(blast_ids))
    write_lines(step_dir / "inter.ID", inter_ids)
    write_lines(step_dir / "union.ID", sorted(hmm_ids.union(blast_ids)))

    merged_pep = Path(config["work_dir"]) / "02_hmmsearch" / "species_all.pep.fasta"
    inter_fa = step_dir / "inter.ID.fa"
    run_command(
        [config["tools"]["seqkit"], "grep", "-f", str(step_dir / "inter.ID"), str(merged_pep), "-o", str(inter_fa)],
        logger,
    )

    pfam_hmm = config["pfam"].get("hmm_file")
    if not pfam_hmm:
        logger.warning("pfam.hmm_file is not set; skipping hmmscan validation")
        write_lines(step_dir / "identify.ID", inter_ids)
        run_command(
            [config["tools"]["seqkit"], "grep", "-f", str(step_dir / "identify.ID"), str(merged_pep), "-o", str(step_dir / "identify.ID.fa")],
            logger,
        )
        return

    tbl = step_dir / "pfam.tbl"
    run_command(
        [
            config["tools"]["hmmscan"],
            "--tblout",
            str(tbl),
            "-o",
            str(step_dir / "pfam.out.txt"),
            str(Path(pfam_hmm).expanduser().resolve()),
            str(inter_fa),
        ],
        logger,
    )
    pfam_id = config["pfam"]["pfam_id"]
    retained = _filter_hmmscan_ids(tbl, pfam_id)
    write_lines(step_dir / "identify.ID", retained)
    run_command(
        [
            config["tools"]["seqkit"],
            "grep",
            "-f",
            str(step_dir / "identify.ID"),
            str(merged_pep),
            "-o",
            str(step_dir / "identify.ID.fa"),
        ],
        logger,
    )


def _read_id_set(path: Path) -> set[str]:
    if not path.exists():
        return set()
    with path.open("r", encoding="utf-8") as fh:
        return {line.strip() for line in fh if line.strip()}


def _filter_hmmscan_ids(tbl_path: Path, target_pfam: str) -> list[str]:
    out: set[str] = set()
    with tbl_path.open("r", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 3:
                continue
            if fields[1].startswith(target_pfam):
                out.add(fields[2])
    return sorted(out)
