from __future__ import annotations

import logging
from pathlib import Path

from ..utils import ensure_dir, run_command, write_lines


def run(config: dict, logger: logging.Logger) -> None:
    step_dir = ensure_dir(Path(config["work_dir"]) / "03_blast")
    db_dir = Path(config["work_dir"]) / "01_database"
    blast_cfg = config["blast"]
    tools = config["tools"]
    ref_cfg = config["reference_species"]

    ref_id_file = Path(ref_cfg["family_member_ids"]).expanduser().resolve()
    if not ref_id_file.exists():
        raise FileNotFoundError(f"Reference family_member_ids not found: {ref_id_file}")

    ref_species_dir = db_dir / ref_cfg["name"]
    pep_fasta = list(ref_species_dir.glob("*.pep.fasta"))
    if not pep_fasta:
        raise FileNotFoundError(f"Reference species pep file missing in: {ref_species_dir}")
    ref_pep = pep_fasta[0]

    ref_fasta = step_dir / "reference_family.pep.fasta"
    run_command(
        [tools["seqkit"], "grep", "-f", str(ref_id_file), str(ref_pep), "-o", str(ref_fasta)],
        logger,
    )
    run_command(
        [tools["makeblastdb"], "-in", str(ref_fasta), "-dbtype", "prot", "-parse_seqids"],
        logger,
    )

    all_hits: set[str] = set()
    for sp in config["species"]:
        name = sp["name"]
        species_pep = list((db_dir / name).glob("*.pep.fasta"))[0]
        out_file = step_dir / f"{name}.blast"
        run_command(
            [
                tools["blastp"],
                "-query",
                str(species_pep),
                "-db",
                str(ref_fasta),
                "-outfmt",
                "6 std qlen slen",
                "-out",
                str(out_file),
                "-evalue",
                str(blast_cfg["evalue"]),
                "-num_threads",
                str(blast_cfg["num_threads"]),
                "-num_alignments",
                str(blast_cfg["num_alignments"]),
            ],
            logger,
        )
        with out_file.open("r", encoding="utf-8") as fh:
            for line in fh:
                cols = line.strip().split("\t")
                if cols:
                    all_hits.add(cols[0])

    write_lines(step_dir / "species.blast.id", sorted(all_hits))
