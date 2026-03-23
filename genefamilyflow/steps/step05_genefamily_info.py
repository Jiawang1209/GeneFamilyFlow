from __future__ import annotations

import logging
from pathlib import Path

from ..utils import ensure_dir, run_command, write_lines


def run(config: dict, logger: logging.Logger) -> None:
    step_dir = ensure_dir(Path(config["work_dir"]) / "05_genefamily_info")
    id_file = Path(config["work_dir"]) / "04_identification" / "identify.ID"
    merged_pep = Path(config["work_dir"]) / "02_hmmsearch" / "species_all.pep.fasta"
    identify_fa = step_dir / "identify.ID.fa"
    run_command(
        [config["tools"]["seqkit"], "grep", "-f", str(id_file), str(merged_pep), "-o", str(identify_fa)],
        logger,
    )

    all_gff = step_dir / "species_all.gff3"
    _concat_gff_files(config, all_gff)
    bed_path = step_dir / "species_all.bed"
    _gff_gene_to_bed(all_gff, bed_path)
    logger.info("Generated bed file: %s", bed_path)

    run_command(
        [
            config["tools"]["Rscript"],
            "R/05_genefamily_info.R",
            "--input_fasta",
            str(identify_fa),
            "--input_bed",
            str(bed_path),
            "--outdir",
            str(step_dir),
        ],
        logger,
    )


def _concat_gff_files(config: dict, out_path: Path) -> None:
    db_dir = Path(config["work_dir"]) / "01_database"
    with out_path.open("w", encoding="utf-8") as fw:
        for sp in config["species"]:
            gff_file = list((db_dir / sp["name"]).glob("*.gff*"))[0]
            with gff_file.open("r", encoding="utf-8") as fr:
                for line in fr:
                    fw.write(line)


def _gff_gene_to_bed(gff_path: Path, out_path: Path) -> None:
    lines: list[str] = []
    with gff_path.open("r", encoding="utf-8") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "gene":
                continue
            chrom, start, end, strand, attr = cols[0], cols[3], cols[4], cols[6], cols[8]
            gene_id = attr.split(";")[0].replace("ID=", "")
            lines.append(f"{chrom}\t{start}\t{end}\t{gene_id}\t.\t{strand}")
    write_lines(out_path, lines)
