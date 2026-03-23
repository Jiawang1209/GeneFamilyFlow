from __future__ import annotations

import logging
from pathlib import Path

from ..utils import ensure_dir, run_command, write_lines


def run(config: dict, logger: logging.Logger) -> None:
    step_dir = ensure_dir(Path(config["work_dir"]) / "02_hmmsearch")
    db_dir = Path(config["work_dir"]) / "01_database"
    pfam = config["pfam"]
    tools = config["tools"]
    hmm_file = pfam.get("hmm_file")
    if not hmm_file:
        raise ValueError("pfam.hmm_file must be provided for step02_hmmsearch")
    hmm_path = Path(hmm_file).expanduser().resolve()
    if not hmm_path.exists():
        raise FileNotFoundError(f"HMM file not found: {hmm_path}")

    round1_id_files: list[Path] = []
    for sp in config["species"]:
        name = sp["name"]
        species_dir = db_dir / name
        pep_files = list(species_dir.glob("*.pep.fasta"))
        if not pep_files:
            raise FileNotFoundError(f"No pep fasta found in {species_dir}")
        pep_fasta = pep_files[0]
        domtbl = step_dir / f"{name}.round1.domtblout"
        out = step_dir / f"{name}.round1.hmmout"
        run_command(
            [
                tools["hmmsearch"],
                "--cut_tc",
                "--domtblout",
                str(domtbl),
                "-o",
                str(out),
                str(hmm_path),
                str(pep_fasta),
            ],
            logger,
        )
        id_file = step_dir / f"{name}.round1.ids"
        _extract_hmm_ids(domtbl, id_file, float(pfam["hmmsearch_evalue"]))
        round1_id_files.append(id_file)

    merged_round1 = step_dir / "round1_ids.txt"
    write_lines(merged_round1, _merge_id_files(round1_id_files))

    merged_fasta = step_dir / "round1_ids.fa"
    run_command(
        [
            tools["seqkit"],
            "grep",
            "-f",
            str(merged_round1),
            _collect_merged_pep(config),
            "-o",
            str(merged_fasta),
        ],
        logger,
    )

    aln_out = step_dir / "round1_ids.aln"
    run_command(
        [
            tools["muscle"],
            "-in",
            str(merged_fasta),
            "-out",
            str(aln_out),
        ],
        logger,
    )
    custom_hmm = step_dir / f"custom_{pfam['pfam_id']}.hmm"
    run_command([tools["hmmbuild"], str(custom_hmm), str(aln_out)], logger)

    round2_id_files: list[Path] = []
    for sp in config["species"]:
        name = sp["name"]
        species_dir = db_dir / name
        pep_fasta = list(species_dir.glob("*.pep.fasta"))[0]
        domtbl = step_dir / f"{name}.round2.domtblout"
        out = step_dir / f"{name}.round2.hmmout"
        run_command(
            [
                tools["hmmsearch"],
                "--domtblout",
                str(domtbl),
                "-o",
                str(out),
                str(custom_hmm),
                str(pep_fasta),
            ],
            logger,
        )
        id_file = step_dir / f"{name}.round2.ids"
        _extract_hmm_ids(domtbl, id_file, float(pfam["hmmsearch_evalue"]))
        round2_id_files.append(id_file)
    write_lines(step_dir / "2st_id", _merge_id_files(round2_id_files))


def _extract_hmm_ids(domtblout: Path, out_path: Path, evalue_cutoff: float) -> None:
    hits: list[str] = []
    with domtblout.open("r", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 7:
                continue
            try:
                evalue = float(fields[6])
            except ValueError:
                continue
            if evalue < evalue_cutoff:
                hits.append(fields[0])
    write_lines(out_path, sorted(set(hits)))


def _merge_id_files(paths: list[Path]) -> list[str]:
    values: set[str] = set()
    for p in paths:
        if not p.exists():
            continue
        with p.open("r", encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if line:
                    values.add(line)
    return sorted(values)


def _collect_merged_pep(config: dict) -> str:
    out = Path(config["work_dir"]) / "02_hmmsearch" / "species_all.pep.fasta"
    db_dir = Path(config["work_dir"]) / "01_database"
    with out.open("w", encoding="utf-8") as fw:
        for sp in config["species"]:
            pep_fasta = list((db_dir / sp["name"]).glob("*.pep.fasta"))[0]
            with pep_fasta.open("r", encoding="utf-8") as fr:
                fw.write(fr.read())
                if not str(fr).endswith("\n"):
                    fw.write("\n")
    return str(out)
