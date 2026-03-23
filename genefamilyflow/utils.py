from __future__ import annotations

import logging
import subprocess
from pathlib import Path
from typing import Iterable

from Bio import SeqIO


def get_logger(name: str, level: int = logging.INFO) -> logging.Logger:
    logger = logging.getLogger(name)
    if not logger.handlers:
        handler = logging.StreamHandler()
        formatter = logging.Formatter(
            fmt="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )
        handler.setFormatter(formatter)
        logger.addHandler(handler)
    logger.setLevel(level)
    return logger


def ensure_dir(path: str | Path) -> Path:
    out = Path(path)
    out.mkdir(parents=True, exist_ok=True)
    return out


def run_command(cmd: list[str], logger: logging.Logger, cwd: str | Path | None = None) -> None:
    logger.info("Executing: %s", " ".join(cmd))
    subprocess.run(cmd, check=True, cwd=cwd)


def write_lines(path: str | Path, lines: Iterable[str]) -> None:
    p = Path(path)
    with p.open("w", encoding="utf-8") as fh:
        for line in lines:
            fh.write(f"{line}\n")


def read_fasta_ids(path: str | Path) -> list[str]:
    ids: list[str] = []
    for record in SeqIO.parse(str(path), "fasta"):
        ids.append(record.id)
    return ids


def clean_gene_id(gene_id: str) -> str:
    """Normalize representative IDs from common suffix patterns."""
    out = gene_id
    if ".MSU" in out:
        out = out.split(".MSU", maxsplit=1)[0]
    if out.endswith(".v3.1"):
        out = out[: -len(".v3.1")]
    if ".p" in out:
        left, right = out.split(".p", maxsplit=1)
        if right.isdigit():
            out = left
    return out


def unique_keep_order(items: Iterable[str]) -> list[str]:
    seen: set[str] = set()
    out: list[str] = []
    for item in items:
        if item in seen:
            continue
        seen.add(item)
        out.append(item)
    return out
