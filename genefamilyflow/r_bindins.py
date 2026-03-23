from __future__ import annotations

import logging
from pathlib import Path

from .utils import run_command


def run_r_script(
    script_path: str | Path, config_path: str | Path, logger: logging.Logger
) -> None:
    script = Path(script_path).expanduser().resolve()
    if not script.exists():
        raise FileNotFoundError(f"R script not found: {script}")
    run_command(["Rscript", str(script), "--config", str(config_path)], logger)
