from __future__ import annotations

import shutil
from pathlib import Path

import click

from .config import load_config
from .steps import STEP_REGISTRY
from .utils import ensure_dir, get_logger


@click.group()
def main() -> None:
    """GeneFamilyFlow command line interface."""


@main.command("init")
@click.option(
    "--output",
    "output_path",
    default="config.yaml",
    show_default=True,
    help="Path to write an initial config file.",
)
def init_command(output_path: str) -> None:
    """Create a default configuration file."""
    root = Path(__file__).resolve().parent.parent
    template = root / "config" / "default_config.yaml"
    output = Path(output_path).expanduser().resolve()
    ensure_dir(output.parent)
    if output.exists():
        raise click.ClickException(f"Config file already exists: {output}")
    shutil.copyfile(template, output)
    click.echo(f"Created config template: {output}")


@main.command("run-step")
@click.argument("step_name")
@click.option("--config", "config_path", default="config.yaml", show_default=True)
def run_step(step_name: str, config_path: str) -> None:
    """Run one pipeline step by step name."""
    cfg = load_config(config_path).data
    logger = get_logger("genefamilyflow")
    if step_name not in STEP_REGISTRY:
        raise click.ClickException(
            f"Unknown step '{step_name}'. Available: {', '.join(STEP_REGISTRY.keys())}"
        )
    logger.info("Running step: %s", step_name)
    STEP_REGISTRY[step_name](cfg, logger)
    logger.info("Completed step: %s", step_name)


@main.command("run")
@click.option("--config", "config_path", default="config.yaml", show_default=True)
@click.option(
    "--from-step",
    "from_step",
    default=None,
    help="Run from this step (inclusive).",
)
def run_pipeline(config_path: str, from_step: str | None) -> None:
    """Run all pipeline steps in order."""
    cfg = load_config(config_path).data
    logger = get_logger("genefamilyflow")
    steps = list(STEP_REGISTRY.keys())
    if from_step:
        if from_step not in STEP_REGISTRY:
            raise click.ClickException(f"Unknown from-step '{from_step}'")
        steps = steps[steps.index(from_step) :]
    for step_name in steps:
        logger.info("Running step: %s", step_name)
        STEP_REGISTRY[step_name](cfg, logger)
        logger.info("Completed step: %s", step_name)
