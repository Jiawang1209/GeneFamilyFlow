# GeneFamilyFlow Architecture

## Goal

GeneFamilyFlow provides a reproducible, configurable framework for multi-species gene family analysis by combining:

- Python for orchestration and command-line control
- External bioinformatics tools for sequence and evolutionary analysis
- R for downstream statistics and visualization

## Core Layout

- `genefamilyflow/cli.py`: CLI entry (`init`, `run-step`, `run`)
- `genefamilyflow/config.py`: YAML loading and validation
- `genefamilyflow/utils.py`: shared execution and parsing helpers
- `genefamilyflow/steps/`: step-by-step analysis modules
- `R/`: plotting/statistical scripts
- `R/report/GeneFamilyFlow_Report.Rmd`: report template
- `config/default_config.yaml`: user-editable configuration template

## Data Flow

1. User prepares `config.yaml` from default template.
2. CLI validates config and resolves step order.
3. Step modules generate outputs in `work/XX_step_name`.
4. R scripts consume step outputs and generate figures/tables.
5. Rmarkdown template collects all results into a final report.

## Execution Model

- `genefamilyflow init --output config.yaml`
- `genefamilyflow run --config config.yaml`
- `genefamilyflow run-step step06_tree --config config.yaml`

## Design Principles

- Keep each step independently runnable.
- Separate orchestration from plotting logic.
- Use config-driven behavior instead of hard-coded species/tool paths.
- Store intermediate files for reproducibility and troubleshooting.
