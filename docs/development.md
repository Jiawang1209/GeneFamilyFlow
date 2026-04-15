# Development Guide

Contributing to GeneFamilyFlow — coding style, testing, and how to add new analysis steps.

## Project Layout

```
GeneFamilyFlow/
├── Snakefile                  # Main workflow (all rules)
├── config/default_config.yaml # Centralized parameters
├── envs/                      # Conda environment definitions
│   ├── genefamily.yaml        # Full env (all tools)
│   └── genefamily-minimal.yaml # Dev/test env
├── scripts/                   # Python utilities (stdlib only)
├── R/                         # R visualization scripts
├── tests/                     # pytest test suite (389 tests, 97% cov)
├── docs/                      # Documentation
├── example/                   # Example data for the tutorial
├── output/                    # Pipeline outputs (gitignored)
└── work/                      # Intermediate files (gitignored)
```

## Coding Style

### Python (`scripts/`)

- **Standard library only** — no third-party imports in `scripts/` to keep tests fast and CI lightweight
- **Full type annotations**: `str | Path`, `list[T]`, `Iterator[T]`
- **Immutable data classes**: `@dataclass(frozen=True)`
- **Argparse CLI**: every script exposes `main(argv) -> int` and is callable as `python scripts/foo.py --help`
- **Streaming for large files**: use generators, not `file.read()`
- **Every script must have a matching `tests/test_*.py`**

Example scaffolding:

```python
from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator


@dataclass(frozen=True)
class Record:
    gene_id: str
    value: float


def parse(path: Path) -> Iterator[Record]:
    with path.open() as fh:
        for line in fh:
            ...


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    args = parser.parse_args(argv)
    ...
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

### R (`R/`)

- **`optparse` CLI** for Snakemake integration
- **tidyverse + ggplot2** style
- **Fail fast** with `stop()` for missing required arguments
- **Output PDF to `--outdir`**

Example:

```r
suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
})

option_list <- list(
  make_option("--input", type = "character", default = NULL),
  make_option("--outdir", type = "character", default = "output")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input)) stop("--input is required")

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
...
ggsave(file.path(opt$outdir, "figure.pdf"), plot = p, width = 10, height = 8)
```

### Snakemake rules

- **One rule per analysis step**, named `stepN_description`
- **Paths from config**, never hardcoded
- **Conditional rules** with `if/else` (e.g. `tree_tool`, `precomputed`)
- **Use `temp()` for large intermediate files**
- **Always `log:` block**

Example:

```python
rule step_example:
    input:
        fasta = f"{WORK_DIR}/input/{{species}}.fa",
    output:
        result = f"{OUT_DIR}/example/{{species}}.txt",
    log:
        f"{LOG_DIR}/example_{{species}}.log",
    threads: config.get("step_example", {}).get("threads", 4)
    shell:
        """
        some_tool -i {input.fasta} -o {output.result} 2>&1 | tee {log}
        """
```

## Testing

### Run the test suite

```bash
pytest tests/ -v
# 389 passed, 1 skipped in ~11s
```

### Coverage

```bash
pytest tests/ --cov=scripts --cov-report=term-missing
```

Target: **80%+ coverage** (currently 389 tests, ~97% coverage).

### `# pragma: no cover` convention

Network calls and subprocess invocations that exist purely to shell out to
an external binary or HTTP endpoint are exempt from coverage. Mark the
function body (or the single call) with `# pragma: no cover - <reason>`
so the coverage report doesn't flag it as untested code. Current
exemptions:

- `scripts/fetch_jaspar_plants.py::fetch_url` — HTTP call to JASPAR
- `scripts/fetch_kegg_pathway_names.py::fetch` — HTTP call to KEGG REST
- `scripts/scan_promoter_fimo.py::run_fimo` — subprocess to the FIMO binary

The end-to-end path for these still gets exercised in the conda-based CI
job where the real tool is available (see `TestScanEndToEnd` in
`tests/test_scan_promoter_fimo.py`). Do not replace them with
subprocess/urllib mocks in the unit suite — that produces brittle tests
without real signal.

**What is allowed:** stubbing our _own_ higher-level wrappers to cover
orchestration logic above the exempt bodies. For example,
`tests/test_scan_promoter_fimo.py::TestScanOrchestration` monkeypatches
`scan_promoter_fimo.run_fimo` (our function, not `subprocess.run`) to
test the tempdir/parse/write flow in `scan()`. Likewise
`tests/test_fetch_kegg_pathway_names.py::TestMainCli` stubs
`fetch_kegg_pathway_names.fetch` (ours, not `urllib.request.urlopen`)
to test the argparse/parse/write flow in `main()`. This keeps the
network/subprocess seam honest while still exercising the Python we own.
Module entrypoints (`if __name__ == "__main__": raise SystemExit(main())`)
are covered via `runpy.run_module(..., run_name="__main__")` with
`sys.argv = [..., "--help"]` so argparse short-circuits before any
external call.

### Test structure

- Each script has a matching `tests/test_{script}.py`
- Use `tests/data/` for fixture files
- AAA pattern: Arrange, Act, Assert
- Descriptive test names: `test_returns_longest_when_multiple_transcripts_present`

Example:

```python
from pathlib import Path
from scripts.parse_hmmsearch import parse_domtblout


def test_parses_valid_domtblout(tmp_path: Path) -> None:
    # Arrange
    content = "geneA\t-\t100\tPF00854.hmm\t...\n"
    fixture = tmp_path / "test.domtblout"
    fixture.write_text(content)

    # Act
    records = list(parse_domtblout(fixture))

    # Assert
    assert len(records) == 1
    assert records[0].gene_id == "geneA"
```

### Continuous Integration

Two workflows live under `.github/workflows/`:

- **`test.yml`** — fast lane. Runs on every push/PR to `main` across a
  Python 3.11/3.12 matrix. Installs only `pytest` + `pytest-cov` so it
  finishes in well under a minute. FIMO, hmmer, blast, MEME, and the R
  stack are all absent, so tests that need external tools (e.g.
  `TestScanEndToEnd` in `tests/test_scan_promoter_fimo.py`) auto-skip.

- **`full-ci.yml`** — heavy lane. Builds the full
  `envs/genefamily.yaml` via `setup-micromamba` and then runs the
  complete pytest suite, a `snakemake -n` dry-run against the default
  config, and a second dry-run with Step 13 enrichment enabled. The
  JASPAR/FIMO end-to-end test unbolts its skip guard because FIMO is
  on PATH inside this env. Triggers:

  | Trigger | When |
  |---|---|
  | `workflow_dispatch` | Manual run from the Actions tab |
  | `schedule` (cron `17 5 * * *`) | Nightly 05:17 UTC |
  | `pull_request` with label `full-ci` | Attach the label to a PR to run on the PR head |

  The PR trigger gates on `github.event.label.name == 'full-ci'` so
  ad-hoc label edits don't trigger spurious runs.

## Adding a New Analysis Step

1. **Add the Python utility** (if needed):
   ```bash
   touch scripts/my_new_tool.py
   touch tests/test_my_new_tool.py
   ```
   Write the test first (TDD), then implement.

2. **Add the R visualization** (if needed):
   ```bash
   touch R/15_my_new_step.R
   ```
   Use `optparse` for CLI parameters.

3. **Add the Snakemake rule** in `Snakefile`:
   - Place in step order (search for the step above it)
   - Read paths from `config["step15"]`
   - Add output to `rule all`'s input list

4. **Add config parameters** to `config/default_config.yaml`:
   ```yaml
   step15:
     r_script: "R/15_my_new_step.R"
     figure_width: 10
     figure_height: 8
   ```

5. **Dry-run and verify the DAG**:
   ```bash
   snakemake --configfile config/default_config.yaml -n -p
   ```

6. **Run the tests**:
   ```bash
   pytest tests/ -v
   ```

7. **Update docs**:
   - `docs/configuration.md` — new parameters
   - `docs/output.md` — new output files
   - `README.md` — step table

## Commit Style

Follow conventional commits:

```
feat: add step15 for GO enrichment analysis
fix: handle empty HMM hits in parse_hmmsearch
docs: add tutorial for custom datasets
test: cover edge cases in parse_pfam_scan
refactor: extract shared FASTA parser
```

## Pre-commit Checklist

- [ ] `pytest tests/ -v` passes
- [ ] `snakemake --configfile config/default_config.yaml -n` dry-run succeeds
- [ ] New code has tests
- [ ] Docs updated (if config/output changed)
- [ ] No hardcoded paths

## Resources

- [Snakemake docs](https://snakemake.readthedocs.io/)
- [Pfam database](http://pfam.xfam.org/)
- [HMMER manual](http://eddylab.org/software/hmmer/Userguide.pdf)
- [JCVI tutorial](https://github.com/tanghaibao/jcvi/wiki)
- [MCScanX manual](https://github.com/wyp1125/MCScanX)
