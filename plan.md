# GeneFamilyFlow — Follow-up Plan

Snapshot date: 2026-04-14
Test suite: 304 passing, 1 skipped, 93% script coverage.

This document captures concrete next steps that emerged while finishing
the Step 10 offline rewrite and the test-coverage backfill. Items are
ordered by expected payoff, not by required sequence — most are
independent and can be picked up in any order.

---

## 1. Validate the Step 10 rewrite end-to-end

The offline JASPAR + FIMO scanner is wired in and unit-tested, but it
has not yet been exercised against real promoter sequences inside the
full pipeline.

- Run `snakemake --configfile config/default_config.yaml -n -p` and
  confirm the DAG includes `step10_fimo_scan` for the default
  `scan_method=jaspar`.
- Then run `snakemake --until step10_promoter -j 4 --use-conda` against
  the `example/` data and inspect:
  - `work/10_promoter/jaspar_fimo/plantCARE_output_jaspar.tab`
  - `output/10_promoter/promoter_elements.pdf`
- Confirm `R/10_promoter.R` consumes the new `.tab` file unchanged
  (this was the explicit design constraint — the R script must not need
  to know which `scan_method` produced the input).
- Open an issue if the R rendering breaks or the heatmap categories
  collapse because JASPAR motif names don't map cleanly onto PlantCARE
  categories; in that case we may need a small motif-to-category lookup
  table alongside `cir_element.desc.20240509.xlsx`.

## 2. Step 13 (GO/KEGG) — finish the offline path

The KEGG pathway-name fixture is already shipped (`build_kegg_names.py`,
`fetch_kegg_pathway_names.py`), but the eggNOG side still requires the
user to provide either a real eggNOG data dir or a precomputed
annotations dir. Two things would push step 13 over the line:

- A **precomputed** eggNOG fixture for the bundled `example/` species
  (just enough rows to render the GO/KEGG figures), so `enabled: true`
  works out of the box without an external download.
- Decide whether the placeholder string check
  (`/path/to/eggnog_data`) belongs in `validate_config.py` or in the
  Snakefile rule itself. Currently it lives in `validate_config.py`
  and the test suite covers it — leave it there unless the rule needs
  the same guard.

## 3. Coverage-gap audit on the remaining "low-cover" scripts

The three scripts still under 95% are all blocked by external
dependencies, not missing tests:

| Script | Cover | What's uncovered |
|---|---|---|
| `scan_promoter_fimo.py` | 69% | `run_fimo` subprocess invocation — needs the FIMO binary |
| `fetch_kegg_pathway_names.py` | 73% | `fetch()` HTTP call (already `# pragma: no cover`) |
| `fetch_jaspar_plants.py` | 78% | `fetch_url()` HTTP call |

Decision: **do not chase these with mocks**. Instead, add an explicit
`# pragma: no cover` to the `fetch_url` body in `fetch_jaspar_plants.py`
to match the convention already in place in `fetch_kegg_pathway_names`,
so the coverage report doesn't make it look like real test debt. For
`scan_promoter_fimo`, the `TestScanEndToEnd` class already skips when
the binary is absent — that is the right shape; consider running it in
CI inside the conda env where FIMO exists.

## 4. CI — actually exercise the conda env

Right now `pytest tests/ -v` runs in CI but FIMO, hmmer, blast, and the
R stack are all skipped or mocked. A nightly job (or a dedicated
workflow gated on a `[full-ci]` PR label) should:

- Build the full `envs/genefamily.yaml` env
- Run `pytest tests/ -m "unit or integration"` so the
  `TestScanEndToEnd::test_fimo_run_against_real_binary` test stops
  being skipped
- Run `snakemake -n -p --configfile config/default_config.yaml`
  inside that env

This is the only way the JASPAR/FIMO path gets continuous coverage
without bloating the unit tests with subprocess mocks.

## 5. Documentation polish

- README's commit history mentions the `currently N tests, X% coverage`
  badge — keep it in sync. (Just bumped to 304 tests, 93% coverage.)
- `docs/development.md` should mention the `# pragma: no cover`
  convention for network/subprocess code so future contributors don't
  get nagged by reviewers about missing branches.
- `docs/output.md` is now accurate for Step 10. Sanity-check the other
  step tables match what the rules actually emit — easy to check by
  running `find output -type f` after a real run.

## 6. Refactor opportunity (optional)

Several `scripts/*.py` modules duplicate the same FASTA / GFF3 parsing
boilerplate (`read_protein_fasta`, `filter_longest_transcript`,
`extract_longest_transcript_from_gff3`, `filter_gff3`). Each has its
own `parse_*` generator with the same error-message style and the same
"missing file / not a file / no records" branches.

Pulling these into a tiny shared module (`scripts/_io.py` with
`stream_fasta()` and `parse_gff3_lines()`) would:

- Eliminate ~60 lines of duplication
- Make any future error-message changes one-shot
- Let the existing tests stay where they are (the helper would live
  behind the same public API)

This is **not urgent**. Bring it up only if a third script wants the
same parser, otherwise the duplication is cheap and clear.

## 7. Things explicitly not in scope

- A new pipeline step. The current 14 steps are stable.
- Migrating from Snakemake 8 to anything else.
- Replacing IQ-TREE / FastTree / MEME / MCScanX with newer tools.
- Adding a web UI.

These are all valid futures, but none of them fall out of the recent
work and none should sneak in as "while I'm here" changes.
