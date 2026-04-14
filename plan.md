# GeneFamilyFlow — Follow-up Plan

Snapshot date: 2026-04-14
Test suite: 306 passing, 1 skipped, 94% script coverage.

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

## 2. Step 13 (GO/KEGG) — finish the offline path  ✅ done 2026-04-14

- **Precomputed eggNOG fixture shipped**:
  `example/13.GO_KEGG/Athaliana.emapper.annotations` and
  `example/13.GO_KEGG/Osativa.emapper.annotations` are hand-curated in
  real eggNOG-mapper `.emapper.annotations` format, with coverage for
  every family ID in the corresponding NPF fixture plus ~40 background
  genes. A dry-run of
  `snakemake --configfile config/default_config.yaml <overlay enabling
  step13> -n` now reports `step13_parse_eggnog: 2` and
  `step13_go_kegg: 2` jobs without error. `tests/test_parse_eggnog.py`
  has a parametrized `TestShippedFixtures` class that guarantees the
  fixtures stay parseable and family-complete.
- **Stale `annotations.tsv` removed** — the 2-row placeholder from
  `example/13.GO_KEGG/` was never consumed by any rule.
- **Placeholder string check stays in `validate_config.py`** — the
  test suite already covers it, and the Snakefile rule only runs when
  `run_eggnog: true`, so a second guard there would be redundant.
- **`docs/configuration.md` Step 13 section rewritten** to match the
  real schema (it still referenced a long-removed `annotation_file`
  key) and to call out the bundled fixture.

Follow-up for a future session: an actual
`snakemake --until step13_go_kegg -j 2 --use-conda` run inside the
full conda env, to confirm `R/13_go_kegg.R` produces non-empty
dotplots against the fixture (the DAG wiring is verified; the PDF
render has not been).

## 3. Coverage-gap audit on the remaining "low-cover" scripts

The three scripts still under 95% are all blocked by external
dependencies, not missing tests:

| Script | Cover | What's uncovered |
|---|---|---|
| `scan_promoter_fimo.py` | 69% | `run_fimo` subprocess invocation — needs the FIMO binary |
| `fetch_kegg_pathway_names.py` | 73% | `fetch()` HTTP call (already `# pragma: no cover`) |
| `fetch_jaspar_plants.py` | 98% | (resolved — `fetch_url` + `_atomic_write` now `# pragma: no cover`) |

Decision: **do not chase these with mocks**. `fetch_jaspar_plants.py`
has been aligned with the `fetch_kegg_pathway_names` convention —
`fetch_url`, `_atomic_write`, and the single call site inside `fetch()`
are marked `# pragma: no cover - network-dependent`. For
`scan_promoter_fimo`, the `TestScanEndToEnd` class already skips when
the binary is absent — that is the right shape; consider running it in
CI inside the conda env where FIMO exists.

## 4. CI — actually exercise the conda env  ✅ done 2026-04-14

- `.github/workflows/full-ci.yml` builds the full
  `envs/genefamily.yaml` via `mamba-org/setup-micromamba@v1`
  (`cache-environment: true` so repeat runs are cheap).
- Triggers: `workflow_dispatch` (manual), nightly cron at 05:17 UTC,
  and `pull_request` with the `full-ci` label attached. Scheduled and
  dispatched runs always fire; PR runs gate on the label name so the
  expensive matrix doesn't run on every push.
- Inside the full env it runs:
  1. `pytest tests/ -v --cov=scripts` — the `TestScanEndToEnd` class
     in `tests/test_scan_promoter_fimo.py` auto-unbolts its skipif
     guard because FIMO is now on PATH, so the real subprocess path
     gets exercised.
  2. `snakemake --configfile config/default_config.yaml -n -p` — full
     default DAG validation against installed tools.
  3. A second dry-run with a Step 13 overlay that flips
     `step13_go_kegg.enabled` to `true`, so the eggNOG fixture wiring
     added in item 2 stays live.
- Kept the existing `test.yml` workflow untouched so every push still
  gets the fast lane.

Follow-up: once the first scheduled run succeeds on `main`, add a
short note to `docs/development.md` pointing at the label trigger
(the rules file notes it but the dev doc doesn't call it out yet).

## 5. Documentation polish  ✅ done 2026-04-14

- README's test/coverage badge is now synced to "306 tests, 94% coverage".
- `docs/development.md` has a `# pragma: no cover` convention section
  listing the three current exemptions (`fetch_jaspar_plants.fetch_url`,
  `fetch_kegg_pathway_names.fetch`, `scan_promoter_fimo.run_fimo`) and
  a CI section explaining the fast-lane vs heavy-lane split.
- `docs/output.md` was audited against every `rule *.output` block in
  `rules/*.smk`. Fixed discrepancies:
  - **Step 4**: dropped phantom `identify.ID.txt` and
    `per_domain/{domain}.ID.txt`. Replaced with the real
    `work/04_identification/*` intermediates (merged IDs, combined
    FASTA, Pfam_scan.out, per-species family.ID).
  - **Step 6**: `phylogenetic_tree.pdf` → real pair
    `phylogenetic_tree_circular.pdf` + `phylogenetic_tree_rectangular.pdf`.
    Moved `alignment.aln` / `tree.nwk` / IQ-TREE log into the
    `work/06_tree/` section where they actually live.
  - **Step 7**: added `meme_info.txt` (produced alongside
    `meme_location.txt`) and moved `meme.html`/`meme.xml`/logos under
    `work/07_motif/meme_out/`.
  - **Step 8**: `kaks.tab.xls` is never under `output/` — routed to
    `example/8.collinearity/kaks.tab.xls` (precomputed) or
    `work/08_jcvi/kaks.tab.xls` (computed). Removed the phantom
    `synteny.pdf` (no rule produces it).
  - **Step 9**: `{species}.collinearity` and `{species}.tandem` live
    under `example/9.mcscanx/` or `work/09_mcscanx/`, never under
    `output/09_circos/`. Added the `{species}.gene_type` classifier
    output for the `precomputed: false` path.
  - **Step 11**: removed phantom `edge_annotation.xlsx` (it's an
    input, not an output) and `node_annotation.tsv` (no rule produces
    it). Only `PPI_network.pdf` is truly emitted.
  - **Step 13 + Step 14**: added entirely — these were missing from
    `docs/output.md` even though the rules have been live for a while.
    Step 13 section documents the enrichment PDF plus
    `{species}.go.tsv` / `.kegg.tsv` / `.universe.txt` intermediates
    and the eggNOG routing. Step 14 section covers
    `qRT_PCR.pdf` + `qRT_PCR_summary.csv`.
  - Updated the top-level directory tree with optional `13_go_kegg/`
    and `14_qrt_pcr/` entries and added a note clarifying that only
    the final FASTA and plotted PDFs are promoted into `output/`.

## 6. Refactor — shared FASTA/GFF3 parsing helpers  ✅ done 2026-04-14

- **New module `scripts/_parsers.py`** centralises the three primitives
  that used to be copy-pasted across four scripts:
  - `stream_fasta(path, *, error_cls=FastaFormatError)` — generator that
    validates headers and sequence blocks incrementally, yielding
    `RawFastaRecord` dataclass instances. Callers can pass their own
    `error_cls` subclass so `pytest.raises` still sees the module-local
    type.
  - `iter_gff3_rows(path, *, error_cls, strict_attributes)` — generator
    that yields `("comment", line)` or `("feature", RawGff3Row)` tuples,
    with the column-count / integer / range checks shared.
  - `parse_gff3_attributes(raw, *, strict, error_cls)` — dual-mode
    attribute parser. Lenient mode preserves the historical
    `filter_gff3` behaviour (silently skip items missing `=`); strict
    mode raises on both missing-`=` and empty-key items, matching the
    `extract_longest_transcript_from_gff3` contract.
- **Named `_parsers`, not `_io`, deliberately** — `_io` collides with
  CPython's builtin `_io` frozen module, so `from _io import …` inside
  a script run via `python3 scripts/foo.py` would resolve to the
  stdlib and silently break at runtime. Both import styles (package
  import `from scripts._parsers` and direct-script
  `from _parsers`) are handled via a small try/except at the top of
  each caller.
- **Four refactored scripts** now delegate their I/O:
  `scripts/read_protein_fasta.py`,
  `scripts/filter_longest_transcript.py`,
  `scripts/filter_gff3.py`,
  `scripts/extract_longest_transcript_from_gff3.py`. Each keeps its
  existing public dataclasses (`ProteinRecord`, `FastaRecord`, both
  `Gff3Feature` variants, `LongestTranscript*`) and its existing
  public functions (`parse_fasta`, `parse_gff3`, `_parse_attributes`,
  `filter_gff3_lines`, `collect_longest_transcripts`, …) as thin
  wrappers so the test contracts in `tests/test_*.py` stay unchanged.
- **Tests still 306 passing / 1 skipped.** The four targeted test
  files (`test_read_protein_fasta`, `test_filter_longest_transcript`,
  `test_filter_gff3`, `test_extract_longest_transcript_from_gff3`) —
  78 tests total — all pass, including the `pytest.raises(match=…)`
  assertions that pin exact error strings. A standalone CLI smoke
  test (`python3 scripts/read_protein_fasta.py` +
  `python3 scripts/filter_longest_transcript.py` with a small FASTA)
  also exercises the `ImportError` fallback branch end-to-end.

Net: roughly ~100 lines of duplicated parsing boilerplate collapsed
into `_parsers.py`, with no public API changes and no test churn.

## 7. Things explicitly not in scope

- A new pipeline step. The current 14 steps are stable.
- Migrating from Snakemake 8 to anything else.
- Replacing IQ-TREE / FastTree / MEME / MCScanX with newer tools.
- Adding a web UI.

These are all valid futures, but none of them fall out of the recent
work and none should sneak in as "while I'm here" changes.
