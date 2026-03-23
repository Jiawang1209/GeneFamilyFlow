# GeneFamilyFlow

GeneFamilyFlow is a configurable framework for multi-species gene family analysis.
It combines Python-based workflow orchestration with R-based statistical analysis
and visualization.

## Features

- Multi-species, config-driven workflow
- CLI with full-pipeline and single-step execution
- HMM + BLAST gene family identification
- Tree, motif, collinearity, MCScanX, promoter, and PPI integration points
- Rmarkdown report template for final result reporting

## Project Structure

- `genefamilyflow/`: Python package and pipeline steps
- `config/default_config.yaml`: configuration template
- `R/`: R scripts for plotting and statistics
- `R/report/GeneFamilyFlow_Report.Rmd`: report template
- `doc/`: architecture, pipeline, and configuration documents
- `example/`: original example datasets and scripts

## Installation

### Conda Environment

```bash
conda env create -f environment.yml
conda activate genefamilyflow
```

### Python Package

```bash
pip install -e .
```

## Quickstart

1. Initialize config template:

```bash
genefamilyflow init --output config.yaml
```

2. Edit `config.yaml` and set:
   - species file paths (`pep_fasta`, `cds_fasta`, `gff3`, optional `genome_fasta`)
   - `reference_species.family_member_ids`
   - `pfam.hmm_file`
   - tool executable names/paths if needed

3. Run one step for validation:

```bash
genefamilyflow run-step step01_database --config config.yaml
```

4. Run full pipeline:

```bash
genefamilyflow run --config config.yaml
```

## Analysis Workflow (Demo)

This repository currently provides a demo-oriented but runnable pipeline layout.
The workflow is adapted from `example/GeneFamilyFlow.md` and organized into
11 executable steps.

### Input Data Requirements

For each species:

- Protein FASTA (`*.pep.fasta`) - required
- CDS FASTA (`*.cds.fasta`) - required
- GFF3/GTF annotation (`*.gff3` or `*.gtf`) - required
- Genome FASTA (`*.fa`/`*.fasta`) - recommended

Reference data:

- Known family member IDs in a reference species (for BLAST seed set)
- Pfam HMM model file (for domain-based candidate detection)

### Step-by-Step Pipeline

#### Step01 `step01_database`

- Standardize per-species source files into `work/01_database/{species}`.
- Ensure all required inputs exist before downstream analysis.

Run:

```bash
genefamilyflow run-step step01_database --config config.yaml
```

#### Step02 `step02_hmmsearch`

- Round1 HMM search on all species proteins.
- Filter hits by E-value threshold.
- Build custom HMM from round1 alignments.
- Round2 HMM search with custom model.
- Output key ID list: `work/02_hmmsearch/2st_id`.

Run:

```bash
genefamilyflow run-step step02_hmmsearch --config config.yaml
```

#### Step03 `step03_blast`

- Build BLAST database from known reference-family proteins.
- BLAST all species proteins against reference DB.
- Output candidate IDs: `work/03_blast/species.blast.id`.

Run:

```bash
genefamilyflow run-step step03_blast --config config.yaml
```

#### Step04 `step04_identification`

- Intersect HMM and BLAST candidates.
- Validate by Pfam domain (`hmmscan`).
- Produce final family sequences:
  - `work/04_identification/identify.ID`
  - `work/04_identification/identify.ID.fa`

Run:

```bash
genefamilyflow run-step step04_identification --config config.yaml
```

#### Step05 `step05_genefamily_info`

- Merge annotation into BED-like gene coordinates.
- Generate family basic statistics (length, MW, pI, hydrophobicity) using R.
- Main outputs under `work/05_genefamily_info/`.

#### Step06 `step06_tree`

- MUSCLE alignment of identified proteins.
- IQ-TREE phylogeny inference.
- R script generates tree figures.

#### Step07 `step07_motif`

- MEME motif discovery on identified proteins.
- Parse motif locations and generate motif structure plot.

#### Step08 `step08_collinearity`

- JCVI-oriented collinearity analysis integration.
- Ka/Ks plotting entry via `R/08_jcvi_kaks.R`.

#### Step09 `step09_mcscanx`

- MCScanX-oriented duplication/synteny integration.
- Circos plotting entry via `R/09_circos.R`.

#### Step10 `step10_promoter`

- Promoter analysis integration (upstream regions, cis-element plotting entry).

#### Step11 `step11_ppi`

- PPI analysis integration (ortholog-based mapping and network plotting entry).

### Optional Downstream Analyses (R)

- `R/12_rnaseq.R`: RNA-seq expression heatmap
- `R/13_go_kegg.R`: GO/KEGG enrichment plotting
- `R/14_qrt_pcr.R`: qRT-PCR visualization

### Run Full Demo Pipeline

```bash
genefamilyflow run --config config.yaml
```

Or resume from a specific step:

```bash
genefamilyflow run --config config.yaml --from-step step05_genefamily_info
```

## CLI Commands

- `genefamilyflow init --output config.yaml`
- `genefamilyflow run-step <step_name> --config config.yaml`
- `genefamilyflow run --config config.yaml [--from-step step05_genefamily_info]`

Available step names:

- `step01_database`
- `step02_hmmsearch`
- `step03_blast`
- `step04_identification`
- `step05_genefamily_info`
- `step06_tree`
- `step07_motif`
- `step08_collinearity`
- `step09_mcscanx`
- `step10_promoter`
- `step11_ppi`

## Documentation

- `doc/architecture.md`
- `doc/pipeline_steps.md`
- `doc/configuration_guide.md`

## Notes

- Current status is a framework demo with complete step interfaces.
- Some advanced biological logic (strict production behavior) is still being
  expanded from the original `example/GeneFamilyFlow.md` protocol.
- You can progressively harden each step using your real datasets and species
  naming rules.
