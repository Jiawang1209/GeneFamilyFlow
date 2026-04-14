# GeneFamilyFlow

[![tests](https://github.com/Jiawang1209/GeneFamilyFlow/actions/workflows/test.yml/badge.svg)](https://github.com/Jiawang1209/GeneFamilyFlow/actions/workflows/test.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Snakemake](https://img.shields.io/badge/snakemake-≥8.0-brightgreen.svg)](https://snakemake.github.io)
[![Python](https://img.shields.io/badge/python-≥3.11-blue.svg)](https://www.python.org/)

Automated multi-species gene family analysis pipeline built with Snakemake.

GeneFamilyFlow integrates HMM search, BLAST, Pfam domain verification, phylogenetic tree construction, motif discovery, synteny analysis, promoter analysis, and PPI network visualization into a single reproducible workflow.

## Features

- **Multi-domain support** — search multiple Pfam domains independently, merge results into one gene family set
- **Multi-species** — analyze any number of species in parallel via `{species}` wildcard
- **Flexible tree building** — choose IQ-TREE (accurate) or FastTree (fast) via config
- **Auto-subfamily clustering** — top-down monophyletic clade splitting, with branch coloring, clade highlights, and radial subfamily labels rendered in both circular and rectangular layouts
- **Precomputed mode** — skip expensive computation steps (JCVI, MCScanX) when you already have results
- **Fully configurable** — all parameters in a single YAML config file
- **Tested** — unit tests for all Python scripts

## Pipeline Overview

| Step | Description | Tools |
|------|-------------|-------|
| 1 | Data preparation (filter longest transcript) | seqkit |
| 2 | HMM search (per domain, 2 rounds) | HMMER, ClustalW |
| 3 | BLAST search (per domain) | BLAST+ |
| 4 | Gene family identification (merge + Pfam verify) | pfam_scan.pl |
| 5 | Gene family statistics | R (Biostrings, Peptides) |
| 6 | Phylogenetic tree + auto-subfamily clustering | MUSCLE + IQ-TREE/FastTree + ggtree |
| 7 | Motif & gene structure | MEME + R (ggtree, gggenes) |
| 8 | JCVI synteny & Ka/Ks | JCVI + R |
| 9 | MCScanX synteny & Circos | MCScanX + R (circlize) |
| 10 | Promoter cis-element analysis | bedtools + FIMO/JASPAR Plants (offline default) + R |
| 11 | PPI network | R (ggNetView) |

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/Jiawang1209/GeneFamilyFlow.git
cd GeneFamilyFlow
```

### 2. Create Conda environment

```bash
# Full environment (all tools)
conda env create -f envs/genefamily.yaml
conda activate GeneFamilyFlow

# Or minimal environment (Python + Snakemake only, for development)
conda env create -f envs/genefamily-minimal.yaml
conda activate GeneFamilyFlow
```

### 3. Install tools not available via Conda

These tools require manual installation:

| Tool | Purpose | Installation |
|------|---------|-------------|
| MCScanX | Synteny detection (Step 9) | [GitHub](https://github.com/wyp1125/MCScanX) — compile from source |
| KaKs_Calculator | Ka/Ks calculation (Step 8) | [SourceForge](https://sourceforge.net/projects/kakscalculator2/) |
| ggNetView | PPI network visualization (Step 11) | `Rscript -e 'install.packages("ggNetView")'` |

Step 10 runs **fully offline by default** using FIMO (already in `envs/genefamily.yaml`) against the JASPAR 2024 CORE plants non-redundant motif bundle (~805 plant TF PFMs), shipped at `example/10.promoter/JASPAR2024_plants.meme`. Refresh from `jaspar.elixir.no` with `python scripts/fetch_jaspar_plants.py`. The legacy `plantcare` (manual web upload) and `local` (literal-substring fallback) modes are still selectable via `step10.scan_method`.

## Usage

### Prepare input data

Place your data in the database directory (default: `example/1.database/`):

```
example/1.database/
├── Sbicolor.pep.fasta        # Protein sequences (required)
├── Athaliana.pep.fasta
├── Osativa.pep.fasta
├── Sbicolor.gff3             # Gene annotations (for JCVI/MCScanX computation)
├── Sbicolor.cds.fasta        # CDS sequences (for JCVI computation)
└── Sbicolor.genome.fasta     # Genome sequences (for promoter extraction)
```

Provide the HMM model for each Pfam domain. BLAST seed sequences are extracted automatically from reference species at runtime (see `step3.reference_species`):

```
example/2.hmmsearch/PF00854.hmm                          # Pfam HMM model
example/3.blast/references/Athaliana.domains.tsv         # 2-col clean table (committed)
example/3.blast/references/Athaliana.pep.fasta           # gene-level proteome (user-provided)
example/3.blast/references/Osativa.domains.tsv           # 2-col clean table (committed)
example/3.blast/references/Osativa.pep.fasta             # gene-level proteome (user-provided)
```

### Configure the pipeline

Edit `config/default_config.yaml`:

```yaml
# Species to analyze
species:
  - Sbicolor
  - Athaliana
  - Osativa

target_species: "Sbicolor"

# Pfam domains to search (add multiple for multi-domain analysis)
step2:
  domains:
    - pfam_id: "PF00854"
      pfam_hmm_file: "example/2.hmmsearch/PF00854.hmm"
      seed_references:
        - Athaliana
    # - pfam_id: "PF00005"
    #   pfam_hmm_file: "path/to/PF00005.hmm"
    #   seed_references: [Athaliana, Osativa]

# Reference species used by rule step03_extract_seed to build the BLAST seed
step3:
  reference_species:
    Athaliana:
      domains_table: "example/3.blast/references/Athaliana.domains.tsv"
      proteome:      "example/3.blast/references/Athaliana.pep.fasta"
    Osativa:
      domains_table: "example/3.blast/references/Osativa.domains.tsv"
      proteome:      "example/3.blast/references/Osativa.pep.fasta"

# Tree tool: "iqtree" (accurate, slow) or "fasttree" (fast, approximate)
step6:
  tree_tool: "iqtree"
  # Auto-split tips into N monophyletic subfamilies (0 disables).
  auto_n_subfamilies: 8
  auto_subfamily_prefix: "Sub"
  # Render "circular", "rectangular", or "both" PDFs
  tree_layout: "both"

# Use precomputed results for synteny (set to false to compute from scratch)
step8:
  precomputed: true
step9:
  precomputed: true
```

### Run the pipeline

```bash
# Dry run — check the DAG without executing
snakemake --configfile config/default_config.yaml -n -p

# Full run
snakemake --configfile config/default_config.yaml -j 8 --use-conda

# Run up to a specific step
snakemake --configfile config/default_config.yaml --until step06_tree_build -j 4

# Run only steps 1-4 (gene family identification)
snakemake --configfile config/default_config.yaml --until step04_pfam_filter -j 8
```

### Run tests

```bash
pytest tests/ -v
```

## Output

Results are organized in the `output/` directory:

```
output/
├── 02_hmmsearch/
│   └── {domain}/
│       ├── {species}.pfam.domtblout  # First-round HMM search results
│       ├── custom.hmm                # Custom-built HMM (second round)
│       └── hmm_2nd_ids.txt           # Final HMM candidate IDs
├── 03_blast/
│   └── {domain}/
│       ├── {species}.blast           # Raw BLASTP hits per species
│       └── blast_ids.txt             # Final BLAST candidate IDs
├── 04_identification/
│   └── identify.ID.clean.fa          # Final gene family sequences
├── 05_genefamily_info/
│   ├── Gene_Information.xlsx          # Gene properties table
│   └── Gene_Information.csv
├── 06_tree/
│   ├── phylogenetic_tree_circular.pdf    # Circular layout with subfamily hilights
│   └── phylogenetic_tree_rectangular.pdf # Rectangular layout with subfamily hilights
├── 07_motif/
│   ├── meme_location.txt              # Motif coordinates
│   └── Tree_Domain_Motif_GeneStructure.pdf  # Composite figure
├── 08_jcvi_kaks/
│   └── 08_jcvi_kaks.pdf              # Ka/Ks distribution
├── 09_circos/
│   └── circos_{species}.pdf           # Circos plots per species
├── 10_promoter/
│   └── promoter_elements.pdf          # Cis-element analysis
└── 11_ppi/
    └── PPI_network.pdf                # PPI network
```

## Documentation

Full documentation lives in [`docs/`](docs/):

| Document | Description |
|----------|-------------|
| [docs/configuration.md](docs/configuration.md) | Complete parameter reference for `config/default_config.yaml` |
| [docs/tutorial.md](docs/tutorial.md) | Step-by-step walkthrough using the bundled example data |
| [docs/output.md](docs/output.md) | Description of every output file and how to interpret it |
| [docs/troubleshooting.md](docs/troubleshooting.md) | Common errors and fixes (install, runtime, external tools) |
| [docs/development.md](docs/development.md) | Contributor guide: coding style, testing, adding new steps |

## Contributing

Contributions are welcome. Please read [docs/development.md](docs/development.md) before submitting a PR. All pull requests must:

- Pass `pytest tests/ -v` (currently 306 tests, 94% coverage)
- Pass `snakemake --configfile config/default_config.yaml -n` dry-run
- Include tests for new Python code
- Update relevant docs when changing config or output

## Citation

If you use GeneFamilyFlow in your research, please cite this repository:

```
Jiawang1209. GeneFamilyFlow: Automated multi-species gene family analysis pipeline.
GitHub repository, https://github.com/Jiawang1209/GeneFamilyFlow
```

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.
