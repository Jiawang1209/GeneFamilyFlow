# GeneFamilyFlow

Automated multi-species gene family analysis pipeline built with Snakemake.

GeneFamilyFlow integrates HMM search, BLAST, Pfam domain verification, phylogenetic tree construction, motif discovery, synteny analysis, promoter analysis, and PPI network visualization into a single reproducible workflow.

## Features

- **Multi-domain support** — search multiple Pfam domains independently, merge results into one gene family set
- **Multi-species** — analyze any number of species in parallel via `{species}` wildcard
- **Flexible tree building** — choose IQ-TREE (accurate) or FastTree (fast) via config
- **Precomputed mode** — skip expensive computation steps (JCVI, MCScanX) when you already have results
- **Fully configurable** — all parameters in a single YAML config file
- **Tested** — 65 unit tests for all Python scripts

## Pipeline Overview

| Step | Description | Tools |
|------|-------------|-------|
| 1 | Data preparation (filter longest transcript) | seqkit |
| 2 | HMM search (per domain, 2 rounds) | HMMER, ClustalW |
| 3 | BLAST search (per domain) | BLAST+ |
| 4 | Gene family identification (merge + Pfam verify) | pfam_scan.pl |
| 5 | Gene family statistics | R (Biostrings, Peptides) |
| 6 | Phylogenetic tree | MUSCLE + IQ-TREE/FastTree |
| 7 | Motif & gene structure | MEME + R (ggtree, gggenes) |
| 8 | JCVI synteny & Ka/Ks | JCVI + R |
| 9 | MCScanX synteny & Circos | MCScanX + R (circlize) |
| 10 | Promoter cis-element analysis | bedtools + PlantCARE + R |
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

PlantCARE (Step 10) is an [online tool](https://bioinformatics.psb.ugent.be/webtools/plantcare/html/) — submit promoter sequences manually and download results.

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

Provide HMM and seed files for each Pfam domain:

```
example/2.hmmsearch/PF00854.hmm       # Pfam HMM model
example/3.blast/PF00854.TAIR.ID.fa    # Seed sequences for BLAST
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
      seed_sequence_file: "example/3.blast/PF00854.TAIR.ID.fa"
    # - pfam_id: "PF00005"
    #   pfam_hmm_file: "path/to/PF00005.hmm"
    #   seed_sequence_file: "path/to/PF00005.seed.fa"

# Tree tool: "iqtree" (accurate, slow) or "fasttree" (fast, approximate)
step6:
  tree_tool: "iqtree"

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
├── 04_identification/
│   └── identify.ID.clean.fa          # Final gene family sequences
├── 05_genefamily_info/
│   ├── Gene_Information.xlsx          # Gene properties table
│   └── Gene_Information.csv
├── 06_tree/
│   └── phylogenetic_tree.pdf          # Phylogenetic tree figure
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

## License

MIT
