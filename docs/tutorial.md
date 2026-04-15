# Tutorial: Running the Example Dataset

This tutorial walks through a complete run of GeneFamilyFlow using the bundled example data (NPF gene family across *Sorghum bicolor*, *Arabidopsis thaliana*, and *Oryza sativa*).

## Prerequisites

- Conda or Mamba installed
- ~10 GB free disk space
- ~30 minutes (with `-j 8`)

## Step 1: Clone and set up environment

```bash
git clone https://github.com/Jiawang1209/GeneFamilyFlow.git
cd GeneFamilyFlow

conda env create -f envs/genefamily.yaml
conda activate GeneFamilyFlow
```

## Step 2: Verify the example data

```bash
ls example/1.database/
# Sbicolor.pep.fasta  Athaliana.pep.fasta  Osativa.pep.fasta
# Sbicolor.gff3       Athaliana.gff3       Osativa.gff3
# Sbicolor.cds.fasta  Athaliana.cds.fasta  Osativa.cds.fasta
```

The default config (`config/default_config.yaml`) points to this directory.

## Step 3: Install external tools (if not using precomputed mode)

The default config uses `precomputed: true` for Steps 8 and 9, so you can skip MCScanX and KaKs_Calculator for the example run.

If you want to run from scratch, see [troubleshooting.md](troubleshooting.md#external-tools).

## Step 4: Download the Pfam database

Step 4 requires a local Pfam-A database for `pfam_scan.pl`:

```bash
mkdir -p ~/pfam_db && cd ~/pfam_db
wget http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
wget http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz
gunzip *.gz
hmmpress Pfam-A.hmm
```

Update `config/default_config.yaml`:

```yaml
step4:
  pfam_database_dir: "/Users/YOUR_USERNAME/pfam_db"
```

## Step 5: Dry-run to verify the DAG

```bash
snakemake --configfile config/default_config.yaml -n -p
```

You should see ~36 jobs scheduled. No errors = DAG is valid.

## Step 6: Run the pipeline

```bash
# Full run with 8 parallel jobs
snakemake --configfile config/default_config.yaml -j 8 --use-conda

# Or run up to a specific step
snakemake --configfile config/default_config.yaml --until step06_tree_build -j 4

# Run only gene family identification (Steps 1-4)
snakemake --configfile config/default_config.yaml --until step04_pfam_filter -j 8
```

## Step 7: Inspect results

```bash
ls output/
# 04_identification/  05_genefamily_info/  06_tree/  07_motif/
# 08_jcvi_kaks/       09_circos/           10_promoter/  11_ppi/
```

See [output.md](output.md) for a description of every output file.

## Step 8: Between-step manual curation (optional)

After Step 6 (tree), you may want to manually curate `tree_group_2.xlsx` to define subfamilies. This file is used by downstream visualization steps.

## Step 9: Re-running with custom data

To analyze your own gene family:

1. Copy input FASTA/GFF3 files to a new directory (e.g. `data/`).
2. Copy `config/default_config.yaml` → `config/my_config.yaml`.
3. Edit `my_config.yaml`:
   - `species`, `target_species`, `database_dir`
   - `step2.domains` (Pfam ID, HMM file, seed sequences)
   - `step9.species_config` (chromosome counts per species)
4. Dry-run: `snakemake --configfile config/my_config.yaml -n`
5. Run: `snakemake --configfile config/my_config.yaml -j 8 --use-conda`

See [configuration.md](configuration.md) for the full parameter reference.

## Running the tests

```bash
pytest tests/ -v
# 389 passed, 1 skipped in ~11s
```
