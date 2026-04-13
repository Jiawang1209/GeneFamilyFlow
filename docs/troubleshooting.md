# Troubleshooting

Common errors and fixes. If your issue isn't listed, check `logs/<rule_name>.log` for the full error message, then [open an issue](https://github.com/Jiawang1209/GeneFamilyFlow/issues).

## Installation

### Conda environment creation fails

**Symptom**: `conda env create -f envs/genefamily.yaml` hangs or fails on dependency resolution.

**Fix**: Use `mamba` instead of `conda` — much faster dependency solver:

```bash
conda install -n base -c conda-forge mamba
mamba env create -f envs/genefamily.yaml
```

Or use the minimal environment for development/testing:

```bash
conda env create -f envs/genefamily-minimal.yaml
```

### `pfam_scan.pl` not found

**Symptom**: `step04_pfam_scan` fails with `command not found`.

**Fix**: Install via bioconda:

```bash
mamba install -n GeneFamilyFlow -c bioconda pfam_scan
```

### MCScanX compilation errors

**Symptom**: `make` fails with C++ errors.

**Fix**: MCScanX requires old C++ standards. Edit `MCScanX/Makefile`:

```makefile
CXXFLAGS = -O3 -fcommon
```

Then:

```bash
cd MCScanX && make
```

Update `config/default_config.yaml`:

```yaml
step9:
  mcscanx_path: "/absolute/path/to/MCScanX/MCScanX"
```

---

## Pfam Database

### Where to download

```bash
mkdir -p ~/pfam_db && cd ~/pfam_db
wget http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
wget http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz
gunzip *.gz
hmmpress Pfam-A.hmm
```

Update config:

```yaml
step4:
  pfam_database_dir: "/Users/YOUR_USER/pfam_db"
```

### `hmmpress` must be run before first use

**Symptom**: `Error: failed to open binary auxfiles`.

**Fix**:

```bash
cd ~/pfam_db && hmmpress Pfam-A.hmm
```

---

## Runtime Errors

### Snakemake: "MissingInputException"

**Symptom**:

```
MissingInputException: Missing input files for rule step02_hmmsearch_pfam:
  example/2.hmmsearch/PF00854.hmm
```

**Fix**: Check `step2.domains[].pfam_hmm_file` paths in your config. Paths are relative to the working directory where you run `snakemake`.

### BLAST: "BLAST Database error"

**Symptom**: `BLAST Database error: No alias or index file found`.

**Fix**: Delete stale BLAST indexes and rerun:

```bash
rm -rf work/03_blast/
snakemake --configfile config/default_config.yaml -j 8
```

### IQ-TREE out of memory

**Symptom**: `step06_tree_build` killed by OOM.

**Fix**: Either reduce bootstrap replicates or switch to FastTree:

```yaml
step6:
  tree_tool: "fasttree"
  fasttree_options: "-gamma"
```

FastTree uses ~10x less memory than IQ-TREE.

### R script fails: `package not found`

**Symptom**: `Error in library(ggtree): there is no package called 'ggtree'`.

**Fix**: Install missing R packages inside the conda env:

```bash
conda activate GeneFamilyFlow
R -e 'BiocManager::install(c("ggtree", "Biostrings"))'
R -e 'install.packages(c("ggplot2", "dplyr", "optparse", "openxlsx"))'
```

---

## External Tools

### KaKs_Calculator not found

**Fix**: Download from [SourceForge](https://sourceforge.net/projects/kakscalculator2/), compile, and add to `PATH`:

```bash
export PATH=$PATH:/path/to/KaKs_Calculator2/bin
```

Or set `step8.precomputed: true` and provide `kaks_file` to skip the computation.

### PlantCARE (online tool)

PlantCARE is a web service and cannot be automated. Workflow:

1. Run pipeline until Step 10 produces `promoter_sequences.fa`.
2. Submit sequences manually at [PlantCARE](https://bioinformatics.psb.ugent.be/webtools/plantcare/html/).
3. Download results into `step10.plantcare_dir`.
4. Rerun: `snakemake --until step10_promoter`.

### `ggNetView` R package

`ggNetView` is not on CRAN. Install manually:

```bash
R -e 'install.packages("ggNetView", repos="http://R-Forge.R-project.org")'
```

---

## Performance

### Pipeline is slow

Check:

1. **Parallel jobs**: Use `-j 8` or higher.
2. **Threads per rule**: Increase `max_threads` in config.
3. **Tree tool**: Switch from `iqtree` to `fasttree` (10-100x faster).
4. **Precomputed mode**: Set `step8.precomputed: true` and `step9.precomputed: true` if you have results.

### Disk space exhausted

**Fix**: Clean intermediate files:

```bash
rm -rf work/
```

Or set:

```yaml
keep_intermediate: false
compress_intermediate: true
```

---

## Dry-run and Debugging

### Inspect the DAG without running

```bash
snakemake --configfile config/default_config.yaml -n -p
```

### Visualize the DAG

```bash
snakemake --configfile config/default_config.yaml --dag | dot -Tpdf > dag.pdf
```

### Force re-run a specific step

```bash
snakemake --configfile config/default_config.yaml --forcerun step06_tree_build -j 4
```

### Unlock a stale working directory

```bash
snakemake --configfile config/default_config.yaml --unlock
```
