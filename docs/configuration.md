# Configuration Reference

Every parameter in `config/default_config.yaml` explained. Copy the default file and edit it for your dataset:

```bash
cp config/default_config.yaml config/my_config.yaml
snakemake --configfile config/my_config.yaml -j 8
```

---

## General Settings

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `species` | list[str] | `[Sbicolor, Athaliana, Osativa]` | Species to analyze. Each name must match FASTA/GFF3 file prefixes in `database_dir`. |
| `target_species` | str | `"Sbicolor"` | Species used for single-species detailed analyses (tree, motif, promoter, PPI). |
| `database_dir` | path | `"example/1.database"` | Input directory containing `{species}.pep.fasta`, `{species}.gff3`, `{species}.cds.fasta`, `{species}.genome.fasta`. |
| `output_dir` | path | `"output"` | Final results directory. |
| `work_dir` | path | `"work"` | Intermediate files directory. |
| `log_dir` | path | `"logs"` | Per-rule log files. |

---

## Step 1: Data Preparation

```yaml
step1:
  filter_gff3: false
  gff3_gene_ids: []
  gff3_chromosomes: []
```

| Key | Description |
|-----|-------------|
| `filter_gff3` | Enable GFF3 subsetting (useful for large genomes). |
| `gff3_gene_ids` | Gene IDs to keep (empty = all). |
| `gff3_chromosomes` | Chromosomes to keep (empty = all). |

---

## Step 2: HMM Search

```yaml
step2:
  domains:
    - pfam_id: "PF00854"
      pfam_hmm_file: "example/2.hmmsearch/PF00854.hmm"
      seed_sequence_file: "example/3.blast/PF00854.TAIR.ID.fa"
  hmmsearch_evalue: 1e-10
  clustering_tool: "clustalw"
  hmmbuild_calibrate: true
```

**Multi-domain**: Add multiple entries to `domains`. Each domain is searched independently; hits are merged in Step 4.

| Key | Description |
|-----|-------------|
| `domains[].pfam_id` | Pfam accession (e.g. `PF00854`). |
| `domains[].pfam_hmm_file` | Path to the Pfam `.hmm` model file. |
| `domains[].seed_sequence_file` | Seed sequences from a reference species (used for BLAST in Step 3). |
| `hmmsearch_evalue` | E-value cutoff for `hmmsearch`. |
| `clustering_tool` | `clustalw` or `muscle` — used to build the custom HMM from round-1 hits. |
| `hmmbuild_calibrate` | Run `hmmbuild` with calibration. |

---

## Step 3: BLAST Search

```yaml
step3:
  blast_evalue: 1e-10
  blast_num_threads: 10
  blast_num_alignments: 10
  blast_outfmt: "6 std qlen slen"
```

Standard BLAST+ parameters. `blast_outfmt: "6 std qlen slen"` adds query/subject length columns to tabular output for downstream filtering.

---

## Step 4: Gene Family Identification

```yaml
step4:
  merge_method: "intersection"
  use_hmmscan: true
  use_pfam_scan: true
  pfam_scan_evalue: 1e-5
  pfam_database_dir: "/data/Liu/LiuYue/Evolution/Pfam"
```

| Key | Description |
|-----|-------------|
| `merge_method` | `intersection` (strict: HMM ∩ BLAST) or `union` (lenient: HMM ∪ BLAST). |
| `use_hmmscan` | Run `hmmscan` against full Pfam database. |
| `use_pfam_scan` | Run `pfam_scan.pl` for domain architecture verification. |
| `pfam_scan_evalue` | E-value cutoff for domain verification. |
| `pfam_database_dir` | **Required** — local Pfam database directory (contains `Pfam-A.hmm`). |

**Important**: `pfam_database_dir` must point to a pre-downloaded Pfam database. See [troubleshooting.md](troubleshooting.md#pfam-database).

---

## Step 5: Gene Family Statistics

```yaml
step5:
  gene_fasta_file: "example/4.identification/identify.ID.clean.fa"
  gene_bed_file: "example/10.promoter/species_10.bed"
  output_formats: [xlsx, csv]
  r_script: "R/05_genefamily_info.R"
```

Computes molecular weight, isoelectric point, hydrophobicity, subcellular localization prediction, and genomic coordinates.

---

## Step 6: Phylogenetic Tree

```yaml
step6:
  alignment_tool: "muscle"
  alignment_threads: 8
  tree_tool: "iqtree"
  iqtree_model: "MFP"
  iqtree_bootstrap: 1000
  iqtree_seed: 12345
  iqtree_threads: 8
  fasttree_options: "-gamma"
  tree_figure_width: 15
  tree_figure_height: 15
```

| Key | Description |
|-----|-------------|
| `alignment_tool` | `muscle` (recommended) or `clustalw`. |
| `tree_tool` | `iqtree` (accurate, slower) or `fasttree` (fast, approximate). |
| `iqtree_model` | `MFP` = automatic model selection; or specify e.g. `JTT+G`. |
| `iqtree_bootstrap` | Ultrafast bootstrap replicates. |
| `fasttree_options` | Extra flags when `tree_tool == fasttree`. |

**Note**: After this step, manual curation of `tree_group_2.xlsx` is required to define subfamilies used by downstream steps.

---

## Step 7: Motif and Gene Structure

```yaml
step7:
  meme_mod: "anr"
  meme_protein: true
  meme_nmotifs: 10
  meme_minw: 10
  meme_maxw: 50
  species_prefix: "Sobic."
```

MEME motif discovery parameters. `species_prefix` filters target-species genes for the composite figure.

---

## Step 8: JCVI Synteny and Ka/Ks

```yaml
step8:
  precomputed: true
  jcvi_database_type: "prot"
  jcvi_evalue: 1e-10
  synteny_minspan: 30
  kaks_file: "example/8.collinearity/kaks.tab.xls"
```

**`precomputed: true`** (default) — skip JCVI/Ka/Ks computation and use the file at `kaks_file`. Set to `false` to run from scratch (requires `{species}.gff3` and `{species}.cds.fasta`).

---

## Step 9: MCScanX Synteny and Circos

```yaml
step9:
  precomputed: true
  mcscanx_path: "~/software/MCScanX/MCScanX"
  circos_chr_style: "roundrect"
  species_config:
    Sbicolor:
      n_chromosomes: 10
      genome_length_file: "example/9.mcscanx/Sbicolor_454_v3.0.1.fa.fai"
      npf_id_file: "example/9.mcscanx/Sbicolor.NPF.id"
```

| Key | Description |
|-----|-------------|
| `precomputed` | Use pre-computed synteny files (set `false` to run MCScanX). |
| `mcscanx_path` | Path to `MCScanX` binary (requires manual install). |
| `species_config` | Per-species chromosome count, `.fai` genome index, and gene ID file. |

---

## Step 10: Promoter Analysis

```yaml
step10:
  compute_promoter: false
  upstream_distance: 2000
  plantcare_dir: "example/10.promoter/PlantCARE_410_SB_plantCARE"
  element_annotation_file: "example/10.promoter/cir_element.desc.20240509.xlsx"
```

| Key | Description |
|-----|-------------|
| `compute_promoter` | Set `true` to extract promoters from `{target_species}.genome.fasta` using `bedtools`. |
| `upstream_distance` | bp upstream of TSS. |
| `plantcare_dir` | Directory of PlantCARE output files (manual submission required). |
| `element_annotation_file` | Cis-element description/category file. |

**PlantCARE is an online tool** — submit sequences at [bioinformatics.psb.ugent.be/webtools/plantcare](https://bioinformatics.psb.ugent.be/webtools/plantcare/html/) and download results into `plantcare_dir`.

---

## Step 11: PPI Network

```yaml
step11:
  ppi_database_file: "example/11.ppi/AraNet.txt"
  ppi_edge_file: "example/11.ppi/edge_annotation.xlsx"
  igraph_layout: "gephi"
```

`igraph_layout` options: `gephi`, `fr` (Fruchterman-Reingold), `circle`, `kk` (Kamada-Kawai).

---

## Optional Steps (13-14)

```yaml
step13_go_kegg:
  enabled: false
  annotation_file: "example/13.GO_KEGG/annotations.tsv"

step14_qrtpcr:
  enabled: false
  expression_data: "example/14.qRT_PCR/qrt_pcr_data.xlsx"
```

Disabled by default. Set `enabled: true` and provide input files to activate.

---

## Runtime Settings

```yaml
max_threads: 8
snakemake_jobs: 8
use_conda: true
conda_env: "envs/genefamily.yaml"
keep_intermediate: true
```

| Key | Description |
|-----|-------------|
| `max_threads` | Max threads per rule. |
| `snakemake_jobs` | Concurrent jobs (equivalent to `-j` flag). |
| `use_conda` | Use conda environments (pass `--use-conda` to Snakemake). |
| `keep_intermediate` | Keep `work/` files after success. |
