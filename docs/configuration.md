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
      seed_references:
        - Athaliana
  hmmsearch_evalue: 1e-10
  clustering_tool: "clustalw"
  hmmbuild_calibrate: true
```

**Multi-domain**: Add multiple entries to `domains`. Each domain is searched independently; hits are merged in Step 4.

| Key | Description |
|-----|-------------|
| `domains[].pfam_id` | Pfam accession (e.g. `PF00854`). |
| `domains[].pfam_hmm_file` | Path to the Pfam `.hmm` model file. |
| `domains[].seed_references` | List of reference species keys (see `step3.reference_species`). The pipeline filters each reference by this `pfam_id` at runtime to build the BLAST seed FASTA. Omit to use every declared reference. |
| `hmmsearch_evalue` | E-value cutoff for `hmmsearch`. |
| `clustering_tool` | `clustalw` or `muscle` — used to build the custom HMM from round-1 hits. |
| `hmmbuild_calibrate` | Run `hmmbuild` with calibration. |

---

## Step 3: BLAST Search

```yaml
step3:
  reference_species:
    Athaliana:
      domains_table: "example/3.blast/references/Athaliana.domains.tsv"
      proteome:      "example/3.blast/references/Athaliana.pep.fasta"
    Osativa:
      domains_table: "example/3.blast/references/Osativa.domains.tsv"
      proteome:      "example/3.blast/references/Osativa.pep.fasta"
  blast_evalue: 1e-10
  blast_num_threads: 10
  blast_num_alignments: 10
  blast_outfmt: "6 std qlen slen"
```

**Seed extraction**. For each domain, `rule step03_extract_seed` invokes `scripts/extract_seed_by_pfam.py` once per reference species, filters the clean 2-column domains TSV by the target `pfam_id`, pulls matching proteins from the gene-level proteome FASTA, and concatenates the per-reference outputs into `output/03_blast/{domain}/seed.fa`. Swapping the target domain requires no manual FASTA preparation — just update `step2.domains`.

| Key | Description |
|-----|-------------|
| `reference_species.<name>.domains_table` | 2-column TSV `gene_id<TAB>pfam_id`, one row per hit. Clean files for Arabidopsis/rice are committed under `example/3.blast/references/`. |
| `reference_species.<name>.proteome` | Gene-level protein FASTA. The first whitespace-delimited token after `>` must be the gene ID (isoform suffix already stripped). Large files are gitignored; drop your own under the same directory. |
| `blast_evalue` / `blast_num_threads` / `blast_num_alignments` / `blast_outfmt` | Standard BLAST+ parameters. `blast_outfmt: "6 std qlen slen"` adds query/subject length columns for downstream filtering. |

---

## Step 4: Gene Family Identification

```yaml
step4:
  merge_method: "intersection"
  domain_combine: "any"
  use_hmmscan: true
  use_pfam_scan: true
  pfam_scan_evalue: 1e-5
  pfam_database_dir: "/data/Liu/LiuYue/Evolution/Pfam"
```

| Key | Description |
|-----|-------------|
| `merge_method` | Per-domain HMM/BLAST merge: `intersection` (strict: HMM ∩ BLAST) or `union` (lenient: HMM ∪ BLAST). |
| `domain_combine` | Cross-domain combination strategy: `any` (union — keep genes with at least one listed domain) or `all` (intersection — require every listed domain). See below. |
| `use_hmmscan` | Run `hmmscan` against full Pfam database. |
| `use_pfam_scan` | Run `pfam_scan.pl` for domain architecture verification. |
| `pfam_scan_evalue` | E-value cutoff for domain verification. |
| `pfam_database_dir` | **Required** — local Pfam database directory (contains `Pfam-A.hmm`). |

**Important**: `pfam_database_dir` must point to a pre-downloaded Pfam database. See [troubleshooting.md](troubleshooting.md#pfam-database).

### Choosing `domain_combine`

A gene family may contain multiple Pfam domains. Use this parameter to express the biological relationship between them:

- **`any` (default)** — A gene is a family member if it has **at least one** of the listed domains. Use this for families with variable architecture, where members may contain different subsets of the domains. Example: NLR immune receptors (some have only NB-ARC, others combine NB-ARC + LRR + TIR).

- **`all`** — A gene is a family member only if it has **every** listed domain. Use this for strict multi-domain families where the architecture is conserved across all members. Example: enzymes that require both a binding domain and a catalytic domain.

When `step2.domains` lists only one domain, `any` and `all` produce the same result.

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

  # "jaspar" (default, recommended) | "local" | "plantcare"
  scan_method: "jaspar"

  # scan_method == "jaspar": offline FIMO + JASPAR 2024 CORE plants (~805 PFMs)
  jaspar_meme:    "example/10.promoter/JASPAR2024_plants.meme"
  jaspar_url:     "https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_plants_non-redundant_pfms_meme.txt"
  fimo_threshold: 1e-4

  # scan_method == "local": offline literal-substring fixture (~142 motifs)
  motif_library:  "example/10.promoter/plantcare_motifs.tsv"

  # scan_method == "plantcare": manual web upload, kept for back-compat
  plantcare_dir:  "example/10.promoter/PlantCARE_410_SB_plantCARE"

  element_annotation_file: "example/10.promoter/cir_element.desc.20240509.xlsx"
  jaspar_element_annotation_file: "example/10.promoter/jaspar_element.desc.xlsx"
```

| Key | Description |
|-----|-------------|
| `compute_promoter` | Set `true` to extract promoters from `{target_species}.genome.fasta` using `bedtools`. |
| `upstream_distance` | bp upstream of TSS. |
| `scan_method` | Cis-element scan strategy. `jaspar` runs FIMO offline against the shipped JASPAR Plants bundle (recommended). `local` runs an offline literal-substring scan against the small fixture motif library. `plantcare` reads a folder of manually-downloaded PlantCARE results. |
| `jaspar_meme` / `jaspar_url` | Cached MEME bundle path and refresh URL. Run `python scripts/fetch_jaspar_plants.py` to materialize or refresh; use `--from-file` for air-gapped clusters. |
| `fimo_threshold` | FIMO p-value threshold (default `1e-4`). |
| `motif_library` | TSV motif library used by `scan_method=local`. |
| `plantcare_dir` | Folder of manually-downloaded PlantCARE outputs (only used by `scan_method=plantcare`). |
| `element_annotation_file` | Cis-element description/category file consumed by `R/10_promoter.R` for `scan_method` in `{local, plantcare}`. |
| `jaspar_element_annotation_file` | TF-family keyed description xlsx used by `R/10_promoter.R` when `scan_method == "jaspar"`. Auto-generated from the MEME bundle via `scripts/build_jaspar_element_desc.py`; override to supply a curated file. |

**Default is fully offline.** FIMO ships with `envs/genefamily.yaml` (`meme>=5.5`) and the JASPAR bundle is committed to the repo, so the pipeline runs end-to-end without any web round-trip. `scan_method=plantcare` remains available for users who specifically need PlantCARE annotations from [bioinformatics.psb.ugent.be/webtools/plantcare](https://bioinformatics.psb.ugent.be/webtools/plantcare/html/).

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

### Step 13 — GO/KEGG enrichment

```yaml
step13_go_kegg:
  enabled: false
  species: [Athaliana, Osativa]
  run_eggnog: false
  eggnog_data_dir: "/path/to/eggnog_data"   # required when run_eggnog: true
  eggnog_cpu: 20
  precomputed_dir: "example/13.GO_KEGG"     # {species}.emapper.annotations lives here
  kegg_names: "example/13.GO_KEGG/kegg_pathway_names.tsv"
  pvalue: 0.05
  qvalue: 0.2
  min_size: 5
  top_n: 20
```

| Key | Description |
|-----|-------------|
| `enabled` | Toggle. When `false`, the rule is skipped entirely. |
| `species` | Subset of top-level `species` to enrich. One PDF per species. |
| `run_eggnog` | `true` — run `emapper.py` per species (requires `eggnog_data_dir`). `false` — read pre-computed annotations from `precomputed_dir`. |
| `precomputed_dir` | Directory with one `{species}.emapper.annotations` TSV per listed species. The repo ships a curated offline fixture here for `Athaliana` and `Osativa` so `enabled: true` + `run_eggnog: false` works out of the box. |
| `kegg_names` | Optional offline KEGG pathway name table (columns: `term`, `name`). Unknown `ko` IDs fall back to the raw code. |
| `pvalue` / `qvalue` | `clusterProfiler::enricher` thresholds. |
| `min_size` | Minimum gene-set size (`minGSSize`). |
| `top_n` | Top N terms per facet in the dotplot. |

> **Note on the bundled fixture**: the two shipped `.emapper.annotations` files are hand-curated for reproducible CI, not biological interpretation. Replace them with real eggNOG-mapper output before drawing conclusions.

### Step 14 — qRT-PCR

```yaml
step14_qrtpcr:
  enabled: false
  expression_data: "example/14.qRT_PCR/表达量_3.xlsx"
```

Both steps are disabled by default. Set `enabled: true` and provide input files to activate.

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
