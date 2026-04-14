# Output File Reference

Every file produced by GeneFamilyFlow and how to interpret it. Organized by pipeline step.

```
output/
├── 02_hmmsearch/
├── 03_blast/
├── 04_identification/
├── 05_genefamily_info/
├── 06_tree/
├── 07_motif/
├── 08_jcvi_kaks/
├── 09_circos/
├── 10_promoter/
├── 11_ppi/
├── 13_go_kegg/      # only when step13_go_kegg.enabled: true
└── 14_qrt_pcr/      # only when step14_qrtpcr.enabled: true
```

Only Step 4's final family FASTA, plus the plotted PDFs from later steps, are
promoted into `output/`. Step 2/3 intermediate ID lists, the tree Newick, MEME
native outputs, JCVI Ka/Ks tables, and MCScanX collinearity files all live under
`work/` or (for precomputed steps 8 and 9) directly under `example/`. See the
per-step tables below for the canonical locations.

---

## Step 2: HMM Search

| File | Format | Description |
|------|--------|-------------|
| `{domain}/{species}.pfam.domtblout` | HMMER domtblout | First-round search: species proteins vs. canonical Pfam HMM. |
| `{domain}/hmm_candidate_ids.txt` | TXT | Union of first-round hits across species (used to build custom HMM). |
| `{domain}/1st_round.aln` | ClustalW | Alignment of first-round hits. |
| `{domain}/custom.hmm` | HMMER | Custom HMM built from first-round hits (captures family-specific signal). |
| `{domain}/{species}.custom.domtblout` | HMMER domtblout | Second-round search with the custom HMM. |
| `{domain}/hmm_2nd_ids.txt` | TXT | **Final HMM candidate IDs** for the domain. |

**Interpretation**: The two-round approach catches distant family members that the canonical Pfam HMM might miss. `hmm_2nd_ids.txt` is the HMM half of the gene family evidence (merged with BLAST in Step 4).

---

## Step 3: BLAST Search

| File | Format | Description |
|------|--------|-------------|
| `{domain}/seed.fa` | FASTA | Seed sequences extracted per domain from `step3.reference_species` by `rule step03_extract_seed`. |
| `{domain}/seed_db.*` | BLAST index | Indexed seed sequence database built from `seed.fa`. |
| `{domain}/{species}.blast` | TSV (BLAST `-outfmt 6`) | Raw BLASTP hits per species (with `qlen`/`slen` columns). |
| `{domain}/blast_ids.txt` | TXT | **Final BLAST candidate IDs** for the domain (E-value filtered). |

**Interpretation**: BLAST provides a sequence-similarity-based second opinion on HMM hits. `blast_ids.txt` is merged with `hmm_2nd_ids.txt` in Step 4 using `step4.merge_method` (`intersection` = strict, `union` = lenient).

---

## Step 4: Gene Family Identification

The only file promoted into `output/` is the final curated FASTA. Every
intermediate ID list lives under `work/04_identification/`:

| File | Format | Description |
|------|--------|-------------|
| `output/04_identification/identify.ID.clean.fa` | FASTA | **Final curated gene family protein sequences** across all species — the canonical input for every downstream step. |
| `work/04_identification/{domain}.merged.ID` | TXT | Per-domain merge of HMM ∩/∪ BLAST IDs (mode from `step4.merge_method`). |
| `work/04_identification/all_domains.ID` | TXT | Union of every per-domain merged ID list, deduped. |
| `work/04_identification/all_domains.ID.fa` | FASTA | Combined candidate FASTA fed into `pfam_scan.pl`. |
| `work/04_identification/Pfam_scan.out` | pfam_scan | Raw `pfam_scan.pl` verification output (only when `step4.use_pfam_scan: true`). |
| `work/04_identification/pfam_scan.id` | TXT | Gene IDs that satisfy `step4.domain_combine` against the Pfam verification. |
| `work/04_identification/{species}.family.ID` | TXT | Per-species family ∩ species proteome — consumed by steps 9, 10, and 13 via `family_ids_for()`. |

**Interpretation**: `identify.ID.clean.fa` is the input for all downstream analyses. Verify the count matches known family size for your species. The `{species}.family.ID` split is what lets downstream steps fan out per species without re-walking the combined FASTA.

---

## Step 5: Gene Family Statistics

| File | Format | Description |
|------|--------|-------------|
| `Gene_Information.xlsx` | XLSX | Gene properties table (human-readable). |
| `Gene_Information.csv` | CSV | Same data, programmatic access. |

**Columns**:
- `Gene_ID`, `Species`, `Chromosome`, `Start`, `End`, `Strand`
- `Protein_Length` (aa)
- `Molecular_Weight` (Da)
- `Isoelectric_Point` (pI)
- `GRAVY` (hydrophobicity)
- `Instability_Index`
- `Subcellular_Localization` (predicted)

---

## Step 6: Phylogenetic Tree

| File | Format | Description |
|------|--------|-------------|
| `output/06_tree/phylogenetic_tree_circular.pdf` | PDF | Circular ggtree layout with auto-split subfamily highlights and branch coloring. |
| `output/06_tree/phylogenetic_tree_rectangular.pdf` | PDF | Rectangular ggtree layout with the same subfamily annotations. |
| `work/06_tree/alignment.muscle` | FASTA | MUSCLE v5 multiple sequence alignment of the family. |
| `work/06_tree/tree.nwk` | Newick | Final tree file — copied from `iqtree_out.treefile` or FastTree stdout, depending on `step6.tree_tool`. |
| `work/06_tree/iqtree_out.*` | IQ-TREE native | Full IQ-TREE run artifacts (`.iqtree`, `.treefile`, `.contree`, `.log`, `.mldist`, etc.) when `tree_tool: iqtree`. |
| `logs/06_tree_build.log` | TXT | Rule-level log with tool stdout/stderr (model selection, bootstrap progress). |

**Interpretation**:
- Bootstrap values ≥ 80 = high confidence
- `step6.tree_layout` controls which PDFs are produced: `circular`, `rectangular`, or `both` (default).
- Subfamilies are auto-split top-down into `step6.auto_n_subfamilies` monophyletic clades; override with `step6.group_file` when you have a manual curation.

---

## Step 7: Motif and Gene Structure

| File | Format | Description |
|------|--------|-------------|
| `output/07_motif/meme_info.txt` | TXT | Per-motif summary parsed from MEME output (width, E-value, consensus, site count). |
| `output/07_motif/meme_location.txt` | TXT | Motif occurrence coordinates in each protein — consumed by `R/07_domain_motif_structure.R`. |
| `output/07_motif/Tree_Domain_Motif_GeneStructure.pdf` | PDF | **Composite figure**: phylogenetic tree + domain architecture + motif positions + exon/intron structure, all aligned by gene. |
| `work/07_motif/meme_out/meme.txt` | MEME native | Full MEME text output — input to `parse_meme_output.py`. |
| `work/07_motif/meme_out/meme.html`, `meme.xml`, `logo*.png` | MEME native | Motif logos, per-motif diagnostics, and interactive browser output (inspect when sanity-checking the discovered motifs). |

**Interpretation**: The composite figure is the main result — it visualizes conservation of domains, motifs, and gene structure across the family. Use `meme_info.txt` to decide which motifs are informative; use `work/07_motif/meme_out/meme.html` to inspect the raw logos if a motif looks suspicious.

---

## Step 8: JCVI Synteny and Ka/Ks

| File | Format | Description |
|------|--------|-------------|
| `output/08_jcvi_kaks/08_jcvi_kaks.pdf` | PDF | Ka/Ks distribution (violin/box plot) per species pair — the only artifact promoted to `output/`. |
| `<kaks_file>` (`example/8.collinearity/kaks.tab.xls` by default) | TSV | Raw Ka/Ks values for every syntenic ortholog pair. Path is controlled by `step8.kaks_file` when `precomputed: true`. |
| `work/08_jcvi/kaks.tab.xls` | TSV | Ka/Ks table when `precomputed: false` — assembled in place by `rule step08_kaks_calc` from per-pair `*.kaks` files. |
| `work/08_jcvi/{sp_a}.{sp_b}.anchors`, `*.anchors.simple` | JCVI | Pairwise anchor and simplified synteny files (only when `precomputed: false`). |

**Interpretation**:
- **Ka/Ks < 1** = purifying (negative) selection
- **Ka/Ks ≈ 1** = neutral evolution
- **Ka/Ks > 1** = positive (diversifying) selection
- Most gene families show Ka/Ks < 0.5
- There is no synteny dotplot rule in this pipeline — if you need one, generate it separately from the `work/08_jcvi/*.anchors` files.

---

## Step 9: MCScanX Synteny and Circos

| File | Format | Description |
|------|--------|-------------|
| `output/09_circos/circos_{species}.pdf` | PDF | Circos plot per species showing chromosome structure, gene family members, and intra-genome synteny. One PDF per entry in `config.species`. |
| `example/9.mcscanx/{species}.collinearity` *or* `work/09_mcscanx/{species}.collinearity` | MCScanX | Collinear block file — shipped under `example/` when `step9.precomputed: true`, or produced under `work/` by `rule step09_mcscanx`. |
| `example/9.mcscanx/{species}.tandem` *or* `work/09_mcscanx/{species}.tandem` | MCScanX | Tandem duplicate gene pairs — same precomputed/computed split as above. |
| `work/09_mcscanx/{species}.gene_type` | MCScanX | Duplicate-type classification from `duplicate_gene_classifier` (only produced when `precomputed: false`). |

**Interpretation**:
- Red lines in Circos = family gene duplication pairs
- Tandem clusters indicate recent local expansion
- Whole-genome segmental duplications show as long-range links
- Only `circos_{species}.pdf` lives under `output/`; the collinearity/tandem tables stay under `example/` or `work/` because `R/09_circos.R` consumes them directly via `family_ids_for()`.

---

## Step 10: Promoter Analysis

| File | Format | Description |
|------|--------|-------------|
| `output/10_promoter/promoter_elements.pdf` | PDF | Cis-element distribution heatmap across promoters. |
| `work/10_promoter/{target}_upstream.fasta` | FASTA | Extracted upstream regions (when `compute_promoter: true`). |
| `work/10_promoter/jaspar_fimo/plantCARE_output_jaspar.tab` | TSV | FIMO + JASPAR Plants hits in PlantCARE 8-column format (when `scan_method=jaspar`). |
| `work/10_promoter/local_plantcare/plantCARE_output_local.tab` | TSV | Literal-substring scan hits in PlantCARE 8-column format (when `scan_method=local`). |

**Interpretation**:
- Rows = family members, columns = cis-element categories
- Heatmap intensity = element count
- Common categories: light-responsive, hormone-responsive (ABA, ABRE, ERE), stress-responsive
- The intermediate `.tab` files are PlantCARE-compatible regardless of `scan_method`, so `R/10_promoter.R` consumes any of them unchanged.

---

## Step 11: PPI Network

| File | Format | Description |
|------|--------|-------------|
| `output/11_ppi/PPI_network.pdf` | PDF | Protein-protein interaction network figure rendered by `R/11_ppi.R`. |

The edge list consumed by this step (`step11.ppi_edge_file`, default
`example/11.ppi/edge_annotation.xlsx`) is an **input**, not an output — curate
it separately if you need annotations.

**Interpretation**:
- Nodes = family members + interaction partners
- Edge thickness = interaction confidence (AraNet score)
- Hub nodes (high degree) = central family members
- `step11.igraph_layout` controls the layout algorithm (`gephi`, `fr`, `circle`, `kk`).

---

## Step 13: GO/KEGG Enrichment (optional)

Enabled by setting `step13_go_kegg.enabled: true`. One PDF per species listed in
`step13_go_kegg.species`.

| File | Format | Description |
|------|--------|-------------|
| `output/13_go_kegg/{species}_enrichment.pdf` | PDF | Faceted GO + KEGG dotplot from `clusterProfiler::enricher`. |
| `work/13_go_kegg/{species}.go.tsv` | TSV | `TERM\tGENE` table for GO over-representation. |
| `work/13_go_kegg/{species}.kegg.tsv` | TSV | `TERM\tGENE` table for KEGG over-representation. |
| `work/13_go_kegg/{species}.universe.txt` | TXT | Background gene universe (every annotated gene in the emapper output). |
| `work/13_go_kegg/eggnog/{species}.emapper.annotations` | TSV | eggNOG-mapper output, only when `run_eggnog: true`. Otherwise the precomputed file under `step13_go_kegg.precomputed_dir/{species}.emapper.annotations` is used directly. |

**Interpretation**:
- Dot size = gene count in the family that maps to the term
- Dot color = adjusted p-value (thresholded by `step13_go_kegg.pvalue` / `qvalue`)
- Repo ships a curated offline fixture under `example/13.GO_KEGG/` so
  `enabled: true` + `run_eggnog: false` works out of the box; replace with real
  eggNOG-mapper output before drawing biological conclusions.

---

## Step 14: qRT-PCR (optional)

Enabled by setting `step14_qrtpcr.enabled: true` and providing
`expression_data`.

| File | Format | Description |
|------|--------|-------------|
| `output/14_qrt_pcr/qRT_PCR.pdf` | PDF | qRT-PCR expression figure across tissues/treatments. |
| `output/14_qrt_pcr/qRT_PCR_summary.csv` | CSV | Summary table — per-gene means, SD, and sample counts. |

---

## Intermediate Files (`work/`)

The `work/` directory contains per-step intermediate files (HMM searches, BLAST hits, alignments). Useful for debugging but safe to delete after a successful run:

```bash
# Clean work directory
rm -rf work/
```

Set `keep_intermediate: false` in config to auto-delete on success.

## Log Files (`logs/`)

Per-rule execution logs. Check here when a step fails:

```bash
ls logs/
# 01_filter_Sbicolor.log
# 02_hmmsearch_pfam_PF00854_Sbicolor.log
# ...
```
