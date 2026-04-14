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
└── 11_ppi/
```

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

| File | Format | Description |
|------|--------|-------------|
| `identify.ID.clean.fa` | FASTA | Final curated gene family protein sequences across all species. |
| `identify.ID.txt` | TXT | Gene IDs that passed HMM + BLAST + Pfam verification. |
| `per_domain/{domain}.ID.txt` | TXT | Per-domain hit IDs before merging. |

**Interpretation**: `identify.ID.clean.fa` is the input for all downstream analyses. Verify the count matches known family size for your species.

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
| `phylogenetic_tree.pdf` | PDF | Rendered tree figure with bootstrap values. |
| `alignment.aln` | FASTA | Multiple sequence alignment (MUSCLE or ClustalW). |
| `tree.nwk` | Newick | Raw tree file (IQ-TREE or FastTree). |
| `tree.contree` | Newick | Consensus tree with bootstrap support (IQ-TREE only). |
| `iqtree.log` | TXT | IQ-TREE run log with model selection details. |

**Interpretation**:
- Bootstrap values ≥ 80 = high confidence
- Use the tree to manually curate subfamilies in `tree_group_2.xlsx`

---

## Step 7: Motif and Gene Structure

| File | Format | Description |
|------|--------|-------------|
| `meme_location.txt` | TXT | Motif coordinates in each protein. |
| `meme.html`, `meme.xml` | MEME native output | Motif logos and statistics. |
| `Tree_Domain_Motif_GeneStructure.pdf` | PDF | **Composite figure**: phylogenetic tree + domain architecture + motif positions + exon/intron structure, all aligned by gene. |

**Interpretation**: The composite figure is the main result — it visualizes conservation of domains, motifs, and gene structure across the family.

---

## Step 8: JCVI Synteny and Ka/Ks

| File | Format | Description |
|------|--------|-------------|
| `08_jcvi_kaks.pdf` | PDF | Ka/Ks distribution (violin/box plot) per species pair. |
| `kaks.tab.xls` | TSV | Raw Ka/Ks values for every ortholog pair. |
| `synteny.pdf` | PDF | JCVI synteny dotplot (when `precomputed: false`). |

**Interpretation**:
- **Ka/Ks < 1** = purifying (negative) selection
- **Ka/Ks ≈ 1** = neutral evolution
- **Ka/Ks > 1** = positive (diversifying) selection
- Most gene families show Ka/Ks < 0.5

---

## Step 9: MCScanX Synteny and Circos

| File | Format | Description |
|------|--------|-------------|
| `circos_{species}.pdf` | PDF | Circos plot per species showing chromosome structure, gene family members, and intra-genome synteny. |
| `{species}.collinearity` | TXT | MCScanX collinear block file. |
| `{species}.tandem` | TXT | Tandem duplicate gene pairs. |

**Interpretation**:
- Red lines in Circos = family gene duplication pairs
- Tandem clusters indicate recent local expansion
- Whole-genome segmental duplications show as long-range links

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
| `PPI_network.pdf` | PDF | Protein-protein interaction network figure. |
| `edge_annotation.xlsx` | XLSX | Edge list with interaction scores and annotations. |
| `node_annotation.tsv` | TSV | Node metadata (Pfam domains, subcellular localization). |

**Interpretation**:
- Nodes = family members + interaction partners
- Edge thickness = interaction confidence (AraNet score)
- Hub nodes (high degree) = central family members

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
