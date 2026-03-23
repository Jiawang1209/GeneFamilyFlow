# GeneFamilyFlow Pipeline Steps

## Overview

The framework follows 11 Python-orchestrated steps plus optional R-only analyses.

## Step01 `step01_database`

- Standardize species input files into `work/01_database/{species}`.
- Required per species: `pep_fasta`, `cds_fasta`, `gff3`.

## Step02 `step02_hmmsearch`

- Run first-round HMM search using configured Pfam HMM.
- Build custom HMM from first-pass hits.
- Run second-round HMM search and produce `2st_id`.

## Step03 `step03_blast`

- Build BLAST database from reference-family proteins.
- BLAST all species proteins against this DB.
- Output merged candidate IDs.

## Step04 `step04_identification`

- Intersect HMM and BLAST candidates.
- Validate by `hmmscan` for target Pfam domain.
- Produce final `identify.ID` and FASTA.

## Step05 `step05_genefamily_info`

- Merge GFF and convert gene entries to BED.
- Prepare input FASTA for family members.
- Run `R/05_genefamily_info.R` for physicochemical statistics.

## Step06 `step06_tree`

- Multiple sequence alignment by MUSCLE.
- Phylogenetic inference by IQ-TREE.
- Run `R/06_tree.R` for tree plotting.

## Step07 `step07_motif`

- Run MEME motif analysis.
- Parse motif positions into tabular outputs.
- Run `R/07_domain_motif_structure.R` for motif-level visualization.

## Step08 `step08_collinearity`

- JCVI-oriented collinearity integration point.
- Run `R/08_jcvi_kaks.R` for Ka/Ks plotting.

## Step09 `step09_mcscanx`

- MCScanX-oriented duplicate/synteny integration point.
- Run `R/09_circos.R` for circos-style visualization.

## Step10 `step10_promoter`

- Promoter analysis integration point.
- Run `R/10_promoter.R` for cis-element plotting.

## Step11 `step11_ppi`

- PPI analysis integration point.
- Run `R/11_ppi.R` for network plotting.

## Optional R-only Analyses

- `R/12_rnaseq.R`
- `R/13_go_kegg.R`
- `R/14_qrt_pcr.R`
