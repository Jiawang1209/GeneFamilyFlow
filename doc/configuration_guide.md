# Configuration Guide

GeneFamilyFlow uses a YAML configuration file. Start from `config/default_config.yaml`.

## Top-level Keys

- `project_name`: project label
- `output_dir`: final outputs directory
- `work_dir`: intermediate working directory
- `log_dir`: logging directory
- `species`: species list
- `reference_species`: BLAST reference set
- `pfam`, `blast`, `tree`, `motif`, `promoter`: parameter blocks
- `tools`: executable names/paths

## Species Block

Each species must define:

- `name`
- `pep_fasta`
- `cds_fasta`
- `gff3`
- optional `genome_fasta`

## Reference Species

- `name`: one species in `species`
- `family_member_ids`: text file containing known family gene IDs

## Pfam Block

- `pfam_id`: e.g. `PF00854`
- `hmm_file`: local HMM file path
- `hmmscan_evalue`
- `hmmsearch_evalue`

## BLAST Block

- `evalue`
- `num_threads`
- `num_alignments`

## Tree Block

- `aligner`: usually `muscle`
- `iqtree_model`: usually `MFP`
- `bootstrap`: e.g. `1000`

## Motif Block

- `nmotifs`
- `minw`
- `maxw`

## Promoter Block

- `upstream_length`: e.g. `2000`

## Tools Block

Defines command names or absolute paths for external tools:

- `hmmsearch`, `hmmbuild`, `hmmscan`
- `blastp`, `makeblastdb`
- `seqkit`, `muscle`, `iqtree`, `meme`
- `Rscript`

## Recommended Workflow

1. Copy default config: `genefamilyflow init --output config.yaml`
2. Fill all species and reference paths.
3. Confirm tool executables are available in environment.
4. Run single step first for validation:
   - `genefamilyflow run-step step01_database --config config.yaml`
5. Run full pipeline after validation.
