#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Biostrings)
  library(Peptides)
  library(tidyverse)
  library(writexl)
})

option_list <- list(
  make_option("--input_fasta", type = "character"),
  make_option("--input_bed", type = "character"),
  make_option("--outdir", type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
pep.fa <- readAAStringSet(opt$input_fasta)
bed <- read_delim(opt$input_bed, col_names = FALSE, delim = "\t", show_col_types = FALSE) %>%
  set_names(c("Chr", "Start", "End", "ID", "Info", "Strand"))

df <- tibble(ID = names(pep.fa), seq = as.character(pep.fa)) %>%
  mutate(
    Length = Peptides::lengthpep(seq),
    MW = Peptides::mw(seq),
    hydrophobicity = Peptides::hydrophobicity(seq),
    pI = Peptides::pI(seq)
  ) %>%
  left_join(bed, by = "ID")

write_xlsx(df, path = file.path(opt$outdir, "Gene_Information.xlsx"))
write_csv(df, file.path(opt$outdir, "Gene_Information.csv"))
