#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(ggplot2)
})

option_list <- list(make_option("--outdir", type = "character"))
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

motif_file <- file.path(opt$outdir, "meme_location.txt")
if (!file.exists(motif_file)) {
  quit(status = 0)
}

motif_df <- read_delim(motif_file, delim = "\t", col_names = c("Gene", "Start", "End", "Motif"), show_col_types = FALSE)
if (nrow(motif_df) == 0) {
  quit(status = 0)
}

p <- ggplot(motif_df, aes(x = Start, xend = End, y = Gene, yend = Gene, color = Motif)) +
  geom_segment(linewidth = 2) +
  theme_bw() +
  labs(title = "Motif Structure", x = "Position", y = "Gene")

ggsave(file.path(opt$outdir, "Tree_Domain_Motif_GeneStructure.pdf"), p, width = 12, height = 10)
