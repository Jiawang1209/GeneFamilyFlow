#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(ape)
  library(ggplot2)
})

option_list <- list(
  make_option("--treefile", type = "character"),
  make_option("--outdir", type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

tree <- read.tree(opt$treefile)
pdf(file.path(opt$outdir, "phylogenetic_tree.pdf"), width = 10, height = 10)
plot(tree, cex = 0.5, main = "GeneFamilyFlow Phylogenetic Tree")
dev.off()
