#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
})

option_list <- list(make_option("--outdir", type = "character"))
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

pdf(file.path(opt$outdir, "promoter_placeholder.pdf"), width = 10, height = 6)
plot.new()
title("Promoter cis-element analysis placeholder")
dev.off()
