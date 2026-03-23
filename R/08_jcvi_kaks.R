#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(ggplot2)
})

option_list <- list(make_option("--outdir", type = "character"))
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

outfile <- file.path(opt$outdir, "8.jcvi_KaKs.pdf")
pdf(outfile, width = 10, height = 6)
plot.new()
title("JCVI Ka/Ks placeholder plot")
dev.off()
