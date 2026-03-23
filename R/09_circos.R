#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(circlize)
})

option_list <- list(make_option("--outdir", type = "character"))
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

pdf(file.path(opt$outdir, "circos_placeholder.pdf"), width = 8, height = 8)
circos.clear()
circos.initialize(factors = "Chr1", xlim = c(0, 1))
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  circos.text(0.5, 0.5, "MCScanX/Circos placeholder")
})
circos.clear()
dev.off()
