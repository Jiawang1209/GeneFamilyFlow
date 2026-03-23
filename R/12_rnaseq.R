#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(pheatmap)
})

set.seed(1)
m <- matrix(rnorm(100), nrow = 10)
rownames(m) <- paste0("Gene", seq_len(10))
colnames(m) <- paste0("Sample", seq_len(10))

pdf("12_rnaseq_placeholder.pdf", width = 8, height = 10)
pheatmap(m, scale = "row")
dev.off()
