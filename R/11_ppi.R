#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(igraph)
  library(ggraph)
  library(ggplot2)
})

option_list <- list(make_option("--outdir", type = "character"))
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

g <- make_ring(10)
p <- ggraph(g, layout = "circle") +
  geom_edge_link(alpha = 0.5) +
  geom_node_point(size = 3, color = "#d6604d") +
  theme_void() +
  ggtitle("PPI placeholder network")

ggsave(file.path(opt$outdir, "PPI1.pdf"), p, width = 8, height = 8)
