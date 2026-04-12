#!/usr/bin/env Rscript
# 11_ppi.R — Protein-protein interaction network visualization (ggNetView)
#
# Inputs:
#   --edge_file       Tab/CSV/XLSX with columns: from, to, weight (or Source, Target, Weight)
#   --node_annotation Optional tab/CSV/XLSX with node attributes (first column = node ID)
#   --outdir          Output directory
#   --layout          Network layout: "gephi", "fr", "circle", "kk" (default "gephi")
#   --top_modules     Max community modules to highlight (default 15)
#   --module_method   Module detection: "Fast_greedy", "Louvain", "Walktrap" (default "Fast_greedy")
#   --node_add        Node size scaling factor (default 7)
#   --shrink          Module layout shrink factor (default 0.75)
#   --seed            Random seed for layout reproducibility
#   --width, --height

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(ggNetView)
  library(readxl)
  library(patchwork)
})

option_list <- list(
  make_option("--edge_file",       type = "character"),
  make_option("--node_annotation", type = "character", default = NULL),
  make_option("--outdir",          type = "character", default = "output/11_ppi"),
  make_option("--layout",          type = "character", default = "gephi"),
  make_option("--top_modules",     type = "integer",   default = 15),
  make_option("--module_method",   type = "character", default = "Fast_greedy"),
  make_option("--node_add",        type = "double",    default = 7),
  make_option("--shrink",          type = "double",    default = 0.75),
  make_option("--seed",            type = "integer",   default = 1115),
  make_option("--output_name",     type = "character", default = "PPI_network.pdf"),
  make_option("--width",           type = "double",    default = 10),
  make_option("--height",          type = "double",    default = 10)
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$edge_file)) stop("--edge_file is required")
if (!file.exists(opt$edge_file)) stop(sprintf("edge_file not found: %s", opt$edge_file))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# ---- Read edge file (auto-detect format) ---------------------------------
read_auto <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "xlsx") {
    read_xlsx(path)
  } else if (ext == "csv") {
    read_csv(path, show_col_types = FALSE)
  } else {
    read_delim(path, delim = "\t", show_col_types = FALSE)
  }
}

edge_raw <- read_auto(opt$edge_file)

# Normalize column names
col_lower <- tolower(colnames(edge_raw))
rename_map <- c(
  from   = "from",   to     = "to",     weight = "weight",
  source = "from",   target = "to"
)
for (orig in names(rename_map)) {
  idx <- which(col_lower == orig)
  if (length(idx) == 1) colnames(edge_raw)[idx] <- rename_map[orig]
}

required_cols <- c("from", "to")
missing <- setdiff(required_cols, colnames(edge_raw))
if (length(missing) > 0) {
  stop(sprintf("Edge file missing columns: %s", paste(missing, collapse = ", ")))
}

edge <- edge_raw %>%
  select(from, to, any_of("weight")) %>%
  filter(from != to) %>%
  distinct()

if (!"weight" %in% colnames(edge)) {
  edge$weight <- 1
}

if (nrow(edge) == 0) stop("No edges after filtering")
message(sprintf("Loaded %d edges, %d unique nodes",
                nrow(edge), n_distinct(c(edge$from, edge$to))))

# ---- Read optional node annotation ---------------------------------------
node_ann <- NULL
if (!is.null(opt$node_annotation) && file.exists(opt$node_annotation)) {
  node_ann <- read_auto(opt$node_annotation)
  message(sprintf("Node annotation: %d rows, columns: %s",
                  nrow(node_ann), paste(colnames(node_ann), collapse = ", ")))
}

# ---- Build graph ---------------------------------------------------------
graph_obj <- build_graph_from_df(
  df             = as.data.frame(edge),
  node_annotation = if (!is.null(node_ann)) as.data.frame(node_ann) else NULL,
  directed       = FALSE,
  top_modules    = opt$top_modules,
  module.method  = opt$module_method,
  seed           = opt$seed
)

# ---- Plot ----------------------------------------------------------------
p <- ggNetView(
  graph_obj    = graph_obj,
  layout       = opt$layout,
  node_add     = opt$node_add,
  center       = FALSE,
  layout.module = "random",
  shrink       = opt$shrink
)

out_pdf <- file.path(opt$outdir, opt$output_name)
ggsave(out_pdf, p, width = opt$width, height = opt$height)

# Write tidy edge/node tables
write_csv(edge, file.path(opt$outdir, "ppi_edges.csv"))
node_ids <- tibble(node = sort(unique(c(edge$from, edge$to))))
if (!is.null(node_ann)) {
  node_ids <- node_ids %>%
    left_join(node_ann, by = setNames(colnames(node_ann)[1], "node"))
}
write_csv(node_ids, file.path(opt$outdir, "ppi_nodes.csv"))

message(sprintf("Wrote %s | edges=%d  nodes=%d  modules=%d",
                out_pdf, nrow(edge), nrow(node_ids), opt$top_modules))
