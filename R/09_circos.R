#!/usr/bin/env Rscript
# 09_circos.R — Generalized MCScanX + Circos visualization for one species
#
# Inputs (all required unless noted):
#   --chr_fai         FASTA index file (samtools faidx): chr\tlength\t...
#   --bed             Gene BED file (chr, start, end, id, ...)
#   --gene_ids        Text file listing family gene IDs (one per line)
#   --gene_type       MCScanX .gene_type output (id\ttype{0..4})
#   --tandem          MCScanX .tandem output (comma-separated pairs)
#   --collinearity    MCScanX .collinearity output
#   --outdir          Output directory
#   --n_chr           Number of chromosomes to plot (default 10)
#   --chr_prefix      Optional string to prepend to chr names (default "")
#   --id_regex_strip  Optional regex to strip from gene IDs before matching (default "\\.\\d+\\.p$|\\.\\d+$")
#   --window_size     Sliding window size in bp (default 500000)
#   --output_name     PDF filename (default circos.pdf)
#   --width, --height

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(circlize)
  library(ComplexHeatmap)
  library(grid)
})

option_list <- list(
  make_option("--chr_fai",        type = "character"),
  make_option("--bed",            type = "character"),
  make_option("--gene_ids",       type = "character"),
  make_option("--gene_type",      type = "character"),
  make_option("--tandem",         type = "character"),
  make_option("--collinearity",   type = "character"),
  make_option("--outdir",         type = "character", default = "output/09_circos"),
  make_option("--n_chr",          type = "integer",   default = 10),
  make_option("--chr_prefix",     type = "character", default = ""),
  make_option("--id_regex_strip", type = "character",
              default = "\\.\\d+\\.p$|\\.v\\d+\\.\\d+$|\\.\\d+$"),
  make_option("--window_size",    type = "double",    default = 500000),
  make_option("--output_name",    type = "character", default = "circos.pdf"),
  make_option("--width",          type = "double",    default = 8),
  make_option("--height",         type = "double",    default = 8)
)
opt <- parse_args(OptionParser(option_list = option_list))

required <- c("chr_fai", "bed", "gene_ids", "gene_type", "tandem", "collinearity")
for (k in required) {
  if (is.null(opt[[k]])) stop(sprintf("--%s is required", k))
  if (!file.exists(opt[[k]])) stop(sprintf("%s not found: %s", k, opt[[k]]))
}
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

strip_suffix <- function(x) str_remove(x, opt$id_regex_strip)
add_prefix   <- function(x) if (nzchar(opt$chr_prefix)) paste0(opt$chr_prefix, x) else x

# ---- Chromosomes ---------------------------------------------------------
chr_df <- read_delim(opt$chr_fai, delim = "\t", col_names = FALSE,
                     show_col_types = FALSE) %>%
  select(Chr = 1, Length = 2) %>%
  slice(seq_len(opt$n_chr)) %>%
  mutate(Chr = add_prefix(Chr), Start = 1) %>%
  select(Chr, Start, Length)

# ---- BED -----------------------------------------------------------------
bed_df <- read_delim(opt$bed, delim = "\t", col_names = FALSE,
                     show_col_types = FALSE) %>%
  select(Chr = 1, Start = 2, End = 3, ID = 4) %>%
  mutate(Chr = add_prefix(Chr))

# ---- Family gene IDs -----------------------------------------------------
family_ids <- read_lines(opt$gene_ids) %>%
  str_trim() %>%
  discard(~ . == "") %>%
  strip_suffix() %>%
  unique()

family_bed <- bed_df %>%
  mutate(ID_strip = strip_suffix(ID)) %>%
  filter(ID_strip %in% family_ids) %>%
  select(Chr, Start, End, ID = ID_strip) %>%
  distinct(ID, .keep_all = TRUE) %>%
  filter(Chr %in% chr_df$Chr)

if (nrow(family_bed) == 0) {
  warning("No family genes matched to BED; chr_prefix or id_regex_strip may need adjustment.")
}

# ---- Gene type -----------------------------------------------------------
genetype_df <- read_delim(opt$gene_type, delim = "\t", col_names = FALSE,
                          show_col_types = FALSE) %>%
  transmute(ID = strip_suffix(X1), Type = as.integer(X2)) %>%
  distinct(ID, .keep_all = TRUE)

gene_type_bed <- family_bed %>%
  left_join(genetype_df, by = "ID") %>%
  replace_na(list(Type = 0L)) %>%
  select(Chr, Start, End, ID, Type)

# ---- Tandem --------------------------------------------------------------
tandem_df <- tryCatch(
  read_delim(opt$tandem, delim = ",", col_names = FALSE, show_col_types = FALSE),
  error = function(e) tibble(X1 = character(), X2 = character())
) %>%
  mutate(X1 = strip_suffix(X1), X2 = strip_suffix(X2)) %>%
  filter(X1 %in% family_bed$ID | X2 %in% family_bed$ID)

# ---- Collinearity --------------------------------------------------------
col_raw <- read_lines(opt$collinearity) %>% discard(~ str_starts(., "#"))
col_df <- tibble(line = col_raw) %>%
  mutate(line = str_replace(line, "^\\s*\\d+-\\s*\\d+:\\s*", "")) %>%
  separate(line, into = c("X1", "X2", "rest"), sep = "\\s+", extra = "merge", fill = "right") %>%
  transmute(X1 = strip_suffix(X1), X2 = strip_suffix(X2)) %>%
  filter(X1 %in% family_bed$ID | X2 %in% family_bed$ID)

# Write gene-pair table for downstream inspection
bind_rows(
  tandem_df %>% mutate(Group = "Tandem"),
  col_df    %>% mutate(Group = "Collinearity")
) %>%
  write_csv(file.path(opt$outdir, "gene_pairs.csv"))

# ---- Sliding window counts -----------------------------------------------
make_windows <- function(chr_df, win) {
  chr_df %>%
    rowwise() %>%
    mutate(starts = list(seq(0, Length, by = win))) %>%
    unnest(starts) %>%
    transmute(Chr, Start = starts, End = pmin(starts + win, Length))
}
win_df <- make_windows(chr_df, opt$window_size)

count_in_windows <- function(windows, genes) {
  genes_by_chr <- split(genes, genes$Chr)
  windows %>%
    mutate(Number = map_int(seq_len(n()), function(i) {
      g <- genes_by_chr[[Chr[i]]]
      if (is.null(g)) return(0L)
      sum(g$Start >= Start[i] & g$End <= End[i])
    }))
}
win_df <- count_in_windows(win_df, family_bed)

# ---- Helper to look up BED position for a list of IDs --------------------
bed_lookup <- family_bed %>% select(ID, Chr, Start, End)

link_bed <- function(pairs) {
  pairs %>%
    inner_join(bed_lookup, by = c("X1" = "ID")) %>%
    rename(Chr1 = Chr, Start1 = Start, End1 = End) %>%
    inner_join(bed_lookup, by = c("X2" = "ID")) %>%
    rename(Chr2 = Chr, Start2 = Start, End2 = End)
}
tandem_link <- link_bed(tandem_df)
col_link    <- link_bed(col_df)

# ---- Plot ----------------------------------------------------------------
pdf_path <- file.path(opt$outdir, opt$output_name)
pdf(pdf_path, width = opt$width, height = opt$height)

circos.clear()
circos.genomicInitialize(
  as.data.frame(chr_df), plotType = NULL,
  axis.labels.cex = 0.4 * par("cex"),
  labels.cex = 0.6 * par("cex"),
  track.height = 0.01,
  major.by = 20000000
)

# Chromosome labels / ideogram
circos.genomicTrackPlotRegion(
  as.data.frame(chr_df), track.height = 0.05, stack = TRUE, bg.border = NA,
  panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value, col = "#dd3497", border = "black", ...)
  }
)

# Gene labels
if (nrow(family_bed) > 0) {
  circos.genomicLabels(
    as.data.frame(family_bed),
    labels.column = 4,
    padding = 0.1,
    connection_height = mm_h(3),
    col = as.numeric(factor(family_bed$Chr)),
    line_col = as.numeric(factor(family_bed$Chr)),
    cex = 0.7,
    side = "outside"
  )
}

# Sliding window density (set ylim explicitly to avoid degenerate-data error)
win_ymax <- max(win_df$Number, 1L)
circos.genomicTrack(
  as.data.frame(win_df %>% select(Chr, Start, End, Number)),
  track.height = 0.1, bg.col = "#f0f0f0", bg.border = NA,
  ylim = c(0, win_ymax),
  panel.fun = function(region, value, ...) {
    circos.genomicLines(region, value, col = "blue", lwd = 0.35, ...)
    circos.yaxis(labels.cex = 0.2, lwd = 0.1,
                 tick.length = convert_x(0.15, "mm"))
  }
)

# Gene type track
color_assign <- colorRamp2(
  breaks = c(0, 1, 2, 3, 4),
  col = c("#00ADFF", "#e66101", "#fdb863", "#b2abd2", "#5e3c99")
)
if (nrow(gene_type_bed) > 0) {
  circos.genomicTrackPlotRegion(
    as.data.frame(gene_type_bed %>% select(Chr, Start, End, Type)),
    track.height = 0.1, stack = TRUE, bg.border = NA,
    panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value,
                         col = color_assign(value[[1]]),
                         border = color_assign(value[[1]]), ...)
    }
  )
}

# Links
link_color <- scales::alpha("black", alpha = 1)
if (nrow(col_link) > 0) {
  circos.genomicLink(
    col_link %>% select(Chr1, Start1, End1) %>% as.data.frame(),
    col_link %>% select(Chr2, Start2, End2) %>% as.data.frame(),
    col = link_color, border = "#66c2a4", lwd = 1
  )
}
if (nrow(tandem_link) > 0) {
  circos.genomicLink(
    tandem_link %>% select(Chr1, Start1, End1) %>% as.data.frame(),
    tandem_link %>% select(Chr2, Start2, End2) %>% as.data.frame(),
    col = link_color, border = "#fec44f", lwd = 1
  )
}

# Legends
gene_legend <- Legend(
  at = c(0, 1, 2, 3, 4),
  labels = c("Singleton", "Dispersed", "Proximal", "Tandem", "WGD/segmental"),
  labels_gp = gpar(fontsize = 8),
  title = "Gene Type", title_gp = gpar(fontsize = 9),
  grid_height = unit(0.4, "cm"), grid_width = unit(0.4, "cm"),
  type = "points", pch = NA,
  background = c("#00ADFF", "#e66101", "#fdb863", "#b2abd2", "#5e3c99")
)
dup_legend <- Legend(
  at = c(1, 2), labels = c("Collinearity", "Tandem"),
  labels_gp = gpar(fontsize = 8),
  grid_height = unit(0.5, "cm"), grid_width = unit(0.5, "cm"),
  type = "lines", pch = NA,
  legend_gp = gpar(col = c("#66c2a4", "#fec44f"), lwd = 1),
  title = "Gene Duplicate", title_gp = gpar(fontsize = 9)
)

pushViewport(viewport(x = 0.9, y = 0.9)); grid.draw(gene_legend); upViewport()
pushViewport(viewport(x = 0.9, y = 0.15)); grid.draw(dup_legend);  upViewport()

circos.clear()
dev.off()

message(sprintf(
  "Wrote %s | family=%d  tandem_pairs=%d  col_pairs=%d",
  pdf_path, nrow(family_bed), nrow(tandem_link), nrow(col_link)
))
