#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(ggtree)
  library(treeio)
  library(ggnewscale)
  library(readxl)
})

option_list <- list(
  make_option("--treefile", type = "character",
              help = "Newick tree file (e.g. from IQ-TREE)"),
  make_option("--group_file", type = "character", default = NULL,
              help = "Optional: tree_group_2.xlsx with columns 'label', 'Group', 'Species', 'rename'"),
  make_option("--species_map", type = "character", default = NULL,
              help = "Comma-separated prefix=name pairs, e.g. 'AT=Arabidopsis_thaliana,Sobic=Sorghum_bicolor'"),
  make_option("--layout", type = "character", default = "circular",
              help = "Tree layout: circular or rectangular [default: circular]"),
  make_option("--outdir", type = "character",
              help = "Output directory"),
  make_option("--width", type = "double", default = 15),
  make_option("--height", type = "double", default = 15)
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$treefile) || is.null(opt$outdir)) {
  stop("Required: --treefile, --outdir")
}

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# --- Parse species map ---
parse_species_map <- function(map_str) {
  if (is.null(map_str) || map_str == "") return(NULL)
  pairs <- strsplit(map_str, ",")[[1]]
  keys <- sub("=.*", "", pairs)
  vals <- sub(".*=", "", pairs)
  setNames(vals, keys)
}

assign_species <- function(ids, sp_map) {
  if (is.null(sp_map)) return(rep(NA_character_, length(ids)))
  species <- rep(NA_character_, length(ids))
  for (prefix in names(sp_map)) {
    species[str_starts(ids, fixed(prefix))] <- sp_map[[prefix]]
  }
  species
}

sp_map <- parse_species_map(opt$species_map)

# --- Load tree ---
tree <- read.newick(opt$treefile, node.label = "support")
tip_labels <- tree$tip.label

# --- Build tip metadata ---
if (!is.null(opt$group_file) && file.exists(opt$group_file)) {
  tip_meta <- read_xlsx(opt$group_file) %>%
    filter(label %in% tip_labels)
  has_groups <- "Group" %in% colnames(tip_meta) && any(!is.na(tip_meta$Group))
  has_rename <- "rename" %in% colnames(tip_meta) && any(!is.na(tip_meta$rename))
} else {
  tip_meta <- tibble(label = tip_labels)
  if (!is.null(sp_map)) {
    tip_meta <- tip_meta %>%
      mutate(Species = assign_species(label, sp_map))
  }
  has_groups <- FALSE
  has_rename <- FALSE
}

# --- Color palettes ---
subfamily_colors <- c(
  "NPF1" = "#78c679", "NPF2" = "#238443", "NPF3" = "#c2a5cf",
  "NPF4" = "#88419d", "NPF5" = "#fec44f", "NPF6" = "#cc4c02",
  "NPF7" = "#7bccc4", "NPF8" = "#0868ac"
)

species_colors <- c(
  "Arabidopsis_thaliana"  = "#de77ae",
  "Sorghum_bicolor"       = "#d6604d",
  "Oryza_sativa"          = "#66bd63",
  "Oryza_sativa_Japonica" = "#66bd63"
)

# --- Build tree plot ---
p <- ggtree(tree,
            branch.length = "none",
            color = "#969696",
            layout = opt$layout) %<+% tip_meta +
  geom_tree(size = 0.3, color = "#969696")

# Bootstrap support
p <- p +
  geom_nodepoint(
    aes(fill = cut(support, c(0, 45, 75, 100))),
    shape = 21, size = 1.5, stroke = 0.2
  ) +
  scale_fill_manual(
    values = c("black", "grey", "white"),
    name = "Bootstrap (%)",
    breaks = c("(75,100]", "(45,75]", "(0,45]"),
    labels = expression(BP >= 75, 45 <= BP * " < 75", BP < 45),
    guide = guide_legend(override.aes = list(size = 3)),
    na.value = NA
  )

# Tip labels
tip_label_col <- if (has_rename) "rename" else "label"
if (has_groups) {
  p <- p +
    geom_tiplab(aes(x = x + 0.5, label = .data[[tip_label_col]],
                    color = Group), size = 2) +
    scale_color_manual(values = subfamily_colors, name = "Subfamily",
                       na.value = "grey50")
} else {
  p <- p +
    geom_tiplab(aes(x = x + 0.5, label = .data[[tip_label_col]]),
                size = 2, color = "grey30")
}

# Species tip points
if ("Species" %in% colnames(tip_meta) && any(!is.na(tip_meta$Species))) {
  all_sp <- unique(na.omit(tip_meta$Species))
  sp_pal <- species_colors[all_sp]
  missing <- setdiff(all_sp, names(species_colors))
  if (length(missing) > 0) {
    extra <- scales::hue_pal()(length(missing))
    sp_pal <- c(sp_pal, setNames(extra, missing))
  }
  sp_pal <- sp_pal[!is.na(names(sp_pal))]

  p <- p +
    new_scale_fill() +
    geom_tippoint(aes(fill = Species), shape = 21, size = 2) +
    scale_fill_manual(values = sp_pal, name = "Species")
}

xlim_val <- if (opt$layout == "circular") 30 else max(50, length(tip_labels) * 0.3)
p <- p + xlim(NA, xlim_val) +
  theme(legend.background = element_blank())

# --- Save ---
ggsave(file.path(opt$outdir, "phylogenetic_tree.pdf"),
       plot = p, width = opt$width, height = opt$height,
       limitsize = FALSE)

cat(sprintf("Tree plotted: %d tips, layout=%s\n", length(tip_labels), opt$layout))
