#!/usr/bin/env Rscript
#
# Composite figure: Tree | Domain | Motif | Gene Structure
# Aligned by gene order from the phylogenetic tree.

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(ggtree)
  library(treeio)
  library(ggnewscale)
  library(gggenes)
  library(aplot)
  library(readxl)
})

option_list <- list(
  make_option("--treefile", type = "character",
              help = "Newick tree file (target species only)"),
  make_option("--group_file", type = "character", default = NULL,
              help = "Optional: tree_group_2.xlsx for subfamily coloring"),
  make_option("--pfam_scan", type = "character", default = NULL,
              help = "Pfam_scan.out file for domain annotation"),
  make_option("--motif_file", type = "character",
              help = "meme_location.txt (Gene\\tStart\\tEnd\\tMotif)"),
  make_option("--gff3", type = "character", default = NULL,
              help = "GFF3 file for gene structure (CDS features)"),
  make_option("--id_regex_strip", type = "character",
              default = "\\.v\\d+(\\.\\d+)?$|\\.\\d+\\.p$|\\.p$",
              help = "Regex to strip from GFF3 transcript IDs to match tree tips"),
  make_option("--outdir", type = "character",
              help = "Output directory"),
  make_option("--width", type = "double", default = 22.5),
  make_option("--height", type = "double", default = 20)
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$treefile) || is.null(opt$motif_file) || is.null(opt$outdir)) {
  stop("Required: --treefile, --motif_file, --outdir")
}

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# --- Color palettes ---
panel_colors <- c(
  "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99",
  "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a",
  "#ffff99", "#b15928", "#ffffb3", "#bebada", "#fb8072"
)

subfamily_colors <- c(
  "NPF1" = "#78c679", "NPF2" = "#238443", "NPF3" = "#c2a5cf",
  "NPF4" = "#88419d", "NPF5" = "#fec44f", "NPF6" = "#cc4c02",
  "NPF7" = "#7bccc4", "NPF8" = "#0868ac"
)

# =========================================================================
# 1. Tree
# =========================================================================
tree <- read.newick(opt$treefile, node.label = "support")

if (!is.null(opt$group_file) && file.exists(opt$group_file)) {
  tip_meta <- read_xlsx(opt$group_file) %>%
    filter(label %in% tree$tip.label)
  has_groups <- "Group" %in% colnames(tip_meta) && any(!is.na(tip_meta$Group))
  has_rename <- "rename" %in% colnames(tip_meta) && any(!is.na(tip_meta$rename))
} else {
  tip_meta <- tibble(label = tree$tip.label)
  has_groups <- FALSE
  has_rename <- FALSE
}

tip_label_col <- if (has_rename) "rename" else "label"

p_tree <- ggtree(tree, branch.length = "none", color = "#969696") %<+% tip_meta +
  geom_tree(size = 0.3, color = "#969696") +
  geom_nodepoint(
    aes(fill = cut(support, c(0, 45, 75, 100))),
    shape = 21, size = 2, stroke = 0.2
  ) +
  scale_fill_manual(
    values = c("black", "grey", "white"),
    name = "Bootstrap (%)",
    breaks = c("(75,100]", "(45,75]", "(0,45]"),
    labels = expression(BP >= 75, 45 <= BP * " < 75", BP < 45),
    guide = guide_legend(override.aes = list(size = 3)),
    na.value = NA
  )

if (has_groups) {
  p_tree <- p_tree +
    geom_tiplab(aes(x = x + 0.25, label = .data[[tip_label_col]],
                    color = Group), size = 2.5) +
    scale_color_manual(values = subfamily_colors, name = "Subfamily",
                       na.value = "grey50")
} else {
  p_tree <- p_tree +
    geom_tiplab(aes(x = x + 0.25, label = .data[[tip_label_col]]),
                size = 2.5, color = "grey30")
}

p_tree <- p_tree +
  xlim(NA, 25) +
  theme(legend.background = element_blank())

# Gene order from tree (top to bottom)
tree_order <- get_taxa_name(p_tree)

# =========================================================================
# 2. Domain panel (optional)
# =========================================================================
p_domain <- NULL
if (!is.null(opt$pfam_scan) && file.exists(opt$pfam_scan)) {
  domain_df <- read.table(opt$pfam_scan, comment.char = "#", header = FALSE,
                          fill = TRUE) %>%
    filter(V1 %in% tree_order) %>%
    mutate(V1 = factor(V1, levels = rev(tree_order), ordered = TRUE))

  p_domain <- ggplot(domain_df, aes(xmin = V4, xmax = V5, y = V1, fill = V7)) +
    geom_gene_arrow(arrowhead_height = unit(3, "mm"),
                    arrowhead_width = unit(0.1, "mm")) +
    labs(x = "", y = "") +
    ggtitle("Domain") +
    scale_fill_manual(values = panel_colors, name = "Domain") +
    theme_genes() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(hjust = 0.5))
}

# =========================================================================
# 3. Motif panel
# =========================================================================
motif_df <- read_delim(opt$motif_file, delim = "\t",
                       col_names = c("Gene", "Start", "End", "Motif"),
                       show_col_types = FALSE) %>%
  filter(Gene %in% tree_order) %>%
  mutate(
    Gene = factor(Gene, levels = rev(tree_order), ordered = TRUE),
    Motif = factor(Motif,
                   levels = paste0("MEME-", 1:15),
                   ordered = TRUE)
  )

p_motif <- ggplot(motif_df, aes(xmin = Start, xmax = End, y = Gene, fill = Motif)) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"),
                  arrowhead_width = unit(0.1, "mm")) +
  labs(x = "", y = "") +
  ggtitle("Motif") +
  scale_fill_manual(values = panel_colors, name = "Motif") +
  theme_genes() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 0.5))

# =========================================================================
# 4. Gene structure panel (optional)
# =========================================================================
p_structure <- NULL
if (!is.null(opt$gff3) && file.exists(opt$gff3)) {
  gff_raw <- read_delim(opt$gff3, delim = "\t", col_names = FALSE,
                        comment = "#", show_col_types = FALSE) %>%
    filter(X3 == "CDS")

  # Extract transcript ID from attributes, map to tree tip IDs
  gff_cds <- gff_raw %>%
    mutate(
      attr_id = str_extract(X9, "(?<=ID=)[^;]+"),
      ID = str_replace(attr_id, opt$id_regex_strip, "")
    ) %>%
    filter(ID %in% tree_order)

  if (nrow(gff_cds) > 0) {
    # Normalize coordinates per gene (relative to first CDS start)
    gff_norm <- gff_cds %>%
      group_by(ID) %>%
      mutate(
        origin = min(X4),
        start = if_else(X7 == "+", X4 - origin, origin + (max(X5) - min(X4)) - (X5 - origin)),
        end   = if_else(X7 == "+", X5 - origin, origin + (max(X5) - min(X4)) - (X4 - origin))
      ) %>%
      mutate(start = abs(start), end = abs(end)) %>%
      ungroup() %>%
      mutate(
        ID = factor(ID, levels = rev(tree_order), ordered = TRUE),
        Type = "CDS"
      )

    p_structure <- ggplot(gff_norm, aes(xmin = start, xmax = end, y = ID, fill = Type)) +
      geom_gene_arrow(arrowhead_height = unit(3, "mm"),
                      arrowhead_width = unit(0.1, "mm")) +
      labs(x = "", y = "") +
      ggtitle("Gene Structure") +
      scale_fill_manual(values = c("CDS" = "#7fcdbb"), name = "Structure") +
      theme_genes() +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            plot.title = element_text(hjust = 0.5))
  }
}

# =========================================================================
# 5. Combine panels with aplot
# =========================================================================
panels <- list(p_motif)
widths <- c(0.75)

if (!is.null(p_domain)) {
  panels <- c(list(p_domain), panels)
  widths <- c(0.75, widths)
}

if (!is.null(p_structure)) {
  panels <- c(panels, list(p_structure))
  widths <- c(widths, 0.75)
}

p_combined <- p_tree
for (i in seq_along(panels)) {
  if (i == 1 && !is.null(p_domain)) {
    p_combined <- panels[[i]] %>% insert_left(p_combined, width = 1.5)
  } else if (i == 1) {
    p_combined <- panels[[i]] %>% insert_left(p_combined, width = 1.5)
  } else {
    p_combined <- p_combined %>% insert_right(panels[[i]], width = widths[i])
  }
}

# Fallback: if aplot chaining didn't work as expected, use simple approach
if (is.null(p_combined)) p_combined <- p_motif

ggsave(file.path(opt$outdir, "Tree_Domain_Motif_GeneStructure.pdf"),
       plot = p_combined,
       width = opt$width, height = opt$height,
       limitsize = FALSE)

cat(sprintf("Composite figure saved: %d genes, %d panels\n",
            length(tree_order),
            1 + !is.null(p_domain) + 1 + !is.null(p_structure)))
