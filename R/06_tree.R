#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(ggtree)
  library(treeio)
  library(ggnewscale)
  library(readxl)
  library(RColorBrewer)
  library(ape)
})

option_list <- list(
  make_option("--treefile", type = "character",
              help = "Newick tree file (e.g. from IQ-TREE)"),
  make_option("--group_file", type = "character", default = NULL,
              help = "Optional: tree_group_2.xlsx with columns 'label', 'Group', 'Species', 'rename'"),
  make_option("--species_map", type = "character", default = NULL,
              help = "Comma-separated prefix=name pairs, e.g. 'AT=Arabidopsis_thaliana,Sobic=Sorghum_bicolor'"),
  make_option("--layout", type = "character", default = "both",
              help = "Tree layout: circular, rectangular, or both [default: both]"),
  make_option("--outdir", type = "character",
              help = "Output directory"),
  make_option("--width", type = "double", default = 15),
  make_option("--height", type = "double", default = 15),
  make_option("--n_subfamilies", type = "integer", default = 8,
              help = paste(
                "Auto-split tips into N monophyletic subfamilies by recursive",
                "top-down clade splitting when no --group_file is given.",
                "Set 0 to disable. [default: 8]"
              )),
  make_option("--subfamily_prefix", type = "character", default = "Sub",
              help = "Label prefix for auto subfamilies [default: Sub]")
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
phy <- if (methods::is(tree, "treedata")) tree@phylo else tree
tip_labels <- phy$tip.label
n_tips <- length(tip_labels)

# --- Top-down monophyletic clade splitter ---
# Starts from the root as a single clade, then repeatedly splits the clade
# with the most tips into its direct children until we have k clades. This
# guarantees every returned clade is monophyletic, unlike cophenetic-distance
# hierarchical clustering which can produce polyphyletic groups that don't
# hilight cleanly on the tree.
split_tree_into_k_clades <- function(phy, k) {
  n_tips <- length(phy$tip.label)
  k <- min(k, n_tips)
  root_node <- n_tips + 1L
  clades <- list(list(node = root_node, tips = phy$tip.label))

  while (length(clades) < k) {
    sizes <- vapply(clades, function(cl) length(cl$tips), integer(1))
    splittable <- sizes > 1L
    if (!any(splittable)) break
    sizes[!splittable] <- 0L
    biggest_idx <- which.max(sizes)
    parent <- clades[[biggest_idx]]

    children <- phy$edge[phy$edge[, 1] == parent$node, 2]
    child_clades <- lapply(children, function(cn) {
      if (cn <= n_tips) {
        list(node = as.integer(cn), tips = phy$tip.label[cn])
      } else {
        sub <- ape::extract.clade(phy, cn)
        list(node = as.integer(cn), tips = sub$tip.label)
      }
    })

    clades[[biggest_idx]] <- NULL
    clades <- c(clades, child_clades)
  }
  clades
}

# --- Build tip metadata + clade nodes for hilight ---
clade_nodes <- NULL

if (!is.null(opt$group_file) && file.exists(opt$group_file)) {
  tip_meta <- read_xlsx(opt$group_file) %>%
    filter(label %in% tip_labels)
  has_groups <- "Group" %in% colnames(tip_meta) && any(!is.na(tip_meta$Group))
  has_rename <- "rename" %in% colnames(tip_meta) && any(!is.na(tip_meta$rename))

  # Resolve MRCA node for each manual group (for geom_hilight).
  if (has_groups) {
    groups_vec <- sort(unique(na.omit(tip_meta$Group)))
    clade_nodes <- vapply(groups_vec, function(g) {
      members <- tip_meta$label[!is.na(tip_meta$Group) & tip_meta$Group == g]
      members <- members[members %in% tip_labels]
      if (length(members) < 2) {
        as.integer(match(members, tip_labels))
      } else {
        as.integer(ape::getMRCA(phy, members))
      }
    }, integer(1))
    names(clade_nodes) <- groups_vec
  }
} else {
  tip_meta <- tibble(label = tip_labels)
  if (!is.null(sp_map)) {
    tip_meta <- tip_meta %>%
      mutate(Species = assign_species(label, sp_map))
  }
  has_rename <- FALSE

  n_clusters <- min(opt$n_subfamilies, n_tips)
  if (n_clusters >= 2) {
    clade_list <- split_tree_into_k_clades(phy, n_clusters)
    # Order clades by first tip position so numbering matches tree order.
    first_tip_idx <- vapply(clade_list, function(cl) {
      min(match(cl$tips, tip_labels))
    }, integer(1))
    clade_list <- clade_list[order(first_tip_idx)]

    group_labels <- rep(NA_character_, n_tips)
    nodes_v <- integer(length(clade_list))
    for (i in seq_along(clade_list)) {
      cl <- clade_list[[i]]
      group_labels[tip_labels %in% cl$tips] <-
        sprintf("%s%d", opt$subfamily_prefix, i)
      nodes_v[i] <- cl$node
    }
    names(nodes_v) <- sprintf(
      "%s%d", opt$subfamily_prefix, seq_along(clade_list)
    )
    clade_nodes <- nodes_v

    tip_meta <- tip_meta %>% mutate(Group = group_labels)
    has_groups <- TRUE
    cat(sprintf(
      "Auto-assigned %d monophyletic subfamilies (%s1..%s%d)\n",
      length(clade_list), opt$subfamily_prefix, opt$subfamily_prefix,
      length(clade_list)
    ))
  } else {
    has_groups <- FALSE
  }
}

# --- Color palette (shared across branches, tip labels, hilight, cladelab) ---
build_subfamily_colors <- function(groups) {
  levels_g <- sort(unique(na.omit(groups)))
  n <- length(levels_g)
  if (n == 0) return(character(0))
  pal <- if (n <= 12L) {
    RColorBrewer::brewer.pal(max(n, 3L), "Set3")[seq_len(n)]
  } else if (n <= 20L) {
    c(
      RColorBrewer::brewer.pal(12L, "Set3"),
      RColorBrewer::brewer.pal(8L, "Set2")
    )[seq_len(n)]
  } else {
    scales::hue_pal()(n)
  }
  setNames(pal, levels_g)
}

subfamily_colors <- if (has_groups) {
  build_subfamily_colors(tip_meta$Group)
} else {
  character(0)
}

species_colors <- c(
  "Arabidopsis_thaliana"  = "#de77ae",
  "Sorghum_bicolor"       = "#d6604d",
  "Oryza_sativa"          = "#66bd63",
  "Oryza_sativa_Japonica" = "#66bd63"
)

# --- Attach group info to tree via groupOTU so branches color by subfamily ---
tree_grouped <- if (has_groups) {
  split_list <- split(tip_meta$label, tip_meta$Group)
  groupOTU(tree, split_list, group_name = "Subfamily")
} else {
  tree
}

# --- Build a single-layout plot ---
build_plot <- function(layout) {
  if (has_groups) {
    p <- ggtree(
      tree_grouped,
      branch.length = "none",
      layout        = layout,
      size          = 0.45,
      aes(color = Subfamily)
    ) %<+% tip_meta +
      scale_color_manual(
        values    = c(subfamily_colors, "0" = "grey60"),
        name      = "Subfamily",
        breaks    = names(subfamily_colors),
        na.value  = "grey60"
      )
  } else {
    p <- ggtree(
      tree_grouped,
      branch.length = "none",
      layout        = layout,
      color         = "grey40",
      size          = 0.45
    ) %<+% tip_meta
  }

  # Subfamily background highlights (drawn under everything else).
  if (has_groups && !is.null(clade_nodes) && length(clade_nodes) > 0L) {
    for (i in seq_along(clade_nodes)) {
      sub_name <- names(clade_nodes)[i]
      p <- p + geom_hilight(
        node   = clade_nodes[[i]],
        fill   = unname(subfamily_colors[sub_name]),
        alpha  = 0.18,
        extend = if (layout == "circular") 2 else 3
      )
    }
  }

  # Bootstrap support dots (new fill scale so it doesn't fight Species).
  p <- p + new_scale_fill() +
    geom_nodepoint(
      aes(fill = cut(support, c(0, 45, 75, 100))),
      shape = 21, size = 1.6, stroke = 0.25, color = "grey30"
    ) +
    scale_fill_manual(
      values = c("black", "grey", "white"),
      name   = "Bootstrap (%)",
      breaks = c("(75,100]", "(45,75]", "(0,45]"),
      labels = expression(BP >= 75, 45 <= BP * " < 75", BP < 45),
      guide  = guide_legend(override.aes = list(size = 3)),
      na.value = NA
    )

  # Tip labels (share branch color scale so label color matches branch color).
  tip_label_col <- if (has_rename) "rename" else "label"
  if (has_groups) {
    p <- p + geom_tiplab(
      aes(x = x + 0.5,
          label = .data[[tip_label_col]],
          color = Subfamily),
      size = 2,
      show.legend = FALSE
    )
  } else {
    p <- p + geom_tiplab(
      aes(x = x + 0.5, label = .data[[tip_label_col]]),
      size = 2, color = "grey30"
    )
  }

  # Species tip points (fresh fill scale, independent of bootstrap fill).
  if ("Species" %in% colnames(tip_meta) && any(!is.na(tip_meta$Species))) {
    all_sp <- unique(na.omit(tip_meta$Species))
    sp_pal <- species_colors[all_sp]
    missing <- setdiff(all_sp, names(species_colors))
    if (length(missing) > 0L) {
      extra <- scales::hue_pal()(length(missing))
      sp_pal <- c(sp_pal, setNames(extra, missing))
    }
    sp_pal <- sp_pal[!is.na(names(sp_pal))]

    p <- p + new_scale_fill() +
      geom_tippoint(aes(fill = Species), shape = 21,
                    size = 2, color = "grey30") +
      scale_fill_manual(values = sp_pal, name = "Species")
  }

  # Subfamily clade labels (outside tip labels, at layout-appropriate offset).
  if (has_groups && !is.null(clade_nodes) && length(clade_nodes) > 0L) {
    for (i in seq_along(clade_nodes)) {
      sub_name <- names(clade_nodes)[i]
      col <- unname(subfamily_colors[sub_name])
      p <- p + geom_cladelab(
        node        = clade_nodes[[i]],
        label       = sub_name,
        align       = TRUE,
        offset      = if (layout == "circular") 4 else 6,
        offset.text = 0.6,
        barsize     = 1.4,
        fontsize    = 3.5,
        barcolor    = col,
        textcolor   = col,
        angle       = if (layout == "circular") "auto" else 0
      )
    }
  }

  # Adaptive xlim to leave room for tip labels + cladelabs.
  xlim_pad <- if (layout == "circular") {
    max(25, n_tips * 0.18)
  } else {
    max(50, n_tips * 0.35)
  }
  p <- p + xlim(NA, xlim_pad) +
    theme(legend.background = element_blank())

  p
}

# --- Render each requested layout ---
layouts <- if (opt$layout == "both") {
  c("circular", "rectangular")
} else {
  opt$layout
}

for (lay in layouts) {
  p <- build_plot(lay)
  out_path <- file.path(opt$outdir,
                        sprintf("phylogenetic_tree_%s.pdf", lay))
  ggsave(out_path, plot = p,
         width = opt$width, height = opt$height,
         limitsize = FALSE)
  cat(sprintf("Tree plotted: %d tips, layout=%s → %s\n",
              n_tips, lay, out_path))
}
