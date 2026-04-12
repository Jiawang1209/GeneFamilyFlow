#!/usr/bin/env Rscript
# 08_jcvi_kaks.R — Ka/Ks distribution plot (boxplot + beeswarm)
#
# Inputs:
#   --kaks_file     Tab-separated file with columns: Seq_1, Seq_2, Ka, Ks, Ka_Ks (and more)
#   --species_map   Comma-separated prefix=species, e.g.
#                   "AT=Arabidopsis_thaliana,LOC_Os=Oryza_sativa,Sobic=Sorghum_bicolor"
#   --outdir        Output directory
#   --width, --height, --output_name   Optional plotting options

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(ggbeeswarm)
})

option_list <- list(
  make_option("--kaks_file",   type = "character", help = "Path to kaks.tab.xls"),
  make_option("--species_map", type = "character",
              help = "Comma-separated prefix=species mapping"),
  make_option("--outdir",      type = "character", default = "output/08_jcvi_kaks"),
  make_option("--output_name", type = "character", default = "08_jcvi_kaks.pdf"),
  make_option("--width",       type = "double",    default = 12),
  make_option("--height",      type = "double",    default = 7.5)
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$kaks_file) || is.null(opt$species_map)) {
  stop("--kaks_file and --species_map are required")
}
if (!file.exists(opt$kaks_file)) {
  stop(sprintf("kaks_file not found: %s", opt$kaks_file))
}
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# ---- Parse species map ---------------------------------------------------
parse_species_map <- function(s) {
  entries <- str_split(s, ",", simplify = TRUE) %>% str_trim()
  entries <- entries[nzchar(entries)]
  kv <- str_split_fixed(entries, "=", 2)
  if (any(!nzchar(kv[, 2]))) {
    stop("Invalid --species_map entry; expected 'prefix=species'")
  }
  tibble(prefix = kv[, 1], species = kv[, 2]) %>%
    arrange(desc(nchar(prefix)))  # longest-prefix match first
}

species_tbl <- parse_species_map(opt$species_map)

assign_species <- function(ids, tbl) {
  out <- rep(NA_character_, length(ids))
  for (i in seq_len(nrow(tbl))) {
    mask <- is.na(out) & str_starts(ids, fixed(tbl$prefix[i]))
    out[mask] <- tbl$species[i]
  }
  out
}

# ---- Load and tidy -------------------------------------------------------
kaks <- read_delim(opt$kaks_file, delim = "\t", show_col_types = FALSE) %>%
  filter(Seq_1 != "Seq_1") %>%
  select(Seq_1, Seq_2, Ka, Ks, Ka_Ks) %>%
  mutate(across(c(Ka, Ks, Ka_Ks), ~ suppressWarnings(as.numeric(.)))) %>%
  mutate(
    Species_1 = assign_species(Seq_1, species_tbl),
    Species_2 = assign_species(Seq_2, species_tbl)
  ) %>%
  filter(!is.na(Species_1), !is.na(Species_2)) %>%
  mutate(
    Pair = map2_chr(Species_1, Species_2, ~ paste(sort(c(.x, .y)), collapse = " vs "))
  ) %>%
  pivot_longer(c(Ka, Ks, Ka_Ks), names_to = "Type", values_to = "Value") %>%
  filter(!is.na(Value), is.finite(Value)) %>%
  mutate(Type = factor(Type, levels = c("Ka", "Ks", "Ka_Ks")))

if (nrow(kaks) == 0) {
  warning("No usable Ka/Ks rows after filtering; writing empty plot.")
}

# ---- Plot ----------------------------------------------------------------
palette <- c("#4eb3d3", "#dd3497", "#807dba", "#66c2a5", "#fc8d62", "#8da0cb")
pair_levels <- sort(unique(kaks$Pair))
fill_vals <- setNames(palette[seq_along(pair_levels)], pair_levels)

p <- ggplot(kaks, aes(x = Pair, y = Value, fill = Pair)) +
  geom_boxplot(outlier.colour = NA, alpha = 0.65, width = 0.6) +
  geom_quasirandom(aes(shape = Type), size = 1.8, alpha = 0.6,
                   varwidth = TRUE, color = "#000000") +
  geom_hline(yintercept = 1, linetype = 2, color = "grey40") +
  facet_wrap(~ Type, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = fill_vals) +
  scale_shape_manual(values = 21:23) +
  labs(x = NULL, y = "Value", fill = "Species pair") +
  theme_bw() +
  theme(
    panel.grid      = element_blank(),
    panel.border    = element_rect(linewidth = 1),
    legend.position = "bottom",
    axis.text.x     = element_text(angle = 30, hjust = 1),
    strip.text      = element_text(size = 13),
    plot.margin     = margin(5, 5, 5, 20, unit = "pt")
  )

out_pdf <- file.path(opt$outdir, opt$output_name)
ggsave(out_pdf, p, width = opt$width, height = opt$height)

# Also write the tidied table for downstream inspection.
write_csv(kaks, file.path(opt$outdir, "kaks_tidy.csv"))

message(sprintf("Wrote %s (%d rows)", out_pdf, nrow(kaks)))
