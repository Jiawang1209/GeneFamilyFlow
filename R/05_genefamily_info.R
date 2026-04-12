#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Biostrings)
  library(Peptides)
  library(tidyverse)
  library(writexl)
})

option_list <- list(
  make_option("--input_fasta", type = "character",
              help = "Protein FASTA file (identify.ID.clean.fa)"),
  make_option("--input_bed", type = "character",
              help = "BED file with gene coordinates (species_10.bed)"),
  make_option("--species_map", type = "character", default = NULL,
              help = "Comma-separated prefix=name pairs, e.g. 'AT=Arabidopsis_thaliana,Sobic=Sorghum_bicolor,LOC_Os=Oryza_sativa'"),
  make_option("--outdir", type = "character",
              help = "Output directory"),
  make_option("--width", type = "double", default = 12.5),
  make_option("--height", type = "double", default = 7.5)
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input_fasta) || is.null(opt$input_bed) || is.null(opt$outdir)) {
  stop("Required: --input_fasta, --input_bed, --outdir")
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

# --- Load data ---
pep.fa <- readAAStringSet(opt$input_fasta)

bed <- read_delim(opt$input_bed, col_names = FALSE, delim = "\t",
                  show_col_types = FALSE) %>%
  set_names(c("Chr", "Start", "End", "ID", "Info", "Strand"))

# --- Compute physicochemical properties ---
df <- tibble(ID = names(pep.fa), seq = as.character(pep.fa)) %>%
  mutate(
    Length = Peptides::lengthpep(seq),
    MW     = Peptides::mw(seq),
    hydrophobicity = Peptides::hydrophobicity(seq),
    pI     = Peptides::pI(seq),
    Species = assign_species(ID, sp_map)
  ) %>%
  left_join(bed, by = "ID")

# --- Write tables ---
df_out <- df %>% select(-seq)
write_xlsx(df_out, path = file.path(opt$outdir, "Gene_Information.xlsx"))
write_csv(df_out, file.path(opt$outdir, "Gene_Information.csv"))

cat(sprintf("Gene info: %d genes written\n", nrow(df_out)))

# --- Boxplot visualization (if species info available) ---
if (!is.null(sp_map) && any(!is.na(df$Species))) {
  df_plot <- df %>%
    filter(!is.na(Species)) %>%
    select(Species, Length, MW, hydrophobicity, pI) %>%
    pivot_longer(cols = -Species, names_to = "Property", values_to = "Value") %>%
    mutate(Property = factor(Property,
                             levels = c("Length", "MW", "hydrophobicity", "pI")))

  n_sp <- n_distinct(df_plot$Species)
  palette <- scales::hue_pal()(n_sp)

  p <- ggplot(df_plot, aes(x = Species, y = Value, fill = Species)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(shape = 21, size = 2, alpha = 0.5, width = 0.2) +
    facet_wrap(~Property, scales = "free", nrow = 2) +
    labs(x = "", y = "") +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(color = "#000000", size = 10),
      axis.text.x = element_text(angle = 15, hjust = 1),
      panel.border = element_rect(linewidth = 1),
      strip.text = element_text(color = "#000000", size = 14),
      strip.background = element_rect(linewidth = 1)
    )

  ggsave(file.path(opt$outdir, "species_gene_info.pdf"),
         plot = p, width = opt$width, height = opt$height)
  cat("Boxplot saved: species_gene_info.pdf\n")
}
