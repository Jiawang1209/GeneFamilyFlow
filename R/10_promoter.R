#!/usr/bin/env Rscript
# 10_promoter.R — Promoter cis-element analysis (PlantCARE)
#
# Inputs:
#   --plantcare_dir   Directory with PlantCARE .tab files
#   --element_desc    Excel file mapping element -> description (category)
#   --gene_ids        Text file of family gene IDs (one per line)
#   --outdir          Output directory
#   --id_regex_strip  Regex to strip suffixes when matching gene IDs (default handles .v3.1, .1.p, .1)
#   --min_count       Min total occurrences to keep an element (default 5)
#   --top_per_cat     Max elements per category to display (default 5)
#   --width, --height

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(readxl)
  library(legendry)
  library(data.table)
})

option_list <- list(
  make_option("--plantcare_dir",  type = "character"),
  make_option("--element_desc",   type = "character"),
  make_option("--gene_ids",       type = "character"),
  make_option("--outdir",         type = "character", default = "output/10_promoter"),
  make_option("--id_regex_strip", type = "character",
              default = "\\.\\d+\\.p$|\\.v\\d+\\.\\d+$|\\.\\d+$"),
  make_option("--min_count",      type = "integer",   default = 5),
  make_option("--top_per_cat",    type = "integer",   default = 5),
  make_option("--output_name",    type = "character", default = "promoter_elements.pdf"),
  make_option("--width",          type = "double",    default = 22),
  make_option("--height",         type = "double",    default = 18)
)
opt <- parse_args(OptionParser(option_list = option_list))

for (k in c("plantcare_dir", "element_desc", "gene_ids")) {
  if (is.null(opt[[k]])) stop(sprintf("--%s is required", k))
}
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

strip_id <- function(x) str_remove(x, opt$id_regex_strip)

# ---- Load element description mapping ------------------------------------
desc_tbl <- read_xlsx(opt$element_desc) %>%
  select(element, description)

# ---- Load family gene IDs ------------------------------------------------
family_ids <- read_lines(opt$gene_ids) %>%
  str_trim() %>% discard(~ . == "") %>%
  strip_id() %>% unique()

# ---- Load and parse PlantCARE .tab files ---------------------------------
tab_files <- list.files(opt$plantcare_dir, pattern = "\\.tab$",
                        full.names = TRUE, recursive = TRUE)
if (length(tab_files) == 0) stop("No .tab files found in --plantcare_dir")

raw <- fread(file = tab_files, header = FALSE, sep = "\t", quote = "",
             fill = TRUE, select = c(1, 2, 8))
colnames(raw) <- c("raw_id", "element", "function_desc")

# Extract gene ID from PlantCARE column 1
# Format: "Chr01_xxx:+_usf:2000Sobic.001G133900.v3.1"
# Gene ID is after the last "usf:NNNN" prefix
promoter <- raw %>%
  as_tibble() %>%
  filter(!str_detect(element, "Unname"), !is.na(element)) %>%
  filter(!str_detect(function_desc, "short_function") | is.na(function_desc)) %>%
  mutate(gene_id = str_extract(raw_id, "[A-Za-z][A-Za-z0-9._]+$")) %>%
  filter(!is.na(gene_id)) %>%
  mutate(gene_id_clean = strip_id(gene_id)) %>%
  filter(gene_id_clean %in% family_ids) %>%
  inner_join(desc_tbl, by = "element")

if (nrow(promoter) == 0) {
  warning("No promoter records matched family genes.")
  quit(status = 0)
}

message(sprintf("Matched %d promoter records for %d genes",
                nrow(promoter), n_distinct(promoter$gene_id_clean)))

# ---- Filter elements by frequency ----------------------------------------
element_counts <- promoter %>%
  count(element, name = "total") %>%
  filter(total >= opt$min_count)

retain_elements <- element_counts %>%
  inner_join(desc_tbl, by = "element") %>%
  group_by(description) %>%
  slice_max(total, n = opt$top_per_cat, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(description, desc(total)) %>%
  pull(element)

# ---- Build count matrix --------------------------------------------------
count_df <- promoter %>%
  filter(element %in% retain_elements) %>%
  count(gene_id_clean, element, name = "count") %>%
  complete(gene_id_clean = family_ids[family_ids %in% unique(.$gene_id_clean)],
           element = retain_elements, fill = list(count = 0L)) %>%
  inner_join(desc_tbl, by = "element") %>%
  mutate(element = factor(element, levels = retain_elements))

# Write tidy CSV for downstream use
write_csv(count_df, file.path(opt$outdir, "promoter_tidy.csv"))

# ---- Dot heatmap plot ----------------------------------------------------
fill_breaks <- c(0, 5, 10, 20, 40, 60, 100, 200, Inf)
fill_colors <- c("#fde0ef", "#f1b6da", "#d9f0a3", "#addd8e",
                 "#78c679", "#41ab5d", "#238443", "#004529")

p_heatmap <- count_df %>%
  mutate(fill_bin = cut(count, breaks = fill_breaks,
                        right = FALSE, include.lowest = TRUE)) %>%
  ggplot(aes(x = interaction(element, description), y = gene_id_clean)) +
  geom_tile(fill = "white", color = "grey80", linewidth = 0.3,
            height = 1, width = 0.5) +
  geom_point(aes(fill = fill_bin), shape = 22, size = 5) +
  geom_text(aes(label = ifelse(count > 0, count, "")), size = 2.5) +
  scale_fill_manual(values = fill_colors, na.value = "white", name = "Count") +
  labs(x = NULL, y = NULL) +
  scale_x_discrete(
    expand = expansion(mult = c(0.01, 0.01)), position = "top",
    guide = guide_axis_nested(
      type = "bracket",
      levels_text = list(NULL, element_text(color = "#d73027", size = 10))
    )
  ) +
  theme_bw() +
  theme(
    axis.text.x       = element_text(color = "black", angle = 90, size = 9, hjust = 0),
    axis.text.y       = element_text(color = "black", size = 10),
    axis.ticks.y      = element_blank(),
    panel.border      = element_blank(),
    legend.position   = "right"
  )

ggsave(file.path(opt$outdir, opt$output_name), p_heatmap,
       width = opt$width, height = opt$height)

# ---- Category stacked bar chart ------------------------------------------
cat_colors <- c("#de77ae", "#fdb863", "#35978f", "#d6604d",
                "#7570b3", "#e7298a", "#66a61e", "#e6ab02")

p_stat <- count_df %>%
  group_by(gene_id_clean, description) %>%
  summarise(total = sum(count), .groups = "drop") %>%
  ggplot(aes(x = total, y = gene_id_clean, fill = description)) +
  geom_bar(stat = "identity", color = "black", alpha = 0.75,
           position = "fill", width = 0.6) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0)),
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = c("0%", "25%", "50%", "75%", "100%")
  ) +
  scale_fill_manual(values = cat_colors) +
  labs(x = NULL, y = NULL, fill = "Category") +
  theme_bw() +
  theme(
    axis.text.x  = element_text(color = "black", size = 12),
    axis.text.y  = element_text(color = "black", size = 10),
    panel.border = element_rect(linewidth = 1),
    axis.ticks.y = element_blank()
  )

ggsave(file.path(opt$outdir, "promoter_category_stat.pdf"), p_stat,
       width = 8, height = opt$height * 0.8)

message(sprintf("Wrote %s and promoter_category_stat.pdf (%d elements retained)",
                opt$output_name, length(retain_elements)))
