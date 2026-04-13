#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(readxl)
  library(ggpubr)
})

option_list <- list(
  make_option("--expression", type = "character",
              help = "Wide-format xlsx: one ID column + one column per group, rows are replicates"),
  make_option("--sheet", type = "character", default = "Sheet5",
              help = "Sheet name in the xlsx [default: Sheet5]"),
  make_option("--id_col", type = "character", default = "ID",
              help = "Gene ID column name [default: ID]"),
  make_option("--group_cols", type = "character", default = NULL,
              help = paste("Comma-separated group column names in display order.",
                           "Default: all non-ID columns")),
  make_option("--comparisons", type = "character", default = NULL,
              help = paste("Comma-separated 'A:B,C:D' pairs for pairwise t-test",
                           "significance marks on the plot. Optional.")),
  make_option("--outdir", type = "character",
              help = "Output directory"),
  make_option("--width", type = "double", default = 12),
  make_option("--height", type = "double", default = 8),
  make_option("--ncol", type = "integer", default = 4,
              help = "Facet columns [default: 4]")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$expression) || is.null(opt$outdir)) {
  stop("Required: --expression, --outdir")
}
if (!file.exists(opt$expression)) {
  stop(sprintf("Expression file not found: %s", opt$expression))
}

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# --- Load wide-format expression table ---
df_wide <- suppressMessages(
  read_xlsx(opt$expression, sheet = opt$sheet)
)

if (!opt$id_col %in% colnames(df_wide)) {
  stop(sprintf("ID column '%s' not found in sheet '%s' (columns: %s)",
               opt$id_col, opt$sheet, paste(colnames(df_wide), collapse = ", ")))
}

group_cols <- if (is.null(opt$group_cols) || opt$group_cols == "") {
  setdiff(colnames(df_wide), opt$id_col)
} else {
  strsplit(opt$group_cols, ",")[[1]]
}

missing_cols <- setdiff(group_cols, colnames(df_wide))
if (length(missing_cols) > 0) {
  stop(sprintf("Group columns not in sheet: %s",
               paste(missing_cols, collapse = ", ")))
}

# --- Reshape wide → long (Gene, Group, Replicate, Expression) ---
df_long <- df_wide %>%
  select(all_of(c(opt$id_col, group_cols))) %>%
  rename(Gene = !!opt$id_col) %>%
  filter(!is.na(Gene)) %>%
  group_by(Gene) %>%
  mutate(Replicate = row_number()) %>%
  ungroup() %>%
  pivot_longer(
    cols = -c(Gene, Replicate),
    names_to = "Group",
    values_to = "Expression"
  ) %>%
  filter(!is.na(Expression)) %>%
  mutate(Group = factor(Group, levels = group_cols))

if (nrow(df_long) == 0) {
  stop("No non-NA expression values after reshape; check sheet and group columns")
}

# --- Per-group summary (mean / SE / n) ---
df_summary <- df_long %>%
  group_by(Gene, Group) %>%
  summarise(
    mean = mean(Expression, na.rm = TRUE),
    sd   = sd(Expression, na.rm = TRUE),
    n    = sum(!is.na(Expression)),
    se   = sd / sqrt(n),
    .groups = "drop"
  )

summary_path <- file.path(opt$outdir, "qRT_PCR_summary.csv")
write_csv(df_summary, summary_path)

# --- Parse t-test comparisons ---
comparison_list <- NULL
if (!is.null(opt$comparisons) && nzchar(opt$comparisons)) {
  pairs_raw <- strsplit(opt$comparisons, ",")[[1]]
  comparison_list <- lapply(pairs_raw, function(pair) {
    parts <- strsplit(pair, ":")[[1]]
    if (length(parts) != 2) {
      stop(sprintf("Bad comparison spec '%s'; expected 'A:B'", pair))
    }
    unknown <- setdiff(parts, group_cols)
    if (length(unknown) > 0) {
      stop(sprintf("Comparison references unknown group(s): %s",
                   paste(unknown, collapse = ", ")))
    }
    parts
  })
}

# --- Build facet bar plot ---
pal <- if (length(group_cols) <= 8L) {
  RColorBrewer::brewer.pal(max(3L, length(group_cols)), "Set2")[seq_along(group_cols)]
} else {
  scales::hue_pal()(length(group_cols))
}
names(pal) <- group_cols

p <- ggplot(df_long, aes(x = Group, y = Expression, fill = Group)) +
  stat_summary(fun = mean, geom = "col",
               color = "grey20", width = 0.72) +
  stat_summary(fun.data = mean_se, geom = "errorbar",
               width = 0.22, color = "grey25", linewidth = 0.5) +
  geom_jitter(width = 0.15, size = 1.1, alpha = 0.65,
              color = "grey20") +
  facet_wrap(~Gene, scales = "free_y", ncol = opt$ncol) +
  scale_fill_manual(values = pal) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x      = element_text(angle = 40, hjust = 1),
    legend.position  = "none",
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text       = element_text(face = "italic", size = 9),
    panel.grid.minor = element_blank()
  ) +
  labs(x = NULL, y = "Relative expression")

if (!is.null(comparison_list) && length(comparison_list) > 0L) {
  p <- p + ggpubr::stat_compare_means(
    comparisons = comparison_list,
    method      = "t.test",
    label       = "p.signif",
    hide.ns     = TRUE,
    size        = 3,
    tip.length  = 0.01
  )
}

out_pdf <- file.path(opt$outdir, "qRT_PCR.pdf")
ggsave(out_pdf, plot = p,
       width = opt$width, height = opt$height, limitsize = FALSE)

cat(sprintf(
  "qRT-PCR plotted: %d genes × %d groups → %s\n",
  dplyr::n_distinct(df_summary$Gene),
  dplyr::n_distinct(df_summary$Group),
  out_pdf
))
cat(sprintf("Summary table → %s\n", summary_path))
