#!/usr/bin/env Rscript
# GO/KEGG over-representation enrichment using eggNOG-mapper annotations.
#
# Inputs are TERM2GENE long-format tables produced by scripts/parse_eggnog.py.
# GO ontology (BP/MF/CC) is resolved via GO.db so a single GO input is split
# into three faceted dotplots. KEGG pathway names are resolved from an offline
# mapping table; unknown codes fall back to the raw ko id.

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(clusterProfiler)
})

option_list <- list(
  make_option("--term2gene_go",   type = "character",
              help = "TERM2GENE TSV for GO (columns: term, gene)"),
  make_option("--term2gene_kegg", type = "character",
              help = "TERM2GENE TSV for KEGG (columns: term, gene)"),
  make_option("--gene_list",      type = "character",
              help = "One gene per line: the target gene set (family members in species)"),
  make_option("--universe",       type = "character",
              help = "One gene per line: background universe"),
  make_option("--kegg_names",     type = "character", default = NULL,
              help = "Optional offline KEGG pathway name TSV (columns: term, name)"),
  make_option("--species",        type = "character", default = "species",
              help = "Species label used in output filename [default %default]"),
  make_option("--outdir",         type = "character", default = ".",
              help = "Output directory [default %default]"),
  make_option("--pvalue",         type = "double",    default = 0.05,
              help = "enricher p-value cutoff [default %default]"),
  make_option("--qvalue",         type = "double",    default = 0.2,
              help = "enricher q-value cutoff [default %default]"),
  make_option("--min_size",       type = "integer",   default = 5,
              help = "Minimum gene set size [default %default]"),
  make_option("--top_n",          type = "integer",   default = 20,
              help = "Top N terms to display per category [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

required <- c("term2gene_go", "term2gene_kegg", "gene_list", "universe")
missing <- required[vapply(required, function(k) is.null(opt[[k]]), logical(1))]
if (length(missing) > 0) {
  stop("Missing required arguments: ", paste(missing, collapse = ", "))
}
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# ---- Load inputs ----------------------------------------------------------
read_term2gene <- function(path) {
  if (!file.exists(path)) return(NULL)
  df <- suppressMessages(readr::read_tsv(path, show_col_types = FALSE))
  if (nrow(df) == 0) return(NULL)
  df
}
read_gene_list <- function(path) {
  if (!file.exists(path)) return(character(0))
  ln <- readLines(path, warn = FALSE)
  unique(trimws(ln[nzchar(trimws(ln))]))
}

go_t2g   <- read_term2gene(opt$term2gene_go)
kegg_t2g <- read_term2gene(opt$term2gene_kegg)
gene_set <- read_gene_list(opt$gene_list)
universe <- read_gene_list(opt$universe)
if (length(universe) > 0) {
  gene_set <- intersect(gene_set, universe)
}

cat(sprintf(
  "[13_go_kegg] species=%s gene_set=%d universe=%d GO_pairs=%s KEGG_pairs=%s\n",
  opt$species, length(gene_set), length(universe),
  if (is.null(go_t2g)) 0 else nrow(go_t2g),
  if (is.null(kegg_t2g)) 0 else nrow(kegg_t2g)
))

# ---- Placeholder plot when nothing to report ------------------------------
empty_plot <- function(msg) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = msg, size = 5) +
    theme_void() +
    xlim(0, 1) + ylim(0, 1)
}

pdf_path <- file.path(opt$outdir, sprintf("%s_enrichment.pdf", opt$species))

if (length(gene_set) == 0) {
  pdf(pdf_path, width = 8, height = 6)
  print(empty_plot(sprintf("No gene-set members found for %s", opt$species)))
  dev.off()
  cat("[13_go_kegg] Empty gene set; wrote placeholder to", pdf_path, "\n")
  quit(status = 0)
}

# ---- GO: split by ontology via GO.db when available -----------------------
split_go_by_ontology <- function(t2g) {
  if (is.null(t2g) || nrow(t2g) == 0) return(list(BP = NULL, MF = NULL, CC = NULL))
  if (!requireNamespace("GO.db", quietly = TRUE)) {
    return(list(BP = t2g, MF = NULL, CC = NULL))  # fallback: treat as one group
  }
  ont_map <- suppressWarnings(
    AnnotationDbi::select(GO.db::GO.db, keys = unique(t2g$term),
                          columns = c("ONTOLOGY"), keytype = "GOID")
  )
  colnames(ont_map) <- c("term", "ontology")
  joined <- dplyr::left_join(t2g, ont_map, by = "term")
  list(
    BP = dplyr::filter(joined, ontology == "BP") %>% dplyr::select(term, gene),
    MF = dplyr::filter(joined, ontology == "MF") %>% dplyr::select(term, gene),
    CC = dplyr::filter(joined, ontology == "CC") %>% dplyr::select(term, gene)
  )
}

# ---- TERM2NAME resolution --------------------------------------------------
go_term2name <- function(terms) {
  if (length(terms) == 0) return(NULL)
  if (!requireNamespace("GO.db", quietly = TRUE)) return(NULL)
  m <- suppressWarnings(
    AnnotationDbi::select(GO.db::GO.db, keys = unique(terms),
                          columns = c("TERM"), keytype = "GOID")
  )
  colnames(m) <- c("term", "name")
  m
}
kegg_term2name <- function(terms) {
  if (length(terms) == 0 || is.null(opt$kegg_names) || !file.exists(opt$kegg_names)) {
    return(NULL)
  }
  tab <- suppressMessages(readr::read_tsv(opt$kegg_names, show_col_types = FALSE))
  colnames(tab)[1:2] <- c("term", "name")
  dplyr::filter(tab, term %in% terms)
}

# ---- Core enrichment wrapper -----------------------------------------------
run_enrichment <- function(t2g, t2n, label) {
  if (is.null(t2g) || nrow(t2g) == 0) {
    return(list(result = NULL, label = label))
  }
  res <- tryCatch(
    clusterProfiler::enricher(
      gene          = gene_set,
      universe      = if (length(universe) > 0) universe else NULL,
      TERM2GENE     = t2g,
      TERM2NAME     = t2n,
      pvalueCutoff  = opt$pvalue,
      qvalueCutoff  = opt$qvalue,
      minGSSize     = opt$min_size,
      pAdjustMethod = "BH"
    ),
    error = function(e) {
      message(sprintf("[13_go_kegg] %s enrichment failed: %s", label, conditionMessage(e)))
      NULL
    }
  )
  list(result = res, label = label)
}

# ---- Build the dotplot for a single enrichment result ---------------------
make_dotplot <- function(res, label, top_n) {
  if (is.null(res) || nrow(as.data.frame(res)) == 0) {
    return(empty_plot(sprintf("%s: no significant terms", label)))
  }
  df <- as.data.frame(res)
  df <- dplyr::arrange(df, p.adjust)
  df <- utils::head(df, top_n)
  df$Label <- ifelse(is.na(df$Description) | df$Description == "",
                     df$ID, df$Description)
  df$Label <- factor(df$Label, levels = rev(df$Label))
  df$GeneRatio <- vapply(df$GeneRatio, function(x) {
    parts <- as.numeric(strsplit(x, "/", fixed = TRUE)[[1]])
    if (length(parts) == 2 && parts[2] > 0) parts[1] / parts[2] else NA_real_
  }, numeric(1))

  ggplot(df, aes(x = GeneRatio, y = Label, size = Count, color = p.adjust)) +
    geom_point() +
    scale_color_gradient(low = "#e45756", high = "#4c78a8") +
    labs(title = label, x = "Gene Ratio", y = NULL,
         size = "Count", color = "p.adjust") +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank())
}

# ---- Run GO (by ontology) + KEGG ------------------------------------------
go_split <- split_go_by_ontology(go_t2g)
go_results <- lapply(names(go_split), function(ont) {
  t2g <- go_split[[ont]]
  if (is.null(t2g) || nrow(t2g) == 0) return(list(result = NULL, label = paste0("GO-", ont)))
  t2n <- go_term2name(unique(t2g$term))
  run_enrichment(t2g, t2n, paste0("GO-", ont))
})

kegg_t2n <- if (!is.null(kegg_t2g)) kegg_term2name(unique(kegg_t2g$term)) else NULL
kegg_result <- run_enrichment(kegg_t2g, kegg_t2n, "KEGG")

# ---- Render ---------------------------------------------------------------
plots <- c(
  lapply(go_results, function(x) make_dotplot(x$result, x$label, opt$top_n)),
  list(make_dotplot(kegg_result$result, kegg_result$label, opt$top_n))
)

pdf(pdf_path, width = 8, height = 7)
for (p in plots) print(p)
dev.off()

# ---- Tabular outputs (one CSV per category) -------------------------------
write_result_csv <- function(res, label) {
  if (is.null(res) || nrow(as.data.frame(res)) == 0) return(invisible(NULL))
  out <- file.path(opt$outdir,
                   sprintf("%s_%s.csv", opt$species, gsub("[^A-Za-z0-9]+", "_", label)))
  readr::write_csv(as.data.frame(res), out)
}
lapply(go_results, function(x) write_result_csv(x$result, x$label))
write_result_csv(kegg_result$result, "KEGG")

cat("[13_go_kegg] Wrote", pdf_path, "\n")
