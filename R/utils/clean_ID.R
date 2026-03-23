#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
args <- commandArgs(trailingOnly = TRUE)
read_delim(file = args[1], col_names = FALSE, delim = "\t", show_col_types = FALSE) %>%
  set_names("label") %>%
  mutate(label = str_remove(label, "\\.\\d+\\.*\\w*$")) %>%
  mutate(label = str_remove(label, "-\\w+$")) %>%
  mutate(label = str_remove(label, "_P\\d{3}$")) %>%
  write.table(file = paste0(args[1], ".clean.out"), quote = FALSE, row.names = FALSE, col.names = FALSE)
