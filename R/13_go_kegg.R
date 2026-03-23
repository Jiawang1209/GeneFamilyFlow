#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

df <- data.frame(
  term = c("GO:0008150", "GO:0003674", "ko00010", "ko04075"),
  count = c(12, 7, 5, 4),
  group = c("BP", "MF", "KEGG", "KEGG")
)

p <- ggplot(df, aes(x = count, y = term, fill = group)) +
  geom_col() +
  theme_bw() +
  labs(title = "GO/KEGG placeholder")

ggsave("13_go_kegg_placeholder.pdf", p, width = 8, height = 6)
