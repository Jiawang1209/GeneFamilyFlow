#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

df <- data.frame(
  gene = rep(c("SbNPF1.1", "SbNPF2.1"), each = 4),
  group = rep(c("ck_shoot", "ln_shoot", "ck_root", "ln_root"), times = 2),
  value = c(1.0, 1.8, 0.9, 1.6, 1.2, 2.1, 0.8, 1.4)
)

p <- ggplot(df, aes(x = group, y = value, fill = group)) +
  geom_col() +
  facet_wrap(~gene, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "qRT-PCR placeholder")

ggsave("14_qrt_pcr_placeholder.pdf", p, width = 10, height = 6)
