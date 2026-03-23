rm(list = ls())

####----load R Package----####
library(tidyverse)
library(readxl)
library(pheatmap)
library(RColorBrewer)

####----load Data----####
df_group2 <- read_xlsx(path = "6.tree/tree_group_2.xlsx") %>%
  dplyr::filter(Species == "Sorghum_bicolor")

ID <- df_group2 %>%
  dplyr::select(label2) %>%
  purrr::set_names("ID")

####----RNA-seq1----####

expr1 <- read_delim(file = "12.rnaseq/表达量分析结果表_Exp_G_RSEM_TPM_20240721_141710.csv",
                    col_names = T,
                    delim = ",") %>%
  dplyr::filter(`Gene Name` %in% ID$ID) %>%
  dplyr::select(-c(1,3)) %>%
  dplyr::left_join(df_group2 %>% dplyr::select(label2, rename), by = c("Gene Name" = "label2")) %>%
  dplyr::select(-`Gene Name`) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(sum = sum(c_across(where(is.numeric)))) %>%
  dplyr::filter(sum != 0) %>%
  dplyr::select(-sum) %>%
  tibble::column_to_rownames(var = "rename") 
  

pdf(file = "12.rnaseq/heatmap1.pdf",
    height = 16.5,
    width = 7)

pheatmap(expr1, cluster_cols = F, scale = "row",
         color = colorRampPalette(brewer.pal(11, "PiYG"))(25))

dev.off()


####----组织表达----####
df_group2_2 <- df_group2 %>%
  dplyr::select(label2, rename)

expr2 <- read_xlsx(path = "12.rnaseq/工作簿.xlsx", col_names = T) %>%
  dplyr::select(-1) %>%
  dplyr::rename(ID = `BTx623 Homolog`) %>%
  dplyr::filter(!is.na(ID)) %>%
  dplyr::left_join(df_group2_2, by = c("ID" = "label2")) %>%
  dplyr::select(-ID) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(Leaf = mean(c(`Leaf 1`, `Leaf 2`, `Leaf 3`, `Leaf 4`)),
                Stem = mean(c(`Stem 1`, `Stem 2`)),
                Inflorescence = mean(c(`Inflorescence 1`, `Inflorescence 2`)),
                Seed = mean(c(`Seed 1`, `Seed 2`, `Seed 3`))
                ) %>%
  dplyr::ungroup() %>%
  dplyr::select(rename, Root, Stem,Leaf,Inflorescence, Seedling,Seed) %>%
  tibble::column_to_rownames(var = "rename")




pdf(file = "12.rnaseq/heatmap2.pdf",
    height = 16.5,
    width = 7)

pheatmap(expr2, cluster_cols = F, scale = "row",
         color = colorRampPalette(brewer.pal(11, "PiYG"))(25))

dev.off()

