rm(list = ls())

library(tidyverse)
library(readxl)
library(writexl)

annotation <- read_xlsx(path = "manuscript_file_annotation/tree_group_2.xlsx", 
                        col_names = T) %>%
  dplyr::arrange(label, Group)

annotation2 <- read_xlsx(path = "5.genefamily_info/Gene_Information.xlsx",
                         col_names = T) %>%
  dplyr::select(1,8:13)

Info <- annotation %>%
  dplyr::left_join(annotation2, by = c("label" = "ID.x")) %>%
  dplyr::select(-label2) %>%
  dplyr::rename(Chr = X1,
                start = X2,
                end = X3) %>%
  dplyr::select(-n)

write_xlsx(Info, path = "manuscript_file_annotation/Gene_Info.xlsx")

Info2 <- Info %>%
  dplyr::filter(str_starts(label, pattern = "So")) %>%
  dplyr::select(label, rename) %>%
  dplyr::mutate(label = str_remove(label, pattern = "\\.\\d\\.\\w$"))

duplicate1 <- read_xlsx(path = "manuscript_file_annotation/Sbicolor_genepair_gene_duplicate.xlsx",
                        sheet = 1) %>%
  dplyr::left_join(Info2, by = c("X1" = "label")) %>%
  dplyr::select(X1, rename, X2, Group) %>%
  dplyr::left_join(Info2, by = c("X2" = "label")) %>%
  dplyr::select(1,2,3,5,4)

duplicate2 <- read_xlsx(path = "manuscript_file_annotation/Sbicolor_genepair_gene_duplicate.xlsx",
                        sheet = 2) %>%
  dplyr::mutate(Seq_1 = str_remove(Seq_1, pattern = "\\.\\d+$"),
                Seq_2 = str_remove(Seq_2, pattern = "\\.\\d+$")) %>%
  dplyr::left_join(Info2, by = c("Seq_1" = "label")) %>%
  dplyr::select(1,6, 2:5) %>%
  dplyr::left_join(Info2, by = c("Seq_2" = "label")) %>%
  dplyr::select(1,2,3,7,4:6)

write_xlsx(duplicate1, path = "manuscript_file_annotation/gene_duplication_20250914.xlsx")
write_xlsx(duplicate2, path = "manuscript_file_annotation/gene_duplication_2_20250914.xlsx")
