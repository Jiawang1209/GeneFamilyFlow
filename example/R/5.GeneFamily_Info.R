rm(list = ls())

####----load R Pakcage----####
library(tidyverse)
library(Biostrings)
library(Peptides)
library(ggfun)
library(ggbeeswarm)
library(patchwork)
library(ggprism)
library(ggtree)
library(tidytree)
library(treeio)
library(aplot)
library(patchwork)
library(showtext)
font_families()
showtext_auto()

####----load Data----####
pep.fa <- readAAStringSet(filepath = "5.genefamily_info/identify.ID.clean.fa")

all_bed <- read_delim(file = "5.genefamily_info/species_10.bed", col_names = F, delim = "\t") %>%
  purrr::set_names(c("Chr","Start","End","ID","Info","Strand")) %>%
  dplyr::mutate(label2 = str_remove(string = ID, pattern = "\\.\\d+$")) %>%
  dplyr::mutate(label2 = str_remove(string = label2, pattern = ".MSU.*")) %>%
  dplyr::mutate(label2 = str_remove(string = label2, pattern = ".v3$"))


pep_info <- data.frame(pep.fa) %>%
  tibble::rownames_to_column(var = "ID") %>%
  dplyr::mutate(Length = Peptides::lengthpep(seq = pep.fa)) %>%  # lengthpep() 计算长度
  dplyr::mutate(MW = mw(seq = pep.fa)) %>%                 # mw() 计算分子量
  dplyr::mutate(hydrophobicity = hydrophobicity(seq = pep.fa)) %>%     # hydrophobicity() 计算疏水性
  dplyr::mutate(pI = pI(seq = pep.fa)) %>%                 # pI() 计算等电点
  as_tibble() %>% 
  dplyr::mutate(Species = case_when(
    str_starts(string = ID, pattern = "AT") ~ "Arabidopsis_thaliana",
    str_starts(string = ID, pattern = "LOC") ~ "Oryza_sativa_Japonica",
    str_starts(string = ID, pattern = "Sobic") ~ "Sorghum_bicolor"
  )) %>%
  dplyr::mutate(label2 = str_remove(string = ID, pattern = "\\.\\d+$")) %>%
  dplyr::mutate(label2 = str_remove(string = label2, pattern = ".MSU.*")) %>%
  dplyr::mutate(label2 = str_remove(string = label2, pattern = "\\.\\d\\.p")) %>%
  dplyr::left_join(all_bed, by = c("label2" = "label2"))

pep_info2 <- pep_info %>%
  dplyr::select(1,8:11,13,14,2:7)


writexl::write_xlsx(pep_info2, path = "5.genefamily_info/Gene_Information.xlsx")



p_stat <- pep_info2 %>%
  dplyr::select(Species, Length, MW, hydrophobicity, pI) %>%
  tidyr::pivot_longer(cols = -Species, values_to = "Value", names_to = "Kind") %>%
  dplyr::mutate(Kind = factor(Kind, levels = c("Length","MW","hydrophobicity","pI"), ordered = T)) %>%
  ggplot() + 
  geom_boxplot(aes(x = Species, 
                    y = Value, 
                    fill = Species), alpha = 0.75) + 
  geom_quasirandom(aes(x = Species, 
                       y = Value, 
                       fill = Species, 
                       group = Species),
                   shape = 21, 
                   size = 3, 
                   alpha = 0.7,
                   varwidth = TRUE) + 
  scale_fill_manual(values = c("Arabidopsis_thaliana" = "#4eb3d3",
                               "Sorghum_bicolor" = "#dd3497",
                               "Oryza_sativa_Japonica" = "#807dba")) + 
  facet_wrap(~Kind, scales = "free", nrow = 2) + 
  labs(x = "", y = "") + 
  theme_bw() + 
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "#000000", size = 10),
        # axis.text.x = element_text(face = "italic"),
        panel.border = element_rect(linewidth = 1),
        strip.text = element_text(color = "#000000", size = 15),
        strip.background = element_rect(linewidth = 1),
        # axis.text.y = element_blank(),
        # axis.ticks.y = element_blank()
  ) 

p_stat


ggsave(filename = "5.genefamily_info/species_gene_info.pdf",
       plot = p_stat,
       height = 7.5,
       width = 12.5)
