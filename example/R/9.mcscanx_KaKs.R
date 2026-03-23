rm(list = ls())

####----load R Package----####
library(tidyverse)
library(ggfun)
library(ggbeeswarm)
library(readxl)
library(Biostrings)

####----load Data----####
file_v <- list.files(path = "9.mcscanx/", pattern = "_genepair.csv", recursive = T, full.names = T)

map(file_v, function(x){read_delim(x, delim = ",", col_names = T)}) %>%
  do.call(rbind,.) -> file_df

file_df %>%
  dplyr::select(X1, X2) %>%
  write.table(file = "9.mcscanx/mcscanx_Ka_Ks.pair.txt",
              sep = "\t",
              quote = F,
              row.names = F,
              col.names = F)

file_df2 <- file_df %>%
  purrr::set_names(c("ID1", "ID2", "Type")) %>%
  dplyr::select(ID1, ID2, Type) %>%
  purrr::set_names(c("Seq_1", "Seq_2", "Duplcate Type")) %>%
  dplyr::mutate(`Duplcate Type` = str_replace(`Duplcate Type`, pattern = "tandem", "Tandem"))

write.table(file_df2,
            file = "9.mcscanx/mcscanx_Ka_Ks.pair2.txt",
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = F)

# 然后去进行kaks的分析，当然了，在此之前要清理一下序列

####----KaKs----####
MCScanX_kaks <- read_delim(file = "9.mcscanx/kaks.tab.xls", 
                           col_names = T, delim = "\t") %>%
  dplyr::filter(Seq_1 != "Seq_1") %>% 
  dplyr::select(1,2,3,4,5) %>%
  dplyr::filter(Ka_Ks != "NaN") %>%
  dplyr::mutate(Species = case_when(
    str_starts(string = Seq_1, pattern = "AT") ~ "Arabidopsis_thaliana",
    str_starts(string = Seq_1, pattern = "LOC") ~ "Oryza_sativa_Japonica",
    str_starts(string = Seq_1, pattern = "Sobic") ~ "Sorghum_bicolor"
  )) %>%
  dplyr::mutate(tmp = str_c(Seq_1, Seq_2)) %>%
  dplyr::left_join(file_df2 %>% dplyr::mutate(tmp = str_c(Seq_1, Seq_2)),
                   by = c("tmp" =  "tmp")) %>%
  dplyr::distinct(tmp, .keep_all = T) %>%
  dplyr::select(-tmp) %>%
  dplyr::select(Ka, Ks, Ka_Ks, Species, `Duplcate Type`) %>%
  tidyr::pivot_longer(cols = c(Ka, Ks, Ka_Ks), values_to = "Value", names_to = "Type") %>%
  dplyr::mutate(Value = as.numeric(Value)) %>%
  dplyr::mutate(Type = factor(Type, levels = c("Ka", "Ks", "Ka_Ks"), ordered = T))

####----Plot----####

ggplot(data = MCScanX_kaks, aes(x = Species, y = Value, fill = Species)) +
  geom_boxplot(outlier.colour = NA, alpha = 0.65, width = 0.6) + 
  geom_quasirandom(aes(fill = Species,color = Species, shape = Type),size = 2, alpha = 0.5,
                   varwidth = TRUE, color = "#000000") + 
  geom_hline(yintercept = 1, linetype = 2) + 
  facet_grid(`Duplcate Type`~Type, scales = "free") + 
  scale_fill_manual(values = c("Arabidopsis_thaliana" = "#4eb3d3",
                                        "Sorghum_bicolor" = "#dd3497",
                                        "Oryza_sativa_Japonica" = "#807dba")) + 
  scale_color_manual(values = c("Arabidopsis_thaliana" = "#4eb3d3",
                                "Sorghum_bicolor" = "#dd3497",
                                "Oryza_sativa_Japonica" = "#807dba")) + 
  scale_shape_manual(values = 21:23) + 
  labs(x = "", y = "") + 
  theme_bw() + 
  theme(panel.grid = element_blank(),
        panel.border = element_rect(linewidth = 1),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12, color = "#000000"),
        strip.text = element_text(size = 15),
        # axis.text.x = element_text(angle = 45, size = 10, hjust = 1, color = "#000000"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(5,5,5,20, unit = "pt")
        )

ggsave(filename = "9.mcscanx/9.mcscanx_Kaks.pdf",
       height = 7.5,
       width = 16)
  

####----sessionInfo----####
sessionInfo()