rm(list = ls())

####----load R Package----####
library(tidyverse)
library(Biostrings)
library(ggfun)
library(ggraph)
library(tidygraph)
library(purrr)
library(igraph)
library(rlang)
library(patchwork)
library(readxl)
library(writexl)
source("R/get_node_edge.R")


####----load Data----####
df_group2 <- read_xlsx(path = "6.tree/tree_group_2.xlsx") %>%
  dplyr::filter(Species == "Sorghum_bicolor")

# load AtNet file
AtNet <- read_delim(file = "./11.ppi/AraNet.txt", delim = "\t", col_names = F) %>%
  rlang::set_names(c("source", "target", "weight"))

# load Glyma to AT
Glyma_to_AT <- read_delim(file = "./11.ppi/Sbicolor_Athaliana.blast", delim = "\t", col_names = F) %>%
  dplyr::select(1,2) %>%
  dplyr::mutate(X2 = str_remove(X2, pattern = "\\.\\d+$"))


Glyma_to_AT %>%
  dplyr::left_join(AtNet, by = c("X2" = "source")) %>%
  dplyr::filter(weight > 4) %>%
  purrr::set_names(c("Glyma_Source", "AT_Source", "Target", "Weight")) -> Glyma_to_AT_2


# load AT to Glyma blast
AT_to_glyma <- read_delim(file = "./11.ppi/Athaliana_Sbicolor.second.blast", delim = "\t", col_names = F) %>%
  dplyr::select(1,2) %>%
  purrr::set_names(c("AT_Target", "Glyma_Target")) %>%
  dplyr::mutate(AT_Target = str_remove(AT_Target, pattern = "\\.\\d+$"))

Glyma_to_AT_2 %>%
  dplyr::left_join(AT_to_glyma, by = c("Target" = "AT_Target")) %>%
  dplyr::filter(Weight > 4) %>%
  dplyr::select(1,5,4) %>%
  dplyr::filter( Glyma_Source != Glyma_Target) %>%
  purrr::set_names(c("Source", "Target", "Weight")) -> ppi_df

table(ppi_df$Source)
table(ppi_df$Target)

table(c(ppi_df$Source,
        ppi_df$Target)) %>%
  length()


ppi_df2 <- ppi_df %>%
  dplyr::mutate(Source = str_remove(Source, pattern = "\\.\\d$"),
                Target = str_remove(Target, pattern = "\\.\\d$"))


####----PPI Analysis----####
# 注释

unique(c(ppi_df$Source,
         ppi_df$Target)) %>%
  as.data.frame() %>%
  purrr::set_names("ID") %>%
  write.table(file = "11.ppi/PPI.geneid.txt", quote = F,row.names = F, col.names = F)

# 然后拿回去进行注释

####----Pfam----####
Pfam_annotation <- read.table(file = "11.ppi/all.ID.pep.select.Pfam_scan.out", comment.char = "#", header = F) %>%
  dplyr::left_join(df_group2 %>% dplyr::select(label, rename),
                   by = c("V1" = "label")) %>%
  dplyr::mutate(rename = ifelse(is.na(rename), V1, rename)) %>%
  dplyr::select(rename, V7) %>%
  dplyr::rename(domain = V7) %>%
  dplyr::group_by(rename) %>%
  dplyr::mutate(domain2  = str_c(domain, collapse = ",")) %>%
  dplyr::ungroup() %>%
  dplyr::select(1,3) %>%
  dplyr::distinct(rename, .keep_all = T) %>%
  dplyr::rename(domain = domain2)

 
edge <- ppi_df %>%
  dplyr::left_join(df_group2 %>% dplyr::select(label, rename),
                   by = c("Source" = "label")) %>%
  dplyr::rename(Source_rename = rename) %>%
  dplyr::left_join(df_group2 %>% dplyr::select(label, rename),
                   by = c("Target" = "label")) %>%
  dplyr::rename(Target_rename = rename) %>%
  dplyr::mutate(Source_rename = ifelse(is.na(Source_rename), Source, Source_rename),
                Target_rename = ifelse(is.na(Target_rename), Target, Target_rename)
                ) %>%
  dplyr::select(Source_rename, Target_rename, Weight) %>%
  purrr::set_names(c("from", "to", "weight"))


edge_annotation <-  ppi_df %>%
  dplyr::left_join(df_group2 %>% dplyr::select(label, rename),
                   by = c("Source" = "label")) %>%
  dplyr::rename(Source_rename = rename) %>%
  dplyr::left_join(df_group2 %>% dplyr::select(label, rename),
                   by = c("Target" = "label")) %>%
  dplyr::rename(Target_rename = rename) %>%
  dplyr::mutate(Source_rename = ifelse(is.na(Source_rename), Source, Source_rename),
                Target_rename = ifelse(is.na(Target_rename), Target, Target_rename)
  ) %>%
  dplyr::left_join(Pfam_annotation, by = c("Source_rename" = "rename")) %>%
  dplyr::rename(Source_domain = domain) %>%
  dplyr::left_join(Pfam_annotation, by = c("Target_rename" = "rename")) %>%
  dplyr::rename(Target_domian = domain)

write_xlsx(edge_annotation,
           path = "11.ppi/edge_annotation.xlsx")


####----Plot----####

ppi_obj <- get_node_edge("edge")

edge_out <- ppi_obj[[1]]
node_out <- ppi_obj[[2]]

graph_df <- as_tbl_graph(edge_out) %>%
  tidygraph::mutate(Popularity = centrality_degree(mode = 'out')) %>%
  tidygraph::left_join(Pfam_annotation, by = c("name" = "rename"))


p1 <- ggraph(graph_df, layout = 'linear', circular = TRUE) + 
  geom_edge_arc(aes(color = value), alpha = 0.5) + 
  geom_node_point(aes(size = Popularity), fill= "#dd3497",color = "#000000",alpha=1, shape = 21) + 
  scale_edge_colour_distiller(palette = "YiGn", direction = 1, name="Weight") + 
  coord_fixed() +
  scale_size(range = c(1,5)) +
  geom_node_text(
    aes(
      x = 1.05 * x,
      y = 1.05 * y,
      label = name,
      angle = -((-node_angle(x, y) + 90) %% 180) + 90,
    ),
    size = 2, hjust = 'outward'
  ) + 
  guides(color = guide_legend(order = 1),
         size = guide_legend(order = 2)) + 
  theme_nothing() +
  theme(
    plot.margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "cm"),
    plot.background = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank()
  ) + 
  coord_equal(xlim=c(-1.2,1.2),
              ylim=c(-1.2,1.2))

p1

ggsave(filename = "./11.ppi/PPI1.pdf", plot = p1, height = 10, width = 10.5)



# ####----Plot2----####
# 
# p2 <- ggraph(graph_df, layout = 'linear', circular = TRUE) + 
#   geom_edge_arc(aes(color = value), alpha = 0.5) + 
#   geom_node_point(aes(size = Popularity, fill = Domain),color = "#000000",alpha=1, shape = 21) + 
#   scale_edge_colour_distiller(palette = "RdPu", direction = 1, name="Weight") + 
#   scale_fill_manual(values = c("p450" = "#f768a1",
#                                "Terpene" = "#4eb3d3")) + 
#   coord_fixed() +
#   scale_size(range = c(5,15)) +
#   geom_node_text(
#     aes(
#       x = 1.1 * x,
#       y = 1.1 * y,
#       label = name,
#       angle = -((-node_angle(x, y) + 90) %% 180) + 90,
#     ),
#     size = 2.5, hjust = 'outward'
#   ) + 
#   guides(color = guide_legend(order = 1),
#          size = guide_legend(order = 2)) + 
#   theme_nothing() +
#   theme(
#     plot.margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "cm"),
#     plot.background = element_blank(),
#     panel.background = element_blank(),
#     panel.border = element_blank()
#   ) + 
#   coord_cartesian(xlim=c(-1.2,1.2),
#                   ylim=c(-1.2,1.2))
# 
# p2
# 
# ggsave(filename = "PPI2.pdf", plot = p2, height = 10, width = 10.5)


# edge_out %>%
#   dplyr::select(from, to) %>%
#   dplyr::left_join(Pfam_annotation, by = c("from" = "ID")) %>%
#   dplyr::select(from,Domain, to) %>%
#   dplyr::rename(from_Domain = Domain) %>%
#   dplyr::left_join(Pfam_annotation, by = c("to" = "ID")) %>%
#   dplyr::rename(to_Domain = Domain) %>%
#   write.csv(file = "PPI_Domain_annotation.csv",quote = F, row.names = F)
