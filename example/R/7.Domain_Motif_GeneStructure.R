rm(list = ls())

####----load R Package----####
library(tidyverse)
library(ggtree)
library(treeio)
library(tidytree)
library(aplot)
library(ggnewscale)
library(ggfun)
library(readxl)
library(writexl)
library(ggalluvial)
library(patchwork)
library(gggenes)
source("R/gene_structure_data.R")

####----Tree----####
tree2 <- treeio::read.newick(file = "6.tree/Sb.NPF.ID.muscle.treefile", 
                             node.label = "support")

df_group2 <- read_xlsx(path = "6.tree/tree_group_2.xlsx") %>%
  dplyr::filter(Species == "Sorghum_bicolor")

p_plot2 <- ggtree(tree2,
                  branch.length = "none",
                  # size = 0.0001,
                  color = "#969696",
                  # layout = "circular"
)  %<+% df_group2 + 
  geom_tree(size = 0.0001,
            color = "#969696") + 
  geom_nodepoint(aes(fill=cut(support, c(0, 45, 75, 100))), 
                 shape=21, size=3, stroke = 0.2) + 
  scale_fill_manual(values=c("black", "grey", "white"), 
                    name='Bootstrap Percentage(BP)', 
                    breaks=c('(75,100]', '(45,75]', '(0,45]'), 
                    labels=expression(BP>=75,45 <= BP * " < 75", BP < 45),
                    guide = guide_legend(override.aes = list(size = 3))) + 
  geom_tiplab(aes(x = x + 0.25, label = rename, color = Group), size = 3) +
  new_scale_fill() + 
  geom_tippoint(aes(fill = Species), shape = 21, size = 3) + 
  scale_fill_manual(values = c("Arabidopsis_thaliana" = "#de77ae",
                               "Sorghum_bicolor" = "#d6604d",
                               "Oryza_sativa_Japonica" = "#66bd63")) + 
  # geom_nodelab(aes(label = node), color = "red") +
  xlim(NA, 25) + 
  scale_color_manual(values = c("NPF1" = "#78c679",
                                "NPF2" = "#238443",
                                "NPF3" = "#c2a5cf",
                                "NPF4" = "#88419d",
                                "NPF5" = "#fec44f",
                                "NPF6" = "#cc4c02",
                                "NPF7" = "#7bccc4",
                                "NPF8" = "#0868ac"
  ),
  name = "Subfamily")

p_plot2

p_plot2_2 <- p_plot2 %>%
  ggtree::rotate(., node = 137) %>%
  ggtree::rotate(., node = 138) %>%
  ggtree::rotate(., node = 96) %>%
  ggtree::rotate(., node = 91) %>%
  ggtree::rotate(., node = 92)

p_plot2_2


get_taxa_name(p_plot2_2) %>%
  as.data.frame() %>%
  purrr::set_names("Rename") %>%
  dplyr::left_join(df_group2, by = c("Rename" = "label")) -> df_group2_2

p_plot2_2 + 
  # NPF1
  geom_strip("Sobic.003G304100.1.p", "Sobic.K022800.1.p", 
             barsize=2, color='#78c679',offset=3, 
             label="NPF1",fontsize=4.5, offset.text=0.5, extend=0.1, alpha=0.2) + 
  # NPF2
  geom_strip("Sobic.009G099000.1.p", "Sobic.003G399200.1.p",
             barsize=2, color='#238443',offset=3,
             label='NPF2',fontsize=4.5, offset.text=0.5, extend=0.1, alpha=0.2) +
  # NPF3
  geom_strip("Sobic.009G136600.1.p", "Sobic.010G109200.1.p",
             barsize=2, color='#c2a5cf',offset=3,
             label='NPF3',fontsize=4.5, offset.text=0.5, extend=0.1, alpha=0.2) +
  # NPF4
  geom_strip("Sobic.007G081300.1.p", "Sobic.003G107100.1.p",
             barsize=2, color='#88419d',offset=3,
             label='NPF4',fontsize=4.5, offset.text=0.5, extend=0.1, alpha=0.2) +
  # NPF5
  geom_strip("Sobic.009G148900.1.p", "Sobic.009G101800.2.p",
             barsize=2, color='#fec44f',offset=3,
             label='NPF5',fontsize=4.5, offset.text=0.5, extend=0.1, alpha=0.2) +
  # NPF6
  geom_strip("Sobic.004G193000.1.p", "Sobic.001G541900.2.p",
             barsize=2, color='#cc4c02',offset=3,
             label='NPF6',fontsize=4.5, offset.text=0.5, extend=0.1, alpha=0.2) +
  # NPF7
  geom_strip("Sobic.001G282100.1.p", "Sobic.010G133100.1.p",
             barsize=2, color='#7bccc4',offset=3,
             label='NPF7',fontsize=4.5, offset.text=0.5, extend=0.1, alpha=0.2) +
  # NPF8
  geom_strip("Sobic.003G392700.1.p", "Sobic.001G111150.1.p",
             barsize=2, color='#0868ac',offset=3,
             label='NPF8',fontsize=4.5, offset.text=0.5, extend=0.1, alpha=0.2)  + 
  theme(#legend.position = c(0.57, 0.46),
    legend.background = element_blank()) -> p_plot2_3

p_plot2_3

tree.sort.id <- get_taxa_name(p_plot2_3)


####----Domain----####
Domain <- read.table(file = "4.identification/Pfam_scan.out", comment = "#", header = F) %>%
  dplyr::filter(V1 %in% tree.sort.id) 


# 二者数据是一样的，那么按照一个去分析即可
Domain1_plot <- Domain %>%
  dplyr::mutate(V1 = factor(V1, levels = rev(tree.sort.id), ordered = T))

domain_p <- ggplot(Domain1_plot, aes(xmin = V4, xmax = V5, y = V1, fill = V7)) +
  geom_gene_arrow(
    arrowhead_height = unit(3, "mm"),
    arrowhead_width = unit(0.1,"mm")
  ) +
  labs(x = "", y = "") +
  ggtitle(label = "Domain Analysis") +
  scale_fill_manual(values = c('#fa9fb5','#1f78b4','#b2df8a','#33a02c','#fb9a99',
                               '#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a',
                               '#ffff99','#b15928','#ffffb3','#bebada','#fb8072'),
                    name = "Domain")+
  theme_genes() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    # legend.background = element_roundrect(color = "#969696"),
    plot.title = element_text(hjust = 0.5)
  )

domain_p


domain_p %>% insert_left(., p_plot2_3)


####----Motif----####
Motif <- read_delim(file = "7.motif_genestructure/meme_location.txt", col_names = F, delim = "\t")

setdiff(tree.sort.id, unique(Motif$X1))


Motif_plot <- Motif %>%
  dplyr::left_join(df_group2 %>% dplyr::select(1,2,3,4,11), by = c("X1" = "ID")) %>%
  dplyr::mutate(X1 = factor(X1,
                                levels = rev(tree.sort.id),
                                ordered = T)) %>%
  dplyr::mutate(X4 = factor(X4, levels = paste("MEME-", 1:15, sep = ""), ordered = T))

# motif
motif_p <- ggplot(Motif_plot, aes(xmin = X2, xmax = X3, y = X1, fill = X4)) +
  geom_gene_arrow(
    arrowhead_height = unit(3, "mm"),
    arrowhead_width = unit(0.1,"mm")
  ) +
  labs(x = "", y = "") + 
  ggtitle(label = "Motif Analysis") + 
  scale_fill_manual(values = c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99',
                               '#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a',
                               '#ffff99','#b15928','#ffffb3','#bebada','#fb8072'),
                    name = "Motif")+
  theme_genes() + 
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    # legend.background = element_roundrect(color = "#969696"),
    plot.title = element_text(hjust = 0.5)
  )

motif_p

motif_p %>% 
  insert_left(domain_p) %>%
  insert_left(p_plot2_3)


####----CDS----####
gff <- read_delim(file = "1.database/Sbicolor.gff3", delim = "\t", col_names = F, skip = 2, comment = "#") %>%
  dplyr::filter(X3 %in% c("CDS")) %>%
  tidyr::separate(col = "X9", sep = ";", into = c("tmp1"), remove = F) %>%
  dplyr::mutate(ID = str_remove(tmp1, pattern = "ID=")) %>%
  dplyr::mutate(ID2 = ID) %>%
  dplyr::mutate(ID3 = str_remove(ID, pattern = ".v3.1.CDS.\\d+$")) %>%
  dplyr::mutate(ID3 = str_c(ID3, ".p", sep = "")) %>%
  dplyr::filter(ID3 %in% tree.sort.id)

gff2 <- gff %>%
  dplyr::mutate(delta = abs(X5-X4)) %>%
  dplyr::group_by(ID3) %>%
  dplyr::summarise(sum = sum(delta))

gff_clean <- gff %>%
  dplyr::filter(ID3 %in% gff2$ID3) %>%
  dplyr::select(-c(9,10,11,12)) %>%
  purrr::set_names(c("seqnames","Orginal","Type","start","end",".","strand","V8","ID")) 

# gff_clean_sub <- gff_clean %>% dplyr::filter(ID %in% unique(tree.sort.id))



gene_structure_out <- gene_structure_data(gffobject = "gff_clean")

gene_structure_sub <- gene_structure_out %>%
  dplyr::mutate(ID = factor(ID, levels = rev(tree.sort.id), ordered = T))

gene_structure_p <- ggplot(gene_structure_sub, aes(xmin = start, xmax = end, y = ID, fill = Type)) +
  geom_gene_arrow(
    arrowhead_height = unit(3, "mm"),
    arrowhead_width = unit(0.1,"mm")
  ) +
  labs(x = "", y = "") + 
  ggtitle(label = "Gene Structure") + 
  scale_fill_manual(values = c('#7fcdbb','#1f78b4','#b2df8a','#33a02c','#fb9a99',
                               '#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a',
                               '#ffff99','#b15928','#ffffb3','#bebada','#fb8072'),
                    name = "Gene Structure")+
  theme_genes() + 
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    # legend.background = element_roundrect(color = "#969696"),
    plot.title = element_text(hjust = 0.5)
  )

gene_structure_p

####----Combine----####
p_combine <- gene_structure_p %>% 
  insert_left(motif_p, width = 0.75) %>%
  insert_left(domain_p, width = 0.75) %>% 
  insert_left(p_plot2_3, width = 1.5)

# p_combine

ggsave(filename = "./7.motif_genestructure/Tree_Domain_Motif_GeneStructure.pdf",
       plot = p_combine,
       height = 22.5,
       width = 20)


####----sessionInfo----####
sessionInfo()
