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

####----load Data----####
bed <- read_delim(file = "5.genefamily_info/species_10.bed", col_names = F, delim = "\t") %>%
  dplyr::mutate(X4 = str_remove(X4, pattern = ".MSUv7.0")) %>%
  dplyr::mutate(X4 = str_remove(X4, pattern = "\\.v3\\.1"))

#####-----species3-----#####
# tree
tree <- treeio::read.newick(file = "6.tree/identify.ID.muscle.treefile", 
                            node.label = "support")

as_tibble(tree) %>%
  dplyr::filter(!is.na(label)) %>%
  dplyr::mutate(Species = case_when(
    str_starts(string = label, pattern = "AT") ~ "Arabidopsis thaliana",
    str_starts(string = label, pattern = "LOC") ~ "Oryza sativa Japonica",
    str_starts(string = label, pattern = "Sobic") ~ "Sorghum_bicolor"
  )) %>%
  dplyr::mutate(label2 = str_remove(string = label, pattern = "\\.\\d+$")) %>%
  dplyr::mutate(label2 = str_remove(string = label2, pattern = "\\.\\d\\.\\w$")) -> tree_Info

# tree info
AT_Info <- read_xlsx(path = "6.tree/AT_Os_Tree_Info.xlsx", col_names = T) %>%
  dplyr::select(1,2) %>%
  dplyr::filter(!is.na(`Gene Number`)) %>%
  dplyr::mutate(Clade = case_when(
    str_detect(`NPF name`, pattern = "NPF1") ~ "NPF1",
    str_detect(`NPF name`, pattern = "NPF2") ~ "NPF2",
    str_detect(`NPF name`, pattern = "NPF3") ~ "NPF3",
    str_detect(`NPF name`, pattern = "NPF4") ~ "NPF4",
    str_detect(`NPF name`, pattern = "NPF5") ~ "NPF5",
    str_detect(`NPF name`, pattern = "NPF6") ~ "NPF6",
    str_detect(`NPF name`, pattern = "NPF7") ~ "NPF7",
    str_detect(`NPF name`, pattern = "NPF8") ~ "NPF8",
  )) %>%
  dplyr::mutate(species = case_when(
    str_detect(`NPF name`, pattern = "At") ~ "Arabidopsis thaliana",
    str_detect(`NPF name`, pattern = "Os") ~ "Oryza sativa Japonica",
  )) %>%
  dplyr::mutate(`Gene Number` = case_when(
    species == "Arabidopsis thaliana" ~ str_to_upper(`Gene Number`),
    .default = as.character(`Gene Number`)
  )) %>%
  dplyr::mutate(`Gene Number` = str_replace(string = `Gene Number`, pattern = "LOC_OS", replacement = "LOC_Os"))

table(AT_Info$Clade,
      AT_Info$species)

spiece_info <- tree_Info %>% 
  dplyr::select(label, label2, Species) %>%
  dplyr::left_join(AT_Info, by = c("label2" = "Gene Number")) %>%
  dplyr::filter(!is.na(label))

# plot tmp
p_tmp <- ggtree(tree, branch.length = "none") %<+% spiece_info +
  geom_tiplab(aes(color = Clade), offset = 0.5) + 
  geom_tree2() + 
  xlim(NA, 50) + 
  geom_nodepoint(aes(fill=cut(support, c(0, 45, 75, 100))), 
                 shape=21, size=2) + 
  scale_fill_manual(values=c("black", "grey", "white"), 
                    guide='legend', name='Bootstrap Percentage(BP)', 
                    breaks=c('(75,100]', '(45,75]', '(0,45]'), 
                    labels=expression(BP>=75,45 <= BP * " < 75", BP < 45)) + 
  geom_nodelab(aes(label = node), size = 5) +
  new_scale_fill() + 
  geom_tippoint(shape = 21, aes(fill = Species), size = 3)

p_tmp

ggsave(filename = "6.tree/tmp.pdf",
       plot = p_tmp,
       height = 80,
       width = 20,
       limitsize = F)

p_tmp2 <- ggtree::rotate(p_tmp, node = 256) %>%
  ggtree::rotate(., node = 349) %>%
  ggtree::rotate(., node = 348) %>%
  ggtree::rotate(., node = 287) %>%
  ggtree::rotate(., node = 288)
 
ggsave(filename = "6.tree/tmp2.pdf",
       plot = p_tmp2,
       height = 80,
       width = 20,
       limitsize = F)

######------分组------######
get_taxa_name(p_tmp2) %>%
  as.data.frame() %>%
  purrr::set_names("ID") %>%
  dplyr::mutate(Group = case_when(
    # NPF1
    ID %in% c("LOC_Os01g55610.1", "Sobic.K022800.1.p") ~ "NPF1",
    # NPF2
    ID %in% c("LOC_Os12g44100.1", "AT3G45680.1") ~ "NPF2",
    # NPF3
    ID %in% c("Sobic.009G136600.1.p", "AT1G68570.1") ~ "NPF3",
    # NPF4
    ID %in% c("LOC_Os01g01360.1", "AT1G59740.1") ~ "NPF4",
    # NPF6
    ID %in% c("Sobic.001G541900.2.p", "AT3G21670.1") ~ "NPF6",
    # NPF8 
    ID %in% c("AT3G54140.1", "AT2G02040.1") ~ "NPF8",
    # NPF7
    ID %in% c("Sobic.001G282000.1.p", "Sobic.010G166850.1.p") ~ "NPF7",
    # NPF5
    ID %in% c("AT5G46040.1", "AT2G38100.1") ~ "NPF5"
  )) %>%
  tidyr::fill(Group) %>%
  dplyr::mutate(Species = case_when(
    str_starts(string = ID, pattern = "AT") ~ "Arabidopsis_thaliana",
    str_starts(string = ID, pattern = "LOC") ~ "Oryza_sativa_Japonica",
    str_starts(string = ID, pattern = "Sobic") ~ "Sorghum_bicolor"
  )) %>% 
  dplyr::mutate(label2 = str_remove(string = ID, pattern = "\\.\\d+\\.*\\w*$"))  %>%
  dplyr::mutate(label = ID) %>%
  dplyr::select(label, everything()) %>%
  dplyr::left_join(bed %>% dplyr::select(1,2,3,4),
                   by = c("label2" = "X4")) -> df_group

write_xlsx(df_group, path = "6.tree/tree_group.xlsx")

df_group %>%
  dplyr::arrange(Group, Species, X1, X2) %>%
  dplyr::group_by(Group, Species) %>%
  dplyr::mutate(n = 1:n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(species2 = case_when(
    str_starts(label2, pattern = "AT") ~ "At",
    str_starts(label2, pattern = "LOC") ~ "Os",
    str_starts(label2, pattern = "Sobic") ~ "Sb",
  )) %>%
  dplyr::mutate(rename = str_c(str_c(species2, Group, sep = ""),
                               n, sep = ".")) -> df_group2

# write_xlsx(df_group2, path = "6.tree/tree_group_2.xlsx")

df_group2 <- read_xlsx(path = "6.tree/tree_group_2.xlsx", col_names = T)


# 对文件进行批量的改名字


# 基于这个数据，我们继续绘图

p_plot <- ggtree(tree,
                 branch.length = "none",
                 # size = 0.0001,
                 color = "#969696",
                 layout = "circular"
)  %<+% df_group2 + 
  geom_tree(size = 0.0001,
            color = "#969696") + 
  geom_nodepoint(aes(fill=cut(support, c(0, 45, 75, 100))), 
                 shape=21, size=2, stroke = 0.2) + 
  scale_fill_manual(values=c("black", "grey", "white"), 
                    name='Bootstrap Percentage(BP)', 
                    breaks=c('(75,100]', '(45,75]', '(0,45]'), 
                    labels=expression(BP>=75,45 <= BP * " < 75", BP < 45),
                    guide = guide_legend(override.aes = list(size = 3))) + 
  geom_tiplab(aes(x = x + 0.5, label = rename, color = Group), size = 2) +
  new_scale_fill() + 
  geom_tippoint(aes(fill = Species), shape = 21) + 
  scale_fill_manual(values = c("Arabidopsis_thaliana" = "#de77ae",
                               "Sorghum_bicolor" = "#d6604d",
                               "Oryza_sativa_Japonica" = "#66bd63")) + 
  xlim(NA, 30) + 
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

p_plot

p_plot_rotate <- ggtree::rotate(p_plot, node = 256) %>%
  ggtree::rotate(., node = 349) %>%
  ggtree::rotate(., node = 348) %>%
  ggtree::rotate(., node = 287) %>%
  ggtree::rotate(., node = 288)

p_plot_rotate

p_plot_rotate + 
  # NPF1
  geom_strip("LOC_Os01g55610.1", "Sobic.K022800.1.p", 
             barsize=2, color='#78c679',offset=4, 
             label="NPF1",fontsize=4.5, offset.text=1, extend=0.1, alpha=0.2, angle = 35, hjust = 0.5) + 
  # NPF2
  geom_strip("LOC_Os12g44100.1", "AT3G45680.1",
             barsize=2, color='#238443',offset=4,
             label='NPF2',fontsize=4.5, offset.text=1, extend=0.1, alpha=0.2, angle = 70) +
  # NPF3
  geom_strip("Sobic.009G136600.1.p", "AT1G68570.1",
             barsize=2, color='#c2a5cf',offset=4,
             label='NPF3',fontsize=4.5, offset.text=1, extend=0.1, alpha=0.2, angle = 5, hjust = 0.5) +
  # NPF4
  geom_strip("LOC_Os01g01360.1", "AT1G59740.1",
             barsize=2, color='#88419d',offset=4,
             label='NPF4',fontsize=4.5, offset.text=1, extend=0.1, alpha=0.2, angle = -25) +
  # NPF5
  geom_strip("AT5G46040.1", "AT2G38100.1",
             barsize=2, color='#fec44f',offset=4,
             label='NPF5',fontsize=4.5, offset.text=1, extend=0.1, alpha=0.2, angle = -35) +
  # NPF6
  geom_strip("Sobic.001G541900.2.p", "AT3G21670.1",
             barsize=2, color='#cc4c02',offset=4,
             label='NPF6',fontsize=4.5, offset.text=1, extend=0.1, alpha=0.2, angle = -65, hjust = 0.5) +
  # NPF7
  geom_strip("Sobic.001G282000.1.p", "Sobic.010G166850.1.p",
             barsize=2, color='#7bccc4',offset=4,
             label='NPF7',fontsize=4.5, offset.text=1, extend=0.1, alpha=0.2, angle = 25) +
  # NPF8
  geom_strip("AT3G54140.1", "AT2G02040.1",
             barsize=2, color='#0868ac',offset=4,
             label='NPF8',fontsize=4.5, offset.text=1, extend=0.1, alpha=0.2, angle = 65) +
  # background
  # NPF1
  geom_strip("LOC_Os01g55610.1", "Sobic.K022800.1.p", 
             barsize=15, color='#78c679',offset=1.5, 
             label="",fontsize=4.5, offset.text=1, extend=0.5, alpha=0.2, angle = 35, hjust = 0.5) + 
  # NPF2
  geom_strip("LOC_Os12g44100.1", "AT3G45680.1",
             barsize=15, color='#238443',offset=1.5,
             label='',fontsize=4.5, offset.text=1, extend=0.5, alpha=0.2, angle = 70) +
  # NPF3
  geom_strip("Sobic.009G136600.1.p", "AT1G68570.1",
             barsize=15, color='#c2a5cf',offset=1.5,
             label='',fontsize=4.5, offset.text=1, extend=0.5, alpha=0.2, angle = 5, hjust = 0.5) +
  # NPF4
  geom_strip("LOC_Os01g01360.1", "AT1G59740.1",
             barsize=15, color='#88419d',offset=1.5,
             label='',fontsize=4.5, offset.text=1, extend=0.5, alpha=0.2, angle = -25) +
  # NPF5
  geom_strip("AT5G46040.1", "AT2G38100.1",
             barsize=15, color='#fec44f',offset=1.5,
             label='',fontsize=4.5, offset.text=1, extend=0.5, alpha=0.2, angle = -35) +
  # NPF6
  geom_strip("Sobic.001G541900.2.p", "AT3G21670.1",
             barsize=15, color='#cc4c02',offset=1.5,
             label='',fontsize=4.5, offset.text=1, extend=0.5, alpha=0.2, angle = -65, hjust = 0.5) +
  # NPF7
  geom_strip("Sobic.001G282000.1.p", "Sobic.010G166850.1.p",
             barsize=15, color='#7bccc4',offset=1.5,
             label='',fontsize=4.5, offset.text=1, extend=0.5, alpha=0.2, angle = 25) +
  # NPF8
  geom_strip("AT3G54140.1", "AT2G02040.1",
             barsize=15, color='#0868ac',offset=1.5,
             label='',fontsize=4.5, offset.text=1, extend=0.5, alpha=0.2, angle = 65) +
  # add label
  geom_tiplab(aes(x = x + 0.5, label = rename),color = "#000000", size = 2) + 
  theme(#legend.position = c(0.57, 0.46),
        legend.background = element_blank()) -> p_plot2

p_plot2

ggsave(filename = "6.tree/p_plot.pdf",
       plot = p_plot2,
       height = 15,
       width = 15)


#####-----sb-----#####
tree2 <- treeio::read.newick(file = "6.tree/Sb.NPF.ID.muscle.treefile", 
                             node.label = "support")

df_group2 <- read_xlsx(path = "6.tree/tree_group_2.xlsx") %>%
  dplyr::filter(Species == "Sorghum_bicolor")

p_plot2 <- ggtree(tree2,
                  branch.length = "none",
                  # size = 0.0001,
                  color = "#969696",
                  layout = "circular"
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
  geom_tiplab(aes(x = x + 0.5, label = rename, color = Group), size = 3) +
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
             barsize=2, color='#78c679',offset=4.5, 
             label="NPF1",fontsize=4.5, offset.text=1, extend=0.1, alpha=0.2, angle = 75, hjust = 0.5) + 
  # NPF2
  geom_strip("Sobic.009G099000.1.p", "Sobic.003G399200.1.p",
             barsize=2, color='#238443',offset=4.5,
             label='NPF2',fontsize=4.5, offset.text=1, extend=0.1, alpha=0.2, angle = -70) +
  # NPF3
  geom_strip("Sobic.009G136600.1.p", "Sobic.010G109200.1.p",
             barsize=2, color='#c2a5cf',offset=4.5,
             label='NPF3',fontsize=4.5, offset.text=1, extend=0.1, alpha=0.2, angle = -40, hjust = 0.5) +
  # NPF4
  geom_strip("Sobic.007G081300.1.p", "Sobic.003G107100.1.p",
             barsize=2, color='#88419d',offset=4.5,
             label='NPF4',fontsize=4.5, offset.text=1, extend=0.1, alpha=0.2, angle = 10) +
  # NPF5
  geom_strip("Sobic.009G148900.1.p", "Sobic.009G101800.2.p",
             barsize=2, color='#fec44f',offset=4.5,
             label='NPF5',fontsize=4.5, offset.text=1, extend=0.1, alpha=0.2, angle = -70) +
  # NPF6
  geom_strip("Sobic.004G193000.1.p", "Sobic.001G541900.2.p",
             barsize=2, color='#cc4c02',offset=4.5,
             label='NPF6',fontsize=4.5, offset.text=1, extend=0.1, alpha=0.2, angle = 50, hjust = 0.5) +
  # NPF7
  geom_strip("Sobic.001G282100.1.p", "Sobic.010G133100.1.p",
             barsize=2, color='#7bccc4',offset=4.5,
             label='NPF7',fontsize=4.5, offset.text=1, extend=0.1, alpha=0.2, angle = -10) +
  # NPF8
  geom_strip("Sobic.003G392700.1.p", "Sobic.001G111150.1.p",
             barsize=2, color='#0868ac',offset=4.5,
             label='NPF8',fontsize=4.5, offset.text=1, extend=0.1, alpha=0.2, angle = 35)  + 
  theme(#legend.position = c(0.57, 0.46),
    legend.background = element_blank()) -> p_plot2_3

p_plot2_3


ggsave(filename = "6.tree/p_plot2.pdf",
       plot = p_plot2_3,
       height = 12.5,
       width = 12.5)
####----sessionInfo----####
sessionInfo()
