rm(list = ls())

####----load R Package----####
library(tidyverse)
library(ggtree)
library(treeio)
library(tidytree)
library(aplot)
library(ggfun)
library(readxl)
library(ggh4x)
library(patchwork)
library(ggnewscale)

####----load Data----####
#####-----tree file-----#####
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



#####-----promoter-----#####
cir_element.desc <- read_xlsx(path = "./10.promoter/cir_element.desc.20240509.xlsx")

list_filt_v <- list.files(path = "10.promoter", recursive = T, pattern = ".tab", full.names = T)

promoter <- data.table::fread(file = list_filt_v,
                        header = F,
                        sep = "\t",
                        quote="") %>%
  dplyr::select(1,2,8) %>%
  purrr::set_names("X1","X2","X3") %>%
  dplyr::filter(!str_detect(X2, pattern = "Unname")) %>%
  dplyr::filter(!is.na(X2)) %>%
  dplyr::filter(!str_detect(X3, "short_function")) %>%
  dplyr::left_join(cir_element.desc, by = c("X2" = "element")) %>%
  dplyr::filter(!is.na(description)) %>%
  dplyr::filter(str_detect(X1, pattern = "Sobic")) %>%
  dplyr::mutate(X1 = str_remove(X1, pattern = "^.*usf:2000"))

setdiff(tree.sort.id %>%
          str_remove(., pattern = "\\.\\d\\.p"),
        promoter$X1 %>% 
          str_remove(., pattern = ".v3.1") %>%
          unique())


remove_circ <- promoter %>%
  pull(X2) %>%
  table() %>%
  as_tibble() %>%
  purrr::set_names("element","count") %>%
  dplyr::arrange(desc(count)) %>%
  left_join(cir_element.desc, 
            by = c("element" = "element")) %>%
  dplyr::filter(count < 5) %>%
  dplyr::pull(element)

# contain
promoter %>%
  pull(X2) %>%
  table() %>%
  as_tibble() %>%
  purrr::set_names("element","count") %>%
  dplyr::arrange(desc(count)) %>%
  left_join(cir_element.desc, 
            by = c("element" = "element")) %>%
  dplyr::filter(!element %in% remove_circ) %>%
  dplyr::arrange(description, desc(count)) %>% 
  group_by(description) %>%
  top_n(5, count) %>%
  pull(element) -> retain_circ_element


retain_circ_element



# promoter
promoter %>%
  dplyr::select(1,2,4) %>%
  dplyr::mutate(X1 = str_remove(X1, pattern = ".v3.1")) %>%
  dplyr::left_join(df_group2 %>% dplyr::select(1,2,3,4,5,11), by = c("X1" = "label2")) %>%
  dplyr::select(5,2,3) %>%
  dplyr::rename(X1 = ID) %>%
  dplyr::filter(X2 %in% retain_circ_element) %>%
  dplyr::group_by(X1, X2) %>%
  dplyr::summarise(count = n()) %>%
  dplyr::ungroup() %>% 
  tidyr::spread(key = X2, value = count) %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  tidyr::gather(key = "X2", value = "count", -X1) %>%
  left_join(cir_element.desc, by = c("X2" = "element")) %>%
  dplyr::mutate(X2 = factor(X2,levels = retain_circ_element,ordered = T)) %>%
  dplyr::mutate(X1 = factor(X1,levels = rev(tree.sort.id),ordered = T)) %>%
  na.omit() %>%
  ggplot(aes(interaction(X2,description), y = X1)) +
  geom_tile(fill = "white",color = "grey80",size = 0.5, height = 1, width = 0.5) +
  geom_point(aes(fill = cut(count, c(0,5,10,20,40,60,80,100,120,140,160,180,200),right = F)),
             shape = 22,size = 7) + 
  geom_text(aes(label = count),size = 3) + 
  scale_fill_manual(values = c("#fde0ef","#f1b6da","#d9f0a3",
                               "#addd8e","#78c679","#7fbc41","#41ab5d",
                               "#238443","#006837","#004529",
                               "#6a3d9a","#ffff99")) + 
  labs(fill = "Count") + 
  guides(x="axis_nested") +
  labs( x = "", y = "") + 
  # coord_fixed() + 
  theme_bw() + 
  theme(axis.text.x=element_text(color="black",angle=90,size=10,hjust= 0),
        axis.text.x.top = element_text(hjust = 0, vjust = 0.5),
        axis.text.y=element_text(color="black",size=12),
        # axis.text.y=element_blank(),
        legend.background = element_roundrect(color = "#808080", linetype = 1),
        axis.ticks.length.x.top = unit(0.25, "cm"),
        axis.ticks.y=element_blank(),
        panel.border=element_blank(),
        # panel.grid = element_blank(),
        ggh4x.axis.nestline.x = element_line(size = 1),
        ggh4x.axis.nesttext.x = element_text(colour = "#d73027",angle =0, hjust = 0.5, size = 10),
        # plot.margin=unit(c(0.2,0.2,0.2,0.2),units=,"cm")
        ) + 
  scale_x_discrete(expand = expansion(mult = c(0.01, 0.01)),position = "top") + 
  guides(fill = guide_legend(direction = "vertical")) + 
  guides(fill = guide_legend(ncol = 1, bycol = F)) -> p.promoter

p.promoter

p.promoter2 = p.promoter + theme(axis.text.y = element_blank())

promoter %>%
  dplyr::select(1,2,4) %>%
  dplyr::mutate(X1 = str_remove(X1, pattern = ".v3.1")) %>%
  dplyr::left_join(df_group2 %>% dplyr::select(1,2,3,4,5,11), by = c("X1" = "label2")) %>%
  dplyr::select(5,2,3) %>%
  dplyr::rename(X1 = ID) %>%
  dplyr::filter(X2 %in% retain_circ_element) %>%
  dplyr::group_by(X1, X2) %>%
  dplyr::summarise(count = n()) %>%
  dplyr::ungroup() %>% 
  tidyr::spread(key = X2, value = count) %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  tidyr::gather(key = "X2", value = "count", -X1) %>%
  left_join(cir_element.desc, by = c("X2" = "element")) %>%
  dplyr::mutate(X2 = factor(X2,levels = retain_circ_element,ordered = T)) %>%
  dplyr::mutate(X1 = factor(X1,levels = rev(tree.sort.id),ordered = T)) %>%
  na.omit() %>%
  write.csv(file = "./10.promoter/Promote.stat.csv")


# 统计
stat_promoter <-promoter %>%
  dplyr::select(1,2,4) %>%
  dplyr::mutate(X1 = str_remove(X1, pattern = ".v3.1")) %>%
  dplyr::left_join(df_group2 %>% dplyr::select(1,2,3,4,5,11), by = c("X1" = "label2")) %>%
  dplyr::select(5,2,3) %>%
  dplyr::rename(X1 = ID) %>%
  dplyr::filter(X2 %in% retain_circ_element) %>%
  dplyr::group_by(X1, X2) %>%
  dplyr::summarise(count = n()) %>%
  dplyr::ungroup() %>% 
  tidyr::spread(key = X2, value = count) %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  tidyr::gather(key = "X2", value = "count", -X1) %>%
  left_join(cir_element.desc, by = c("X2" = "element")) %>%
  dplyr::mutate(X2 = factor(X2,levels = retain_circ_element,ordered = T)) %>%
  dplyr::mutate(X1 = factor(X1,levels = rev(tree.sort.id),ordered = T)) %>%
  na.omit()  %>%
  dplyr::group_by(X1, description) %>%
  dplyr::summarise(stat = sum(count)) %>%
  dplyr::ungroup() %>%
  ggplot(aes(x = stat, y = X1, fill = description)) + 
  geom_bar(stat = "identity", color = "black", alpha = 0.75, position = "fill", width = 0.5) + 
  labs(x = "", y = "", fill = "Description") + 
  scale_x_continuous(expand = expansion(mult = c(0, 0)),
                     breaks = c(0, 0.25, 0.5, 0.75, 1.00),
                     labels = c("0%", "25%", "50%", "75%", "100%")) + 
  scale_fill_manual(values = c( "#de77ae","#fdb863","#35978f","#d6604d")) + 
  theme_bw() + 
  theme(axis.text.x=element_text(color="black",size=15),
        legend.background = element_roundrect(color = "#808080", linetype = 1),
        axis.ticks.length.x.top = unit(0.25, "cm"),
        axis.ticks.y=element_blank(),
        panel.border = element_rect(linewidth = 1.25),
        # plot.margin=unit(c(0.2,0.2,0.2,0.2),units=,"cm"),
        axis.text.y = element_blank())

stat_promoter


####----Plot----####
p_combine <- stat_promoter %>% 
  insert_left(p_plot2_3, width = 1.5) %>% 
  insert_right(p.promoter2, width = 1.75)

# p_combine

ggsave(filename = "./10.promoter/tree_promoter_stat.pdf",
       plot = p_combine,
       height = 25,
       width = 22)

