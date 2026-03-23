rm(list = ls())

####----load R Package----####
library(tidyverse)
library(AnnotationForge)
library(clusterProfiler)
library(magrittr)
library(org.Sbean.eg.db, lib = "R_Library")
library(ggh4x)
library(readxl)
library(writexl)
library(ggfun)

####----load eggnogmapper----####

eggno <- read_delim(file = "13.GO_KEGG/MM_0ln9gnah.emapper.annotations.tsv", 
                    col_names = T,
                    delim = "\t",
                    comment = "##")

# annotation
emapper <- eggno %>%
  dplyr::select(
    GID = `#query`,
    GO = GOs,
    KO = KEGG_ko,
    Pathway = KEGG_Pathway
  )

#gene_info
gene_info <- emapper %>%
  dplyr::select(GID) %>%
  dplyr::mutate(Gene_Name = GID)

#gene2go
gene2go <- emapper %>%
  dplyr::select(GID, GO) %>%
  separate_rows(GO, sep = ",", convert = F) %>%
  dplyr::filter(!is.na(GO)) %>%
  dplyr::filter(str_detect(GO, pattern = "GO")) %>%
  dplyr::mutate(EVIDENCE = "IEA")

####----create package----####
AnnotationForge::makeOrgPackage(gene_info=gene_info,
                                go=gene2go,
                                maintainer='RPython <yueliu1115@163.com>',
                                author='RPython',
                                outputDir=".",
                                tax_id=20250712,
                                genus='S.',
                                species='bean',
                                goTable="go",
                                version="2.0")

dir.create("R_annotation_Package")

pkgbuild::build('./org.Sbean.eg.db/', dest_path = "R_annotation_Package")

dir.create("R_Library", recursive = T)

install.packages('./R_annotation_Package/org.Sbean.eg.db_2.0.tar.gz',
                 repos = NULL,
                 lib = "R_Library")

####----enrichment----####
df_group2 <- read_xlsx(path = "6.tree/tree_group_2.xlsx") %>%
  dplyr::filter(Species == "Glycine max")

ID <- df_group2 %>%
  dplyr::select(label2) %>%
  purrr::set_names("ID")


# 获取KEGG
get_path2name <- function(){
  keggpathid2name.df <- clusterProfiler:::kegg_list(db = "pathway")
  keggpathid2name.df <- keggpathid2name.df %>%
    purrr::set_names("path_id","path_name") %>%
    dplyr::mutate(path_id = str_replace(path_id, pattern="map", replacement = "ko")) # %>%
  return(keggpathid2name.df)
}
pathway2name <- get_path2name()


eggno <- read_delim(file = "13.GO_KEGG/MM_0ln9gnah.emapper.annotations.tsv", 
                    col_names = T,
                    delim = "\t",
                    comment = "##")

emapper <- eggno %>%
  dplyr::select(
    GID = `#query`,
    GO = GOs,
    KO = KEGG_ko,
    Pathway = KEGG_Pathway
  )

pathway2gene <- dplyr::select(emapper, Pathway, GID) %>%
  separate_rows(Pathway, sep = ",", convert = F) %>%
  filter(str_detect(Pathway, "ko"))

####----Enrichment Analysys----####
GO_out <- enrichGO(gene = ID$ID,
                   OrgDb = org.Sbean.eg.db,
                   keyType = 'GID',
                   ont = 'ALL')


GO_out_df <- as.data.frame(GO_out)

write_xlsx(GO_out_df, path = "13.GO_KEGG/GO_out_df.xlsx")

KEGG_out <- enricher(ID$ID,
                     TERM2GENE = pathway2gene,
                     TERM2NAME = pathway2name,
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.05)
KEGG_out_df <- as.data.frame(KEGG_out)

write_xlsx(KEGG_out_df, path = "13.GO_KEGG/KEGG_out_df.xlsx")

####----Plot----####
plot_df <- read_xlsx(path = "13.GO_KEGG/Plot.xlsx") %>%
  tidyr::separate(col = GeneRatio, sep = "/", into = c("n1", "n2")) %>%
  dplyr::mutate(GeneRatio = as.numeric(n1)/as.numeric(n2)) %>%
  dplyr::mutate(ONTOLOGY = factor(ONTOLOGY, levels = rev(c("BP", "CC", "MF", "KEGG")), ordered = T)) %>%
  dplyr::arrange(ONTOLOGY, Count) %>%
  dplyr::mutate(Description = factor(Description, levels = Description, ordered = T))


plot <- plot_df %>%
  ggplot() + 
  geom_point(data = plot_df %>% dplyr::filter(ONTOLOGY == "KEGG"),
             aes(x = Count, y = interaction(Description, ONTOLOGY), fill = p.adjust, size = Count), shape = 21) + 
  scale_fill_gradient(low = "#a1d99b", high = "#238b45", name = "KEGG p.adjust",
                      guide = guide_colorbar(order = 4)) + 
  ggnewscale::new_scale_fill() + 
  geom_point(data = plot_df %>% dplyr::filter(ONTOLOGY == "MF"),
             aes(x = Count, y = interaction(Description, ONTOLOGY), fill = p.adjust, size = Count), shape = 21) + 
  scale_fill_gradient(low = "#a6bddb", high = "#0570b0", name = "MF p.adjust",
                      guide = guide_colorbar(order = 3)) + 
  ggnewscale::new_scale_fill() + 
  geom_point(data = plot_df %>% dplyr::filter(ONTOLOGY == "CC"),
             aes(x = Count, y = interaction(Description, ONTOLOGY), fill = p.adjust, size = Count), shape = 21) + 
  scale_fill_gradient(low = "#fdd49e", high = "#d7301f", name = "CC p.adjust",
                      guide = guide_colorbar(order = 2)) + 
  ggnewscale::new_scale_fill() +
  geom_point(data = plot_df %>% dplyr::filter(ONTOLOGY == "BP"),
             aes(x = Count, y = interaction(Description, ONTOLOGY), fill = p.adjust, size = Count), shape = 21) + 
  scale_fill_gradient(low = "#8c96c6", high = "#8c6bb1",  name = "BP p.adjust",
                      guide = guide_colorbar(order = 1)) + 
  guides(y = "axis_nested",
         y.sec = guide_axis_manual(breaks = 1:26,
                                   labels = plot_df$Description)) + 
  ggtitle(label = "GO and KEGG annotation") + 
  labs(x = "Count", y = "Description") + 
  scale_size(range= c(3, 7),
             guide = guide_legend(override.aes = list(fill = "#000000"))) + 
  theme_bw() + 
  theme(
    ggh4x.axis.nestline.y = element_line(size = 3, color = c("#74c476", "#41b6c4", "#f46d43", "#9e9ac8")),
    ggh4x.axis.nesttext.y = element_text(colour = c("#74c476", "#41b6c4", "#f46d43", "#9e9ac8")),
    # legend.background = element_roundrect(color = "#969696"),
    panel.border = element_rect(size = 0.5),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm"),
    axis.text = element_text(color = "#000000", size = 11),
    axis.text.y = element_text(color = c(rep("#41ae76",2), rep("#225ea8",6), rep("#fc4e2a",3), rep("#88419d",15))
    ),
    axis.text.y.left = element_blank(),
    axis.ticks.length.y.left = unit(10, "pt"),
    axis.ticks.y.left = element_line(color = NA),
    axis.title = element_text(color = "#000000", size = 15),
    plot.title = element_text(color = "#000000", size = 20, hjust = 0.5)
  ) + 
  coord_cartesian(clip = "off") 

plot

ggsave(filename = "./13.GO_KEGG/plot.pdf",
       plot = plot,
       height = 10,
       width = 11)


####----sessionInfo----####
sessionInfo()
