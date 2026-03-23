rm(list = ls())
library(tidyverse)
library(circlize)
library(ComplexHeatmap)
library(readxl)

dir.create("9.mcscanx/Sorghum_bicolor", recursive = T)

####---- load ChrTAome file
Glyma_Chr <- read_delim(file = './9.mcscanx/Sbicolor_454_v3.0.1.fa.fai',
                     delim = '\t', col_names = F)  %>%
  dplyr::select(1,2) %>%
  dplyr::slice(1:10) %>%
  dplyr::mutate(tmp = 1) %>%
  dplyr::select(1,3,2) %>%
  purrr::set_names(c("X1","X2","X3")) %>% 
  dplyr::mutate(X1 = str_c('Sb',X1))

####---- load bed file
Glyma_bed <- read_delim(file = './10.promoter/species_10.bed', 
                     col_names = F, delim = '\t') %>%
  dplyr::select(1:4) %>%
  dplyr::filter(str_detect(X4, pattern = "Sob")) %>%
  dplyr::mutate(X4 = str_remove(X4, pattern = ".v3.1")) %>%
  dplyr::mutate(X1 = str_c('Sb',X1))

####----AT_TPS ID
df_group2 <- read_xlsx(path = "6.tree/tree_group_2.xlsx") %>%
  dplyr::select(label, rename)

Glyma_TPS.id <- read.table(file = "./9.mcscanx/Sbicolor.NPF.id", sep = " ", header = T) %>%
  set_names("ID") %>% 
  dplyr::left_join(df_group2, by = c("ID" = "label")) %>%
  dplyr::mutate(ID = str_remove(ID, pattern = "\\.\\d\\.p"))

Glyma_TPSs_bed <- left_join(Glyma_TPS.id, Glyma_bed, by = c('ID' = 'X4')) %>%
  dplyr::select(3:5,1,2) %>%
  set_names(c("Chr","Start","End","ID","Rename"))

####---- load gene type file
Glyma_genetype <- read_delim(file = './9.mcscanx/Sbicolor.gene_type',
                          col_names = F, delim = '\t') %>%
  dplyr::mutate(X1 = str_remove(X1, pattern = "\\.\\d+$"))

####---- tandem
TPSs_tandem <- read_delim(file = "./9.mcscanx/Sbicolor.tandem",
                          col_names = F, delim = ",") %>%
  dplyr::mutate(X1 = str_remove(X1, pattern = "\\.\\d+$")) %>%
  dplyr::mutate(X2 = str_remove(X2, pattern = "\\.\\d+$")) %>%
  dplyr::filter(X1 %in% Glyma_TPS.id$ID | X2 %in% Glyma_TPS.id$ID)

####---- collinearity
TPSs_collinearity <- read_delim(file = "./9.mcscanx/Sbicolor.collinearity",
                                col_names = F, delim = "\t", comment = "#")  %>%
  dplyr::select(X2, X3) %>%
  purrr::set_names(c("X1", "X2")) %>%
  dplyr::mutate(X1 = str_remove(X1, pattern = "\\.\\d+$")) %>%
  dplyr::mutate(X2 = str_remove(X2, pattern = "\\.\\d+$")) %>%
  dplyr::filter(X1 %in% Glyma_TPS.id$ID | X2 %in% Glyma_TPS.id$ID)

rbind(TPSs_tandem %>% dplyr::mutate(Group = "Tandem"),
      TPSs_collinearity %>% dplyr::mutate(Group = "Collinearity")) %>%
  write.csv(file = "./9.mcscanx/Sorghum_bicolor/Sbicolor_genepair.csv", quote = F, row.names = F)

unique(c(TPSs_tandem$X1,TPSs_tandem$X2,TPSs_collinearity$X1, TPSs_collinearity$X2)) %>%
  as.data.frame() %>%
  purrr::set_names("ID") %>%
  write.csv(file = "./9.mcscanx/Sorghum_bicolor/Sbicolor_genepair.ID.csv", quote = F, row.names = F)

#对基因组进行滑动窗口 500kb为一个窗口
Chr_No = c()
Chr_start = c()
Chr_end = c()

for (i in seq_along(Glyma_Chr$X1)) {
  print(i)
  print(Glyma_Chr$X1[i])
  
  chr_length <- Glyma_Chr[i,] %>% as.character() %>% .[3] %>% as.numeric()
  start = seq(0, chr_length, by = 500000)
  end = start + 500000
  windows_length = length(start)
  Chr_no <- rep(Glyma_Chr$X1[i] ,windows_length)
  
  Chr_No <- c(Chr_No, Chr_no)
  Chr_start <- c(Chr_start, start)
  Chr_end <- c(Chr_end, end)
  
}

Chr_window <- data.frame(
  Chr = Chr_No,
  Start = Chr_start,
  End = Chr_end
)

#对TPSs基因进行统计 统计窗口内部有多少个基因 分染色体去找是一个好的方法
WIN_NUM <- NULL
for (Chr_tmp in unique(Chr_No)) {
  # Chr_tmp = 'Chr1'
  print(Chr_tmp)
  #基于染色体，筛选基因组窗口文件
  Chr_window_sub <- Chr_window[Chr_window$Chr == Chr_tmp,]
  #筛选染色体，筛选TPSs文件
  Glyma_TPSs_bed_sub <- Glyma_TPSs_bed %>% dplyr::filter(Chr == Chr_tmp)
  
  if (dim(Glyma_TPSs_bed_sub)[1] == 0) {
    WIN_NUM <- c(WIN_NUM, rep(0, dim(Chr_window_sub)[1]))
  }else{
    for (i in 1:dim(Chr_window_sub)[1]) {
      win_num = 0
      Chr_window_sub[i,]$Start
      Chr_window_sub[i,]$End
      for (j in 1:dim(Glyma_TPSs_bed_sub)[1]) {
        if (Glyma_TPSs_bed_sub[j,]$Start > Chr_window_sub[i,]$Start &
            Glyma_TPSs_bed_sub[j,]$End < Chr_window_sub[i,]$End) {
          win_num <- win_num + 1
        }
      }
      WIN_NUM <- c(WIN_NUM, win_num)
    }
  }
}

Chr_window$Number_Count <- WIN_NUM
range(Chr_window$Number_Count)

#对TPSs基因的复制事件类型进行统计
Glyma_TPS.id %>% 
  dplyr::select(ID) %>%
  left_join(Glyma_genetype,
            by = c('ID' = 'X1')) %>%
  left_join(Glyma_TPSs_bed, by=c('ID' = 'ID')) %>%
  dplyr::select(3:5, 1,2) %>%
  set_names(c("Chr","Start","End","ID","Type"))  %>%
  dplyr::mutate(Type = case_when(
    Type == 0 ~ 0,
    Type == 1 ~ 1,
    Type == 2 ~ 2,
    Type == 3 ~ 3,
    Type == 4 ~ 4,
  )) -> gene_type



#tandem bed
TPSs_tandem_bed1 <- TPSs_tandem %>% 
  dplyr::select(1) %>%
  dplyr::left_join(Glyma_bed, by = c('X1' = 'X4')) %>%
  dplyr::select(2,3,4,1) %>%
  dplyr::filter(X1.y !='Un')
TPSs_tandem_bed2 <- TPSs_tandem %>% 
  dplyr::select(2) %>%
  dplyr::left_join(Glyma_bed, by = c('X2' = 'X4')) %>%
  dplyr::select(2,3,4,1) %>%
  dplyr::filter(X1 !='Un')

#WGD bed
TPSs_collin_bed1 <- TPSs_collinearity %>% 
  dplyr::filter(!c(str_detect(X1, pattern = 'U')|str_detect(X2,pattern = 'U'))) %>%
  dplyr::select(1) %>%
  dplyr::left_join(Glyma_bed, by = c('X1' = 'X4')) %>%
  dplyr::select(2,3,4,1) 

TPSs_collin_bed2 <- TPSs_collinearity %>% 
  dplyr::filter(!c(str_detect(X1, pattern = 'U')|str_detect(X2,pattern = 'U'))) %>%
  dplyr::select(2) %>%
  dplyr::left_join(Glyma_bed, by = c('X2' = 'X4')) %>%
  dplyr::select(2,3,4,1)

pdf(file = './9.mcscanx/Sorghum_bicolor/circos_Sbicolor.pdf',
    height = 8,
    width = 8)
####---- 开始画染色体
circos.genomicInitialize(Glyma_Chr, 
                         plotType = NULL, 
                         axis.labels.cex = 0.4*par("cex"),
                         labels.cex = 0.6*par("cex"),
                         track.height = 0.01, #设置的轨道的宽度
                         major.by = 20000000 #设置刻度大小5Mb为一个主刻度
)

# ####----添加label 因为这个label含有700多个基因，太多了，所以不能先添加它
circos.genomicLabels(Glyma_TPSs_bed %>% dplyr::filter(Chr != "Sbsuper_10"),
                     labels.column = 5,
                     padding = 0.1,
                     connection_height = mm_h(3),
                     col = as.numeric(factor(Glyma_TPSs_bed[[1]])),
                     line_col = as.numeric(factor(Glyma_TPSs_bed[[1]])),
                     cex = 0.85,
                     side = "outside")


####---- 染色体填充颜色
circos.genomicTrackPlotRegion(
  Glyma_Chr, track.height = 0.05, stack = TRUE, bg.border = NA,
  panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value, col = '#dd3497', border = 'black', ...)
  } )

# # ####----添加label 因为这个label含有700多个基因，太多了，所以不能先添加它
# circos.genomicLabels(AT_TPSs_bed[sample(1:700,100),],
#                      labels.column = 4,
#                      side = 'inside')

####----添加TPSs滑动窗口的统计
circos.genomicTrack(
  Chr_window, 
  track.height = 0.1, 
  bg.col = '#f0f0f0',
  bg.border = NA,
  panel.fun = function(region, value, ...){
    circos.genomicLines(region, value, col='blue', lwd=0.35,...)
    circos.lines(c(0, 0.5, 1),
                 c(0, 0.5, 1),
                 col = 'blue2',
                 lwd = 0.15,
                 lty = 2)
    circos.yaxis(labels.cex = 0.2,
                 lwd = 0.1,
                 tick.length = convert_x(0.15, 'mm'))
  }
)

####----添加复制事件
color_assign <- colorRamp2(breaks = c(0,1,2,3,4), 
                           col = c("#00ADFF", "#e66101","#fdb863", "#b2abd2", "#5e3c99"))

circos.genomicTrackPlotRegion(
  gene_type, track.height = 0.1, stack = TRUE, bg.border = NA,
  panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value, col = color_assign(value[[1]]),
                       border = color_assign(value[[1]]), ...)
  } )


####----添加连线Link
colors_tmp <- scales::alpha('black', alpha = 1)
#colline
circos.genomicLink(TPSs_collin_bed1 %>% dplyr::slice(1:9),
                   TPSs_collin_bed2 %>% dplyr::slice(1:9), 
                   col = colors_tmp,
                   border = '#66c2a4',
                   lwd = 1)
#tandem
circos.genomicLink(TPSs_tandem_bed1, TPSs_tandem_bed2,
                   col = colors_tmp,
                   border = '#fec44f',
                   lwd = 1)

gene_legend <- Legend(
  at = c(0, 1, 2, 3, 4), 
  labels = c("Singleton","Dispersed", "Proximal", "Tandem", "WGD/segmental"), 
  labels_gp = gpar(fontsize = 8),
  title = "Gene Type", 
  title_gp = gpar(fontsize = 9), 
  grid_height = unit(0.4, "cm"), 
  grid_width = unit(0.4, "cm"), 
  type = "points", pch = NA, 
  background = c("#00ADFF", "#e66101","#fdb863", "#b2abd2", "#5e3c99"))

Duplicate_legend <- Legend(
  at = c(1, 2), 
  labels = c('Collinearity', 'Tandem'),
  labels_gp = gpar(fontsize = 8), 
  grid_height = unit(0.5, 'cm'), 
  grid_width = unit(0.5, 'cm'), 
  type = c('lines', 'lines'), 
  pch = NA, 
  legend_gp = gpar(col = c('#66c2a4', '#fec44f'), lwd = 1),
  title = 'Gene Duplicate', 
  title_gp = gpar(fontsize = 9))

y_coord <- 0.9
x_coord <- 0.9

pushViewport(viewport(x = x_coord, y = y_coord))
grid.draw(gene_legend)
y_coord <- y_coord
upViewport()

pushViewport(viewport(x = x_coord, y = y_coord - 0.75))
grid.draw(Duplicate_legend)
y_coord <- y_coord 
upViewport()

####----clear
circos.clear()

dev.off()
