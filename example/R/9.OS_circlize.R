rm(list = ls())
library(tidyverse)
library(circlize)
library(ComplexHeatmap)

dir.create("9.mcscanx/Os", recursive = T)

####---- load Chrosome file
OS_Chr <- read_delim(file = './9.mcscanx/Osativa_323_v7.0.fa.fai',
                     delim = '\t', col_names = F) %>%
  dplyr::select(1,2) %>%
  dplyr::slice(1:12) %>%
  dplyr::mutate(tmp = 1) %>%
  dplyr::select(1,3,2) %>%
  purrr::set_names(c("X1","X2","X3")) %>% 
  dplyr::mutate(X1 = str_c('Os',X1))
####---- load bed file
OS_bed <- read_delim(file = '9.mcscanx/Osativa.bed', 
                     col_names = F, delim = '\t') %>%
  dplyr::select(1:4) %>%
  dplyr::mutate(X1 = str_c('Os',X1)) %>%
  dplyr::mutate(X4 = str_remove(X4, pattern = ".MSUv7.0"))

####----AT_HSP ID
df_group2 <- read_xlsx(path = "6.tree/tree_group_2.xlsx") %>%
  dplyr::select(label, rename)

OS_TPS.id <- read.table(file = "./9.mcscanx/Os.NPF.id", sep = " ", header = T) %>%
  set_names("ID") %>%
  dplyr::left_join(df_group2, by = c("ID" = "label"))

OS_TPS_bed <- left_join(OS_TPS.id, OS_bed, by = c('ID' = 'X4')) %>%
  dplyr::select(3:5,1,2) %>%
  set_names(c("Chr","Start","End","ID","Rename"))

####---- load gene type file
OS_genetype <- read_delim(file = './9.mcscanx/Osativa.gene_type',
                                 col_names = F, delim = '\t')

####---- tandem
TPS_tandem <- read_delim(file = "./9.mcscanx/Osativa.tandem",
                          col_names = F, delim = ",") %>%
  dplyr::filter(X1 %in% OS_TPS.id$ID | X2 %in% OS_TPS.id$ID)

####---- collinearity
TPS_collinearity <- read_delim(file = "./9.mcscanx/Osativa.collinearity",
                                col_names = F, delim = "\t",comment = "#")  %>%
  dplyr::select(X2, X3) %>%
  purrr::set_names(c("X1", "X2")) %>%
  dplyr::filter(X1 %in% OS_TPS.id$ID | X2 %in% OS_TPS.id$ID)


rbind(TPS_tandem %>% dplyr::mutate(Group = "Tandem"),
      TPS_collinearity %>% dplyr::mutate(Group = "Collinearity")) %>%
  write.csv(file = "./9.mcscanx/Os/OS_genepair.csv", quote = F, row.names = F)

unique(c(TPS_tandem$X1,TPS_tandem$X2,TPS_collinearity$X1, TPS_collinearity$X2)) %>%
  as.data.frame() %>%
  purrr::set_names("ID") %>%
  write.csv(file = "./9.mcscanx/Os/OS_genepair.ID.csv", quote = F, row.names = F)

#对基因组进行滑动窗口 500kb为一个窗口
Chr_No = c()
Chr_start = c()
Chr_end = c()

for (i in seq_along(OS_Chr$X1 )) {
  print(i)
  print(OS_Chr$X1[i])
  
  chr_length <- OS_Chr[i,] %>% as.character() %>% .[3] %>% as.numeric()
  start = seq(0, chr_length, by = 500000)
  end = start + 500000
  windows_length = length(start)
  Chr_no <- rep(OS_Chr$X1[i] ,windows_length)
  
  Chr_No <- c(Chr_No, Chr_no)
  Chr_start <- c(Chr_start, start)
  Chr_end <- c(Chr_end, end)
  
}

Chr_window <- data.frame(
  Chr = Chr_No,
  Start = Chr_start,
  End = Chr_end
)

#对CRGs基因进行统计 统计窗口内部有多少个基因 分染色体去找是一个好的方法
WIN_NUM <- NULL
for (Chr_tmp in paste0('OsChr',1:12)) {
  # Chr_tmp = 'OsChr6'
  print(Chr_tmp)
  #基于染色体，筛选基因组窗口文件
  Chr_window_sub <- Chr_window[Chr_window$Chr == Chr_tmp,]
  #筛选染色体，筛选CRGs文件
  OS_HSPs_bed_sub <- OS_TPS_bed %>% dplyr::filter(Chr == Chr_tmp)
  
  if (dim(OS_HSPs_bed_sub)[1] == 0) {
    WIN_NUM <- c(WIN_NUM, rep(0, dim(Chr_window_sub)[1]))
  }else if (dim(OS_HSPs_bed_sub)[1] != 0) {
    for (i in 1:dim(Chr_window_sub)[1]) {
      win_num = 0
      Chr_window_sub[i,]$Start
      Chr_window_sub[i,]$End
      for (j in 1:dim(OS_HSPs_bed_sub)[1]) {
        if (OS_HSPs_bed_sub[j,]$Start > Chr_window_sub[i,]$Start &
            OS_HSPs_bed_sub[j,]$End < Chr_window_sub[i,]$End) {
          win_num <- win_num + 1
        }else{
          win_num <- win_num
        }
      }
      WIN_NUM <- c(WIN_NUM, win_num)
    }
  }
  
  
}

Chr_window$Number_Count <- WIN_NUM
range(Chr_window$Number_Count)
# mean(Chr_window$Number_Count) -> mean_count

#对CRGs基因的复制事件类型进行统计
OS_TPS.id %>% 
  # dplyr::select(-Type) %>%
  left_join(OS_genetype,
            by = c('ID' = 'X1')) %>%
  left_join(OS_TPS_bed, by=c('ID' = 'ID')) %>%
  dplyr::select(4,5,6, 1,3) %>%
  set_names(c("Chr","Start","End","ID","Type")) %>%
  dplyr::mutate(Type = case_when(
    Type == 0 ~ 0,
    Type == 1 ~ 1,
    Type == 2 ~ 2,
    Type == 3 ~ 3,
    Type == 4 ~ 4,
  )) -> gene_type



#tandem bed
TPS_tandem_bed1 <- TPS_tandem %>%
  dplyr::select(1) %>%
  dplyr::left_join(OS_bed, by = c('X1' = 'X4')) %>%
  dplyr::select(2,3,4,1)
TPS_tandem_bed2 <- TPS_tandem %>%
  dplyr::select(2) %>%
  dplyr::left_join(OS_bed, by = c('X2' = 'X4')) %>%
  dplyr::select(2,3,4,1)

#WGD bed
TPS_collin_bed1 <- TPS_collinearity %>% 
  dplyr::select(1) %>%
  dplyr::left_join(OS_bed, by = c('X1' = 'X4')) %>%
  dplyr::select(2,3,4,1)
TPS_collin_bed2 <- TPS_collinearity %>% 
  dplyr::select(2) %>%
  dplyr::left_join(OS_bed, by = c('X2' = 'X4')) %>%
  dplyr::select(2,3,4,1)

pdf(file = './9.mcscanx/Os/circos_Os.pdf',
    height = 8,
    width = 8)
####---- 开始画染色体
circos.genomicInitialize(OS_Chr, 
                         plotType = NULL, 
                         axis.labels.cex = 0.4*par("cex"),
                         labels.cex = 0.6*par("cex"),
                         track.height = 0.01, #设置的轨道的宽度
                         major.by = 20000000 #设置刻度大小5Mb为一个主刻度
)

# ####----添加label 因为这个label含有700多个基因，太多了，所以不能先添加它
circos.genomicLabels(OS_TPS_bed,
                     labels.column = 5,
                     padding = 0.1,
                     connection_height = mm_h(3),
                     col = as.numeric(factor(OS_TPS_bed[[1]])),
                     line_col = as.numeric(factor(OS_TPS_bed[[1]])),
                     cex = 0.85,
                     side = "outside")

####---- 染色体填充颜色
circos.genomicTrackPlotRegion(
  OS_Chr, track.height = 0.05, stack = TRUE, bg.border = NA,
  panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value, col = '#807dba', border = 'black', ...)
  } )

####----添加CRGs滑动窗口的统计
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
circos.genomicLink(TPS_collin_bed1, TPS_collin_bed2, 
                   col = colors_tmp,
                   border = '#66c2a4',
                   lwd = 1)
#tandem
circos.genomicLink(TPS_tandem_bed1, TPS_tandem_bed2,
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

# ######################
# pdf(file = './OS_HSPs/circos_202209013_outside.pdf',
#     height = 5,
#     width = 5)
# ####---- 开始画染色体
# circos.genomicInitialize(OS_Chr, 
#                          # plotType = c('axis', 'labels'), 
#                          plotType = NULL, 
#                          axis.labels.cex = 0.4*par("cex"),
#                          labels.cex = 0.6*par("cex"),
#                          track.height = 0.05, #设置的轨道的宽度
#                          major.by = 20000000 #设置刻度大小20Mb为一个主刻度
# )
# 
# # ####----添加label 
# circos.genomicLabels(OS_HSPs_bed,
#                      labels.column = 4,
#                      padding = 0.4,
#                      connection_height = mm_h(3),
#                      col = as.numeric(factor(OS_HSPs_bed[[1]])),
#                      line_col = as.numeric(factor(OS_HSPs_bed[[1]])),
#                      cex = 0.5,
#                      side = 'outside')
# 
# 
# ####---- 染色体填充颜色
# circos.genomicTrackPlotRegion(
#   OS_Chr, track.height = 0.05, stack = TRUE, bg.border = NA,
#   panel.fun = function(region, value, ...) {
#     circos.genomicRect(region, value, col = '#238b45', border = 'black', ...)
#   } )
# circos.track(track.index = get.current.track.index(), 
#              bg.border = NA,
#              panel.fun = function(x, y) {
#                circos.genomicAxis(h = 'top',
#                                   direction = 'outside',
#                                   labels.cex = 0.2*par("cex"),
#                                   major.by = 20000000)
#              })
# 
# ####----添加CRGs滑动窗口的统计
# circos.genomicTrack(
#   Chr_window, 
#   track.height = 0.1, 
#   bg.col = '#f0f0f0',
#   bg.border = NA,
#   panel.fun = function(region, value, ...){
#     circos.genomicLines(region, value, col='blue', lwd=0.35,...)
#     circos.lines(c(0, 1, 2, 3, 4, 5),
#                  c(0, 1, 2, 3, 4, 5),
#                  col = 'blue2',
#                  lwd = 0.15,
#                  lty = 2)
#     circos.yaxis(labels.cex = 0.2, 
#                  lwd = 0.1,
#                  tick.length = convert_x(0.15, 'mm'))
#   }
# )
# 
# ####----添加复制事件
# color_assign <- colorRamp2(breaks = c(0,1,2,3), 
#                            col = c('#fb8072', '#80b1d3',
#                                    '#fdb462', '#b3de69'))
# 
# circos.genomicTrackPlotRegion(
#   gene_type, track.height = 0.1, stack = TRUE, bg.border = NA,
#   panel.fun = function(region, value, ...) {
#     circos.genomicRect(region, value, col = color_assign(value[[1]]),
#                        border = color_assign(value[[1]]), ...)
#   } )
# 
# 
# ####----添加连线Link
# colors_tmp <- scales::alpha('black', alpha = 1)
# #colline
# circos.genomicLink(HSPs_collin_bed1, HSPs_collin_bed2, 
#                    col = colors_tmp,
#                    border = '#c6dbef',
#                    lwd = 1)
# #tandem
# circos.genomicLink(HSPs_tandem_bed1, HSPs_tandem_bed2, 
#                    col = colors_tmp,
#                    border = '#fec44f',
#                    lwd = 1)
# 
# gene_legend <- Legend(
#   at = c(0, 1, 2, 3), 
#   labels = c('HSP20','HSP60', 'HSP70', 'HSP90'), 
#   labels_gp = gpar(fontsize = 8),
#   title = 'Gene Type', 
#   title_gp = gpar(fontsize = 9), 
#   grid_height = unit(0.4, 'cm'), 
#   grid_width = unit(0.4, 'cm'), 
#   type = 'points', pch = NA, 
#   background = c('#fb8072', '#80b1d3','#fdb462', '#b3de69'))
# 
# Duplicate_legend <- Legend(
#   at = c(1, 2), 
#   labels = c('Collinearity', 'Tandem'),
#   labels_gp = gpar(fontsize = 8), 
#   grid_height = unit(0.5, 'cm'), 
#   grid_width = unit(0.5, 'cm'), 
#   type = c('lines', 'lines'), 
#   pch = NA, 
#   legend_gp = gpar(col = c('#c6dbef', '#fec44f'), lwd = 1),
#   title = 'Gene Duplicate', 
#   title_gp = gpar(fontsize = 9))
# 
# y_coord <- 0.9
# x_coord <- 0.9
# 
# pushViewport(viewport(x = x_coord, y = y_coord))
# grid.draw(gene_legend)
# y_coord <- y_coord
# upViewport()
# 
# pushViewport(viewport(x = x_coord, y = y_coord - 0.75))
# grid.draw(Duplicate_legend)
# y_coord <- y_coord 
# upViewport()
# 
# ####----clear
# circos.clear()
# 
# dev.off()