rm(list = ls())

####----load R Package----####
library(tidyverse)
library(readxl)
library(ggh4x)
library(agricolae)
library(ggpubr)
library(patchwork)
library(ggsignif)
library(showtext)
font_families()
showtext_auto()

####----load Data----####
Expr <- read_xlsx(path = "14.qRT_PCR/表达量-刘_2.xlsx", col_names = T, sheet = 5)

df_group2 <- read_xlsx(path = "6.tree/tree_group_2.xlsx") %>%
  dplyr::filter(Species == "Sorghum_bicolor") %>%
  dplyr::select(label, rename)

Expr1 <- bind_rows(
  Expr %>%
    dplyr::select(1:3) %>%
    tidyr::pivot_longer(cols = -ID, names_to = "Kind", values_to = "Value"),
  Expr %>%
    dplyr::select(1,4:5) %>%
    tidyr::pivot_longer(cols = -ID, names_to = "Kind", values_to = "Value")
)%>%
  dplyr::left_join(df_group2, by = c("ID" = "label")) %>%
  dplyr::select(4,2,3) %>%
  dplyr::mutate(Kind = factor(Kind, levels = c("ck地上", "ln地上",
                                               "ck地下", "ln地下")))


Expr1 %>%
  dplyr::group_by(rename, Kind) %>%
  dplyr::summarise(mean = mean(Value),
                   sd = sd(Value),
                   se = sd / sqrt(n())) %>%
  dplyr::ungroup() %>%
  ggplot() + 
  geom_bar(aes(x = rename, y = mean, fill = Kind), stat = "identity", position = position_dodge(0.85),
           width = 0.75) + 
  geom_errorbar(aes(x = rename, y = mean, ymin = mean - se, ymax = mean + se, group = Kind), position = position_dodge(0.85),
                width = 0.25) + 
  facet_wrap(~rename, scale = "free")

my_comparisons <- list(c("ck地上", "ln地上"), c("ck地下", "ln地下"))


t.test.function <- function(data){
  
  # data = "Expr1"
  df = get(data)
  y.pos = c()
  out_list <- list()
  
  for (id in unique(df$rename)) {
    # id = "SbNPF1.1"
    df_sub <- df %>%
      dplyr::filter(rename == id)
    
    df_sub_stat <- df_sub %>%
      dplyr::group_by(rename, Kind) %>%
      dplyr::summarise(mean = mean(Value), 
                       sd = sd(Value),
                       se = sd/sqrt(n())) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(position = mean + sd)
    
    
    y.pos_sub <- max(df_sub$Value)
    
    y.pos <- c(y.pos, y.pos_sub)
    
    aov_out <- aov(Value ~ Kind, data = df_sub)
    
    aov_out_letters <- LSD.test(aov_out, "Kind", p.adj = "bonferroni")
    
    letters <- aov_out_letters$groups %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = "Group") %>%
      dplyr::mutate(id = id) %>%
      dplyr::mutate(Group = factor(Group, levels = c("ck地上", "ln地上",
                                                   "ck地下", "ln地下"))) %>%
      dplyr::arrange(Group) %>%
      dplyr::mutate(position = df_sub_stat$position * 1.3)
      
    
    out_list[[id]] <- letters
    
  }
  
  out_df <- do.call(rbind, out_list)
  
  
  return(out_df)
}

t.test_out <- t.test.function(data = "Expr1") %>%
  dplyr::mutate(Group = factor(Group, levels = c("ck地上", "ln地上",
                                               "ck地下", "ln地下"))) %>%
  dplyr::rename(rename = id)

Expr1 %>%
  ggplot(aes(x = rename, y = Value)) +
  geom_bar(aes(x = rename, y = Value, fill = Kind), stat = "summary", fun = "mean",  position = position_dodge(0.85),
           width = 0.75) + 
  stat_summary(aes(group = Kind), fun.data = "mean_se", geom="errorbar", width = 0.2, position = position_dodge(0.85)) + 
  # stat_compare_means(comparisons = my_comparisons) 
  # geom_signif(
  #   comparisons = list(c("ck地上", "ln地上"), 
  #                      c("ck地下", "ln地下")),
  #   map_signif_level = TRUE, textsize = 6
  # ) + 
  geom_text(data = t.test_out, aes(x = rename, y = position, label = groups, group = Group),
            position = position_dodge(0.85),
            size = 6
      ) + 
  facet_wrap(~rename, scale = "free") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) + 
  theme_bw() + 
  theme(
    axis.text = element_text(size = 12.5, color = "#000000"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text = element_text(size = 15, color = "#000000"),
    panel.background = element_rect(color = "#000000", linewidth = 1)
  )

ggsave(filename = "14.qRT_PCR/qRT_20250806.pdf",
       height = 12,
       width = 14)
  

  
  

# 
# 
# df_plot <- bind_rows(
#   Expr1 %>%
#     dplyr::group_by(rename) %>%
#     dplyr::summarise(mean = mean(Value1),
#                      sd = sd(Value1),
#                      se = sd / sqrt(n())) %>%
#     dplyr::mutate(Kind = "aboveground"),
#   Expr2 %>%
#     dplyr::group_by(rename) %>%
#     dplyr::summarise(mean = mean(Value2),
#                      sd = sd(Value2),
#                      se = sd / sqrt(n())) %>%
#     dplyr::mutate(Kind = "underground")
# ) %>%
#   bind_cols(
#     bind_rows(letters, letters2) %>%
#       dplyr::select(3)
#   ) 
# 
# df_plot %>%
#   ggplot() + 
#   geom_bar(aes(x = rename, y = mean,  fill = Kind), 
#            stat = "identity",
#            width = 0.7,
#            show.legend = F) + 
#   geom_errorbar(aes(x = rename, ymin = mean-se, ymax = mean+se), width = 0.15) + 
#   geom_text(data = df_plot %>% dplyr::filter(Kind == "aboveground"),
#             aes(x = rename, y = mean + se  + 0.1, label = groups), size = 5) + 
#   geom_text(data = df_plot %>% dplyr::filter(Kind == "underground"),
#             aes(x = rename, y = mean + se  + 0.75, label = groups), size = 5) + 
#   # facet_wrap(~Kind, scale = "free") + 
#   facet_wrap2(~ Kind, strip = strip_themed(
#     background_x = elem_list_rect(fill = c("#31a354", "#d95f0e"))),
#     scale = "free"
#     ) + 
#   labs(x = "ID", y = "Relative expression (fold change)") + 
#   scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + 
#   scale_fill_manual(values = c("#31a354", "#d95f0e")) + 
#   theme_bw() + 
#   theme(
#     axis.text.y = element_text(size = 12.5, color = "#000000"),
#     axis.text.x = element_text(size = 10, color = "#000000", angle = 30, hjust = 0.5, vjust = 0.5),
#     axis.title = element_text(size = 15, color = "#000000"),
#     strip.text = element_text(size = 15, color = "#000000"),
#     strip.background = element_rect(fill = c("#31a354", "#d95f0e")),
#     panel.background = element_rect(linewidth = 1, color = "#000000")
#   )
# 
# ggsave(filename = "14.qRT_PCR/qRT_PCR.pdf",
#        height = 5,
#        width = 10)
# 
# 
