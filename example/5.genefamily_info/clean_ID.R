args <- commandArgs(trailingOnly = TRUE)
library(tidyverse)
read_delim(file = args[1], col_names = F, delim = "\t") %>%
        set_names(c('label')) %>%
        dplyr::mutate(label = str_remove(string = label, pattern = '\\.\\d+\\.*\\w*$')) %>%
        dplyr::mutate(label = str_remove(string = label, pattern = '-\\w+$')) %>%
        dplyr::mutate(label = str_remove(string = label, pattern = '_P\\d{3}$')) %>%
        write.table(file = paste0(args[1],'.clean.out'), quote=F, row.names=F, col.names=F)
