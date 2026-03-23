# 自己编写一个函数
gene_structure_data <- function(gffobject){
  # gffobject = "Malus_domestica_gff_clean"
  gff_tmp <- get(gffobject)
  gene_id <- table(gff_tmp$ID) %>% names()
  
  scale_gff <- lapply(gene_id, function(x){
    # x = "MD00G1000800"
    gff_tmp_sub <- gff_tmp %>%
      dplyr::filter(ID == x) 
    
    if (gff_tmp_sub$strand[1] == "+") {
      start_location <- gff_tmp_sub %>%
        dplyr::pull(start) %>%
        .[1]
      
      gff_tmp_sub_scale <- gff_tmp_sub %>%
        dplyr::mutate(start = start - start_location,
                      end = end - start_location)
      
      
    }else if (gff_tmp_sub$strand[1] == "-") {
      start_location <- gff_tmp_sub %>%
        dplyr::pull(start) %>%
        .[1]
      
      gff_tmp_sub_scale <- gff_tmp_sub %>%
        dplyr::mutate(start = start_location - start,
                      end = start_location - end) %>%
        dplyr::mutate(start = abs(start),
                      end = abs(end))
    }
    
    return(gff_tmp_sub_scale)
    
  }
  ) %>%
    do.call(rbind, .)
  return(scale_gff)
  
}
