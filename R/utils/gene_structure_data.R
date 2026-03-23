gene_structure_data <- function(gffobject) {
  gff_tmp <- get(gffobject)
  gene_id <- table(gff_tmp$ID) %>% names()
  scale_gff <- lapply(gene_id, function(x) {
    gff_tmp_sub <- gff_tmp %>% dplyr::filter(ID == x)
    if (gff_tmp_sub$strand[1] == "+") {
      start_location <- gff_tmp_sub %>% dplyr::pull(start) %>% .[1]
      gff_tmp_sub_scale <- gff_tmp_sub %>%
        dplyr::mutate(start = start - start_location, end = end - start_location)
    } else {
      start_location <- gff_tmp_sub %>% dplyr::pull(start) %>% .[1]
      gff_tmp_sub_scale <- gff_tmp_sub %>%
        dplyr::mutate(start = abs(start_location - start), end = abs(start_location - end))
    }
    gff_tmp_sub_scale
  }) %>% do.call(rbind, .)
  scale_gff
}
