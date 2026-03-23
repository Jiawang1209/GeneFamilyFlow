get_node_edge <- function(ppi_data){
  ppi <- get(ppi_data)
  
  edge_file <- ppi %>% dplyr::select(1,2,3) %>% purrr::set_names(c("from", "to", "value"))
  
  node_file <- c(edge_file %>% pull(from),
               edge_file %>% pull(to)) %>%
    table() %>%
    as.data.frame() %>%
    purrr::set_names(c("node", "node_size"))
  
  return(list(EdgeObject = edge_file,
              NodeObject = node_file)
         )
}
