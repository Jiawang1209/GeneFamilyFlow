get_node_edge <- function(ppi_data) {
  ppi <- get(ppi_data)
  edge_file <- ppi %>% dplyr::select(1, 2, 3) %>% purrr::set_names(c("from", "to", "value"))
  node_file <- c(edge_file %>% dplyr::pull(from), edge_file %>% dplyr::pull(to)) %>%
    table() %>%
    as.data.frame() %>%
    purrr::set_names(c("node", "node_size"))
  list(EdgeObject = edge_file, NodeObject = node_file)
}
