
extractGenesFromGeneSet <- function(gs) {
  gs_list <- list()
  for (geneset in unique(gs$gs_name)) {
    gs_list[[geneset]] <- gs %>% filter(gs$gs_name == geneset) %>%
      pull(gene_symbol)
  }
  return(gs_list)
}
