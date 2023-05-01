require(dplyr)

#' Get functional term data
#'
#' @param species Species name.
#'
#' @return A list of objects for fast functional enrichment.
get_functional_terms <- function(species, all_genes) {
  SPECIES <- list(
    human = list(
      go = "goa_human",
      re = "Homo sapiens",
      kg = "hsa",
      wi = "Homo sapiens"
    ),
    yeast = list(
      go = "sgd",
      re =  "Saccharomyces cerevisiae",
      kg = "sce"
    )
  )

  # Download terms and mappings
  sp <- SPECIES[[species]]
  cat("Loading GO term data\n")
  go <- fenr::fetch_go(species = sp$go)
  cat("Loading Reactome data\n")
  re <- fenr::fetch_reactome(sp$re)
  cat("Loading KEGG data\n")
  kg <- fenr::fetch_kegg(sp$kg)

  # Reduce large GO data
  go$mapping <- go$mapping |>
    dplyr::rename(gene_id = gene_synonym) |>
    dplyr::filter(gene_id %in% all_genes)
  go$terms <- go$terms |>
    dplyr::filter(term_id %in% unique(go$mapping$term_id))

  list(
    go = go,
    re = re,
    kg = kg
  )
}



