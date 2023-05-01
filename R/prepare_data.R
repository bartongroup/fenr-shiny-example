require(fenr)
require(dplyr)
require(purrr)
require(rlang)

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

  go$mapping <- go$mapping |>
    dplyr::rename(gene_id = gene_synonym)

  terms <- list(
    go = go,
    re = re,
    kg = kg
  )

  # Prepare for fenr
  ontologies <- names(terms)
  purrr::map(ontologies, function(ont) {
    trm <- terms[[ont]]
    fenr::prepare_for_enrichment(trm$terms, trm$mapping, all_genes, feature_name = "gene_id")
  }) |>
    rlang::set_names(ontologies)

}



