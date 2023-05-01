require(ggplot2)
require(dplyr)
require(rlang)
require(fenr)

#' Volcano plot
#'
#' @param d Tibble with x, y and FDR
#' @param fdr_limit FDR limit below which point will be plotted black
#'
#' @return A ggplot object
plot_volcano <- function(d, fdr_limit = 0.05) {
  x <- y <- FDR <- NULL
  sres <- d |>
    dplyr::filter(FDR <= fdr_limit)
  g <- ggplot(d, aes(x, y)) +
    theme_classic() +
    theme(
      panel.grid = element_blank(),
      text = element_text(size = 18)
    ) +
    geom_point(size = 0.2) +
    geom_vline(xintercept = 0, colour = "grey70") +
    labs(x = expression(log[2]~FC), y = expression(-log[10]~P)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03)))
  if (nrow(sres) > 0) {
    g <- g + geom_hline(yintercept = -log10(max(sres$PValue)), colour = "red", linetype = "dashed", alpha = 0.2)
  }
  g
}

#' MA plot
#'
#' @param d Tibble with x, y and FDR
#' @param fdr_limit FDR limit below which point will be plotted black
#'
#' @return A ggplot object
plot_ma <- function(d, fdr_limit = 0.05) {
  x <- y <- NULL
  ggplot(d, aes(x, y)) +
    theme_classic() +
    theme(
      panel.grid = element_blank(),
      text = element_text(size = 18)
    ) +
    geom_point(data = d[d$FDR > fdr_limit,], size = 0.1, colour = "grey50") +
    geom_point(data = d[d$FDR <= fdr_limit,], size = 0.2, colour = "black") +
    geom_hline(yintercept = 0, colour = "grey70") +
    labs(x = expression(log[2]~CPM), y = expression(log[2]~FC))
}



#' Helper function for Shiny app
#'
#' @param de Differential expression results
#' @param input App input
#'
#' @return A tibble with x, y
get_xy_data <- function(de, input) {
  logFC <- PValue <- logCPM <- NULL
  if (input$plot_type == "volcano") {
    xy_data <- de |>
      dplyr::mutate(x = logFC, y = -log10(PValue))
  } else if (input$plot_type == "ma") {
    xy_data <- de  |>
      dplyr::mutate(x = logCPM, y = logFC)
  }
  xy_data
}

#' Helper function for Shiny app
#'
#' @param de Differential expression results
#' @param term_data All functional term data
#' @param input App input
#' @param max_points Maximum number of points allowed to select
#'
#' @return A tibble with functional enrichment results
enrichment_table <- function(de, term_data, input, max_points = 10000) {
  p_adjust <- n_with_sel <- NULL

  xy_data <- get_xy_data(de, input)
  terms <- term_data[[input$ontology]]
  gene2symbol <- rlang::set_names(de$gene_symbol, de$gene_id)
  sel <- NULL
  fe <- NULL
  if (!is.null(input$plot_brush)) {
    brushed <- brushedPoints(xy_data, input$plot_brush)
    sel <- brushed$gene_id
    n <- length(sel)
    if (n > 0 && n <= max_points) {
      fe <- fenr::functional_enrichment(de$gene_id, sel, terms, feat2name = gene2symbol)
      if(!is.null(fe))
        fe <- dplyr::filter(fe, p_adjust < 0.05 & n_with_sel > 2)
    } else if (n > 0) {
      fe <- data.frame(Error = paste('only', max_points, 'points can be selected.'))
    }
  }
  fe
}

#' Helper function for Shiny app
#'
#' @param de Differential expression results
#' @param input App input
#'
#' @return Data for volcano or ma plot
main_plot <- function(de, input) {
  xy_data <- get_xy_data(de, input)
  if (input$plot_type == "volcano") {
    plot_volcano(xy_data)
  } else if (input$plot_type == "ma") {
    plot_ma(xy_data)
  }
}


#' Convert ontology data into fenr format
#'
#' @param terms A list of lists, each containing mapping and terms
#' @param all_genes A character vector with gene background
#'
#' @return A list of fenr objects
prepare_for_fenr <- function(terms, all_genes) {
  ontologies <- names(terms)
  purrr::map(ontologies, function(ont) {
    trm <- terms[[ont]]
    fenr::prepare_for_enrichment(trm$terms, trm$mapping, all_genes, feature_name = "gene_id")
  }) |>
    rlang::set_names(ontologies)
}

