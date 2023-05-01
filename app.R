library(shiny)
library(readr)

source("R/func.R")

de <- readr::read_rds("data/yeast_de.rds")
term_data <- readr::read_rds("data/term_data.rds") |>
  prepare_for_fenr()


ui <- function() {
  shinyUI(fluidPage(
    tags$style("table{font-size: 11px; background-color: #EAF5FF}"),
    titlePanel("fenr example"),
    p("Select a group of genes in the plot to see their functional enrichment."),

    fluidRow(
      column(5,
             radioButtons("plot_type", "Plot type:", choices = c("Volcano" = "volcano", "MA" = "ma"), inline = TRUE),
             plotOutput("main_plot", height = "480px", width = "100%", brush = "plot_brush", hover = "plot_hover")
      ),
      column(7,
             radioButtons("ontology", "Ontology:", choices = c("GO" = "go", "Reactome" = "re", "KEGG" = "kg"), inline = TRUE),
             div(style = 'height: 480px; overflow-y: scroll', tableOutput("enrichment")),
      )
    )
  ))
}

server <- function(input, output, session) {
  # Prevents RStudio from crashing when Shiny window closed manually
  session$onSessionEnded(function() {
    stopApp()
  })

  output$enrichment <- renderTable({
    enrichment_table(de, term_data, input)
  })

  output$main_plot <- renderPlot({
    main_plot(de, input)
  })
}

# Run the application
shinyApp(ui = ui, server = server)

