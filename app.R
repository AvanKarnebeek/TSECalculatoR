library(shiny)

options(shiny.maxRequestSize = 100 * 1024^2)

# Install and load the TSECalculatoR package
if (!requireNamespace("TSECalculatoR", quietly = TRUE)) {
  devtools::install_github("yourusername/TSECalculatoR")
}
library(TSECalculatoR)
library(dplyr)

# UI component
ui <- fluidPage(
  titlePanel("TSECalculatoR"),

  sidebarLayout(
    sidebarPanel(
      h3("About TSECalculatoR"),
      p("For the Bachelor End Thesis Evaluating the T Cell-to-Stroma Enrichment (TSE) Score as a Transcriptomic Predictive Biomarker for Immunotherapy Response, the R package TSECalculatoR was created."),
      p("This package extracts TSE scores, categories, and other biomarkers from a count matrix. It consists of three primary functions:"),
      tags$ul(
        tags$li("calculate_tsescore(countmatrix): Uses GSVA to extract TSE score, global T cell, and stromal cell signatures."),
        tags$li("calculate_tseprofexh(countmatrix): Extracts TSE score, global T cell and stromal cell signatures, proliferation score, and T cell exhaustion score."),
        tags$li("TSE_classify(countmatrix): Classifies patients into TSE categories and provides correlations to centroids with p-values (adapted from Rijnders et al.).")
      ),
      p("The TSECalculatoR can be found and downloaded at",
        tags$a(href = "https://github.com/AvanKarnebeek/TSECalculatoR", "https://github.com/AvanKarnebeek/TSECalculatoR")),

      h4("Upload a Count Matrix"),
      fileInput("file", "Choose CSV File", accept = ".csv"),
      actionButton("calculate", "Calculate TSE Metrics"),
      br(), br(),
      downloadButton("download", "Download Results")
    ),

    mainPanel(
      h3("Results"),
      verbatimTextOutput("summary"),
      tableOutput("table")
    )
  )
)

# Define server
server <- function(input, output, session) {
  # Reactive value to store results
  results <- reactiveVal()

  observeEvent(input$calculate, {
    req(input$file)

    # Load the count matrix
    countmatrix <- read.csv(input$file$datapath, row.names = 1)

    numeric_matrix <- as.matrix(apply(countmatrix, 2, as.numeric))

    rownames(numeric_matrix) <- rownames(countmatrix)

    # Perform calculations using the TSECalculatoR package
    tse_scores <- TSECalculatoR::calculate_tsescore(numeric_matrix)
    tse_prof_exh <- TSECalculatoR::calculate_tseprofexh(numeric_matrix)
    tse_classification <- TSECalculatoR::TSE_classify(numeric_matrix)

    # Combine results into a single data frame
    combined_results <- tse_scores %>%
      left_join(tse_prof_exh %>% dplyr::select(sampleId, Proliferation, ExhTCell), by = "sampleId") %>%
      left_join(tse_classification, by = "sampleId")

    # Store the results
    results(combined_results)
  })

  # Display summary
  output$summary <- renderPrint({
    req(results())
    summary(results())
  })

  # Display table with combined results
  output$table <- renderTable({
    req(results())
    results()
  })

  # Download results
  output$download <- downloadHandler(
    filename = function() {
      paste("TSE_results", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(results(), file, row.names = FALSE)
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)






