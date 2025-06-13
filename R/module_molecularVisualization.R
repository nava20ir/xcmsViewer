library(ChemmineR)
molecularVisualization_UI <- function(id) {
ns <- NS(id) 
shinyjs::useShinyjs()
shiny::tagList(
fluidRow(
  plotOutput(ns("molplot")))
)
}

molecularVisualization <- function(input, output, session, dat, featureSelected=reactive(NULL)) {
  
  maTab <- reactive({
    req( length(featureSelected()$featureid) == 1 ) 
    an <- featureSelected()$annot
    #print('showing colnames an')
    an$ms1$smiles[1]

  })


output$molplot <- renderPlot({
  req(maTab())

  # Convert SMILES to SDF safely
  sdf <- tryCatch({
    smiles <- maTab()
    if (is.null(smiles) || smiles == "") stop("Empty SMILES")
    sdf_list <- smiles2sdf(smiles)
    if (length(sdf_list) == 0 || is.null(sdf_list[[1]])) stop("Invalid conversion")
    sdf_list[[1]]
  }, error = function(e) {
    NULL
  })

  # Plot only if we have a valid SDF
  if (inherits(sdf, "SDF")) {
    plot(sdf)
  } else {
    plot.new()
    text(0.5, 0.5, "Invalid or missing SMILES", cex = 1.5)
  }
})

}

