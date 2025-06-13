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

  sdf <- tryCatch({
    s <- maTab()
    if (is.null(s) || !nzchar(s)) stop("Empty SMILES")
    sdf <- smiles2sdf(s)[[1]]
    sdf <- gen2D(sdf)

    coords <- atomblock(sdf)
    # Defensive check: ensure X and Y are finite
    if (any(!is.finite(coords[, 1:2]))) stop("Non-finite coordinates")

    sdf
  }, error = function(e) {
    message("Error: ", e$message)
    NULL
  })

  if (inherits(sdf, "SDF")) {
    plot(sdf)
  } else {
    plot.new()
    text(0.5, 0.5, "Cannot render molecule", cex = 1.5)
  }
})

}

