library(ChemmineR)
library(ChemmineOB)

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

  req(maTab())

  print("Starting molplot rendering...")

  sdf <- tryCatch({
    smi <- maTab()
    print(paste("SMILES input:", smi))

    if (is.null(smi) || !nzchar(smi)) stop("Empty SMILES")

    sdf <-  ChemmineR::smiles2sdf(smi)[[1]]
    #sdf <- gen2D(sdf)
    #coords <- atomblock(sdf)

    #print("Atomblock coordinates:")
    #print(coords)

    # Defensive check for finite coordinates
    #if (any(!is.finite(as.numeric(coords[, 1:2])))) {
    #  stop("Non-finite coordinates found")
    #}

    sdf
  }, error = function(e) {
    print(paste("Caught error in sdf processing:", e$message))
    NULL
  })

  if (inherits(sdf, "SDF")) {
    tryCatch({
      print("Plotting molecule...")
      plot(sdf)
    }, error = function(e) {
      print(paste("Plotting error:", e$message))
      plot.new()
      text(0.5, 0.5, "Plot error", cex = 1.5)
    })
  } else {
    plot.new()
    text(0.5, 0.5, "Invalid SMILES", cex = 1.5)
  }
})

}

