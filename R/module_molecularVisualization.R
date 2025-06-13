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
    #print('running the module')
    sdf <- tryCatch(smiles2sdf(maTab())[[1]], error = function(e) print('sth is wrong'))
    if (!is.null(sdf)) {
      plot(sdf)
    } else {
      plot.new()
      text(0.5, 0.5, "Invalid SMILES", cex = 1.5)
    }
  })

}

