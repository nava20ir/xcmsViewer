library(shiny)
library(ChemmineR)

ui <- fluidRow(
  
  titlePanel("SMILES Viewer"),
  textInput("smiles", "Enter SMILES:", "CCO"),
  plotOutput("molplot")
)

server <- function(input, output) {
  output$molplot <- renderPlot({
    req(input$smiles)
    
    sdf <- tryCatch(smiles2sdf(input$smiles)[[1]], error = function(e) NULL)
    if (!is.null(sdf)) {
      plot(sdf)
    } else {
      plot.new()
      text(0.5, 0.5, "Invalid SMILES", cex = 1.5)
    }
  })
}

shinyApp(ui, server)