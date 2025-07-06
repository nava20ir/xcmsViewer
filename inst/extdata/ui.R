library(xcmsViewer)
library(omicsViewer)
#source("landingPage.R")
 source("landingPage.R")
#source("/home/shiny/app/landingPage.R")
ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      #loading-content {
        position: fixed;
        top: 0; left: 0;
        width: 100%; height: 100%;
        background: white;
        z-index: 9999;
        display: flex;
        align-items: center;
        justify-content: center;
      }
      .spinner {
        border: 8px solid #f3f3f3;
        border-top: 8px solid #3498db;
        border-radius: 50%;
        width: 60px;
        height: 60px;
        animation: spin 1s linear infinite;
      }
      @keyframes spin {
        0% { transform: rotate(0deg); }
        100% { transform: rotate(360deg); }
      }
    ")),
    tags$script(HTML("
      $(document).on('shiny:connected', function() {
        $('#loading-content').fadeOut(500);
      });
    "))
  ),
  
  div(id = "loading-content", div(class = "spinner")),		
  uiOutput("uis"),
  uiOutput("aout")
)
