library(xcmsViewer)
library(omicsViewer)
#source("landingPage.R")
 source("landingPage.R")
#source("/home/shiny/app/landingPage.R")
ui <- fluidPage(
  uiOutput("uis"),
  uiOutput("aout")
)