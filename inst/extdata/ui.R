library(xcmsViewer)
library(omicsViewer)
source("landingPage.R")
# source("landingPage.R")

ui <- fluidPage(
  uiOutput("uis"),
  uiOutput("aout")
)