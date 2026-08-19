require(shiny)
require(stringr)
require(glue)
require(shinyjs)
require(cowplot)
require(ggplot2)
require(shinyFiles)
require(phenomenalist)

# Modify ID
source('/srv/shiny-server/Phenoptics-Menu/utils/RunPhenomenalist-shiny/RunPhenomenalist-shiny.R')
source('/srv/shiny-server/Phenoptics-Menu/utils/RunPhenomenalist-shiny/phenomenalist-utils-shiny.R')
options(shiny.maxRequestSize=4000*1024^2) 
jsResetCode <- "shinyjs.resetClick = function() {history.go(0)}" # Define the js method that resets the page

# Define UI for application that accepts analysis inputs
ui <- fluidPage(
  sidebarLayout( 
    sidebarPanel(
      fileInput("file1", "Choose RDS File",
                accept = c('.rds')
      ),
      checkboxGroupInput('rb0','Select Clustering',choices = "",inline = T,selected = NULL),
      tags$hr(),
      uiOutput("clusters"),
      
    ),
    mainPanel(
      
      plotOutput("plots"),
      plotOutput('plots_new2'),
      plotOutput('plots_new'),
      actionButton("render", "Update Annotations"),
      downloadButton(
        outputId = "phenomenalist_download",
        label = "Download Results",
        icon = icon("file-download")
      ),
      useShinyjs(),                                           # Include shinyjs in the UI
      extendShinyjs(text = jsResetCode, functions = "resetClick"), # Add the js code to the page
      actionButton("reset_button", "Reset Page")
    )
  )
)
