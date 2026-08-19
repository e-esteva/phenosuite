require(shiny)
require(stringr)
require(glue)
require(shinyjs)
require(shinyFiles)
require(reticulate)
require(circlize)
require(ComplexHeatmap)
require(dplyr)
require(reshape2)
require(shinycssloaders)

# spatial IX - circos

options(shiny.maxRequestSize=330*1024^2) 
jsResetCode <- "shinyjs.resetClick = function() {history.go(0)}" # Define the js method that resets the page


ui <- fluidPage(
  tags$head(tags$style(".shiny-notification {position: fixed; 
                                             opacity: 1 ;
                       top: 80% ;
                       left: 40% ;
                       height: 70px;
                       width: 300px}")),
  sidebarLayout(
    sidebarPanel(
      fileInput("file1", "Choose CSV File",
                accept = c(
                  "text/csv",
                  "text/comma-separated-values,text/plain",
                  ".csv")
      ),
      checkboxGroupInput('celltype_selection','Celltypes to Analyze',choices = "",inline = T,selected = ""),
      useShinyjs(),                                           # Include shinyjs in the UI
      extendShinyjs(text = jsResetCode, functions = "resetClick"), # Add the js code to the page
      actionButton("reset_button", "Reset Page")
      ,
      tags$hr()),
    mainPanel(
      textInput('run_label','Results Name'),
      numericInput('instr_res','Instrument Resolution (microns/pixel)',value = 0.377),
      numericInput('radius_i','Inner Radius (microns)',value = 3),
      numericInput('radius_o','Outer Radius (microns)',value = 30),
      
      selectInput('compute_significance','Compute Effect Size',choices = c('Yes','No')),
      actionButton("Run", "Initiate Analysis"),
      conditionalPanel('input.Run>0',downloadButton(
        outputId = "spatial_dynamics_download",
        label = "Download Results",
        icon = icon("file-download")
      )),
      plotOutput("plots") %>% withSpinner(color="#0dc5c1")
      
      
    )
  )
  
)
