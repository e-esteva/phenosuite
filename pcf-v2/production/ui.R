require(shiny)
require(shinyFiles)
require(stringr)
require(glue)
require(shinyjs)
require(cowplot)
require(reticulate)
require(reshape2)
require(ggpubr)
require(ggsignif)
require('shinycssloaders')
require(dplyr)

options(shiny.maxRequestSize = 100 * 1024^2)
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
      fileInput("vectras", "Choose vectra annotation file(s)",
                accept = c(
                  "text/csv",
                  "text/comma-separated-values,text/plain",
                  ".csv"),multiple = TRUE
      ),
      checkboxGroupInput("celltypes_", "Celltypes", choices = "", inline = T, selected = ""),
      #checkboxGroupInput("cell1", "Select Central Celltype", choices = "", inline = T, selected = NULL),
      selectInput('ref_selection','Select Reference Celltype',choices = "",selected = ""),
      textInput("run_label", "Results Name"),
      numericInput("instr_res", "Instrument Resolution (microns/pixel)", value = 0.377),
      numericInput("radius", "Radius (microns)", value = 30),
      actionButton("Run", "Initiate Analysis"),
      tags$hr(),
      conditionalPanel('input.Run>0',downloadButton(
        outputId = "pcf_download",
        label = "Download Results",
        icon = icon("file-download")
      ))
      ,
      tags$hr(),
      actionButton("reset_button", "Reset Page"),
      useShinyjs(), # Include shinyjs in the UI
      extendShinyjs(text = jsResetCode, functions = "resetClick"), # Add the js code to the page
      
      
      tags$hr()
    ),
    mainPanel(
      plotOutput("plots") %>% withSpinner(color="#0dc5c1")
    )
  )
)
