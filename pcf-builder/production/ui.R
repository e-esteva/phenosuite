require(shiny)
require(stringr)
require(glue)
require(shinyjs)
require(cowplot)
require(ggplot2)
require(shinyFiles)
require(phenomenalist)
require(shinycssloaders)
# Modify ID

options(shiny.maxRequestSize=1000*1024^2) 
jsResetCode <- "shinyjs.resetClick = function() {history.go(0)}" # Define the js method that resets the page

# Define UI for application that accepts analysis inputs
ui <- fluidPage(
  sidebarLayout( 
    sidebarPanel(
      fileInput("PCFs", "Choose CSV File",
                accept = c(
                  "text/csv",
                  "text/comma-separated-values,text/plain",
                  ".csv"),multiple = TRUE
      ),
      #actionButton('edit_clusters',label = 'Edit Cluster Names'),
      uiOutput("samples") %>% withSpinner(color="#0dc5c1"),
      selectInput('ref_selection','Select Reference Group',choices = "",selected = ""),
      selectInput("celltype_to_analyze", "Available Celltypes", choices = "", selected = ""),
      textInput("run_label", "Results Name"),
      
      tags$hr(),
      
      
    ),
    mainPanel(
      plotOutput("plot") %>% withSpinner(color="#0dc5c1"),
      
      actionButton("confirm_pcf", "Confirm"),
      conditionalPanel(condition='input.confirm_pcf>0',downloadButton(
        outputId = "pcf_download",
        label = "Download Results",
        icon = icon("file-download")
      )),
      
      useShinyjs(),                                           # Include shinyjs in the UI
      extendShinyjs(text = jsResetCode, functions = "resetClick"), # Add the js code to the page
      actionButton("reset_button", "Reset Page")
    )
  )
)
