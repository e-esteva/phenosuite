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
require('shinycssloaders')



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
      actionButton('edit_clusters',label = 'Edit Cluster Names'),
      conditionalPanel(condition="input.edit_clusters>0",uiOutput("clusters")),
      checkboxGroupInput('celltype_selection','Celltypes to Analyze',choices = "",inline = T,selected = ""),
      useShinyjs(),                                           # Include shinyjs in the UI
      extendShinyjs(text = jsResetCode, functions = "resetClick"), # Add the js code to the page
      actionButton("reset_button", "Reset Page")
      ,
      tags$hr()),
    mainPanel(
      textInput('run_label','Results Name'),
      selectInput('color_scheme','Select Color Scheme',choices = c('Continuous','Divergent'),selected = 'Continuous'),
      selectInput('scale','Scale Sector Width',choices = c(T,F),selected = T),
      numericInput('label_size','Label Size',0.9),
      plotOutput("circos")%>% withSpinner(color="#0dc5c1"),
      
      
      selectInput('self_interactions','Include Self-Interactions',choices = c('No','Yes'),selected = 'No'),
      selectInput('view','Select View',choices = c('Global','Unimodal')),
      #conditionalPanel(condition="input.view=='Unimodal'",checkboxGroupInput('celltype_selection','Select Celltypes to View',choices = "",inline = T,selected = "")),
      conditionalPanel(condition="input.view=='Unimodal'",selectInput('ref_selection','Select Reference Celltype',choices = "",selected = "")),
      actionButton("Run", "Confirm Circos"),
      conditionalPanel(condition='input.Run>0',downloadButton(
        outputId = "circos_download",
        label = "Download Results",
        icon = icon("file-download")
      )),
      plotOutput("plots")
    )
  )
  
)
