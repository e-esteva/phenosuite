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
      fileInput("log_odds", "Choose log-odds files",
                accept = c(
                  "text/csv",
                  "text/comma-separated-values,text/plain",
                  ".csv"),multiple = TRUE
      ),
      # New dropdown for celltype operation
      selectInput('celltype_operation', 'Celltype Operation', 
                  choices = c('Union', 'Intersect'), 
                  selected = 'Union'),
      
      # Conditional cluster editing - only shows after celltype_operation is selected
      conditionalPanel(condition="input.celltype_operation != ''",
                       actionButton('edit_clusters',label = 'Edit Cluster Names')),
      conditionalPanel(condition="input.edit_clusters>0 && input.celltype_operation != ''",
                       uiOutput("clusters")),
      	
      checkboxGroupInput('celltype_selection','Celltypes to Analyze',choices = "",inline = T,selected = ""),
      useShinyjs(),# Include shinyjs in the UI
      extendShinyjs(text = jsResetCode, functions = "resetClick"), # Add the js code to the page
      actionButton("reset_button", "Reset Page"),
      
      tags$hr()),
    mainPanel(
      textInput('run_label','Results Name'),
      selectInput('color_scheme','Select Color Scheme',choices = c('Continuous','Divergent'),selected = 'Continuous'),
      selectInput('scale','Scale Sector Width',choices = c(T,F),selected = T),
      numericInput('label_size','Label Size',0.9),
      plotOutput("circos")%>% withSpinner(color="#0dc5c1"),
      
      selectInput('action','Select Operation',choices = c('Integrate','Harmonize')),
      selectInput('view','Select View',choices = c('Global','Unimodal')),
      conditionalPanel(condition="input.view=='Unimodal'",selectInput('ref_selection','Select Reference Celltype',choices = "",selected = "")),
      
      #actionButton("Run", "Render Circos"),
      #plotOutput("plots"),  
      
      
      #conditionalPanel(condition='input.Run>0',downloadButton(
      #  outputId = "circos_download",
      #  label = "Download Results",
      #  icon = icon("file-download")
      #))
      fluidRow(
  	column(6,
    	# Action Button
    	actionButton("Run", "Render Circos")
  	),
  	column(6,
    	conditionalPanel(condition='input.Run>0',downloadButton(
      		outputId = "circos_download",
      		label = "Download Results",
      		icon = icon("file-download")
      	))

  	)
      )
      
    )
  )
  
)
