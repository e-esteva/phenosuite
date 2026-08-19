require(shiny)
require(stringr)
require(glue)
require(shinyjs)
require(shinyFiles)
require(phenomenalist)
require(cowplot)
require(shinycssloaders)


options(shiny.maxRequestSize=1000*1024^2) 
jsResetCode <- "shinyjs.resetClick = function() {history.go(0)}" # Define the js method that resets the page


ui <- fluidPage(
  tags$head(tags$style(".shiny-notification {position: fixed; 
                                             opacity: 1 ;
                       top: 80% ;
                       left: 40% ;
                       height: 90px;
                       width: 300px}")),
  sidebarLayout( 
    sidebarPanel(
      fileInput("file1", "Choose RDS File",
                accept = c('.rds')
      ),
      checkboxGroupInput('rb0','Select Clustering',choices = "",inline = T,selected = NULL),
      conditionalPanel(condition="!is.null(input.rb0)",checkboxGroupInput('target_cluster','Target Cluster(s)',choices='',inline=T,selected = NULL)),
      
      
      
      tags$hr(),
      
      sliderInput("cluster_res_range", "Clustering Resolution Range:",
                  min = 1, max = 10,
                  value = c(1,3)),
      
     
      
      
      
      
    
    
    ),
    mainPanel(
      
      plotOutput("plots_clusters")%>% withSpinner(color="#0dc5c1"),
      tags$hr(),
      textInput('run_label','Results Name'),
      selectInput('view','full or binary view of sub-clusters',choices=c('full','binary'),selected='full'),
      tags$hr(),
      conditionalPanel(condition="!is.null(input.run_label)",actionButton('run','Initiate Sub-Clustering')),
      tags$hr(),
      # Include shinyjs in the UI
      conditionalPanel(condition='input.run>0',downloadButton(
        outputId = "sub_clustering_download",
        label = "Download Results",
        icon = icon("file-download")
      )),
      tags$hr(),
      useShinyjs(), 
      extendShinyjs(text = jsResetCode, functions = "resetClick"), # Add the js code to the page
      actionButton("reset_button", "Reset Page")
      
     )
  )
)

