require(shiny)
require(stringr)
require(glue)
require(shinyjs)
require(shinyFiles)
require(phenomenalist)
require(openai)
require(cowplot)
require(dplyr)
library(lubridate)

# GPT-powered phenotyping module


api_key <- Sys.getenv("OPENAI_API_KEY")
models_list=read.csv('/srv/shiny-server/Phenoptics-Menu/utils/gpts/gpt_models-updated.csv')
date_=Sys.Date()
if(date_ - max(as.Date(models_list$created)) > 90){ 
	models_list=openai::list_models(openai_api_key=OPENAI_API_KEY)
	models_data=models_list$data %>% subset(owned_by=='system')
	models_data$created=as_datetime(models_data$created)
	models_data = models_data %>% arrange(desc(created))
	models.gpt=str_subset(string=models_data$id,pattern='gpt')
	models.data=models_data[match(models.gpt,models_data$id),]

	write.csv(models.data,glue('/srv/shiny-server/Phenoptics-Menu/utils/gpts/gpt_models-updated.csv'))
}else{
	models.gpt=models_list$id

	
}
rm(OPENAI_API_KEY)

options(shiny.maxRequestSize=1000*1024^2) 
jsResetCode <- "shinyjs.resetClick = function() {history.go(0)}" # Define the js method that resets the page

ui <- fluidPage(
  sidebarLayout( 
    sidebarPanel(
      fileInput("file1", "Choose RDS File",
                accept = c('.rds')
      ),
      checkboxGroupInput('rb0','Select Clustering',choices = "",inline = T,selected = NULL),
      tags$hr(),
      
      #shinyDirButton("directory", "Select output folder", "Please select a folder"),
      #verbatimTextOutput("directorypath"),
      textInput('API_KEY','openAI API key'),
      
      selectInput("prompt_alg", "Prompt Algorithm:",
                  c("Symmetric single choice response" = "v1",
                    "Asymmetric 2-choice response" = "v2")),
      conditionalPanel(condition="input.prompt_alg=='v2'",textInput('tissue','Tissue')),
      selectInput('model','Select GPT model:',models.gpt),
     
      actionButton('run','Initialize Phenotyping'),
      downloadButton(
        outputId = "phenomenalist_download",
        label = "Download Results",
        icon = icon("file-download")
      ),
      useShinyjs(),                                           # Include shinyjs in the UI
      extendShinyjs(text = jsResetCode, functions = "resetClick"), # Add the js code to the page
      actionButton("reset_button", "Reset Page")
    ),
    mainPanel(
      
      plotOutput("plots_ai"),
      plotOutput('plots_ai3'),
      plotOutput('plots_ai2')
      
      
    )
  )
)
