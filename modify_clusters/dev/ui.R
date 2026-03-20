require(shiny)
require(shinyjs)

options(shiny.maxRequestSize=4000*1024^2) 
jsResetCode <- "shinyjs.resetClick = function() {history.go(0)}"

# Define UI for application that accepts analysis inputs
ui <- fluidPage(
  sidebarLayout( 
    sidebarPanel(
      fileInput("file1", "Choose RDS File",
                accept = c('.rds')
      ),
      checkboxGroupInput('rb0','Select Clustering',choices = "",inline = T,selected = NULL),
      tags$hr(),
      selectInput("use_template", "Use Phenotyping Template?", 
                  choices = c("No", "Yes"), 
                  selected = "No"),
      conditionalPanel(
        condition = "input.use_template == 'Yes'",
        fileInput("template_file", "Upload Phenotyping Template (CSV)",
                  accept = c('.csv', 'text/csv', 'text/comma-separated-values'))
      ),
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
      useShinyjs(),
      extendShinyjs(text = jsResetCode, functions = "resetClick"),
      actionButton("reset_button", "Reset Page")
    )
  )
)
