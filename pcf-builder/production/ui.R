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
      selectizeInput("sample_order", "Plot order (drag to reorder, top → bottom)",
                     choices = NULL, selected = NULL, multiple = TRUE,
                     options = list(plugins = list("drag_drop"))),
      textInput("run_label", "Results Name"),
      
      tags$hr(),

      # ---- Heatmap options -------------------------------------------------
      # Only affect the "Mean interaction heatmap" tab. The violin panel is
      # unchanged.
      tags$strong("Heatmap options"),
      checkboxInput("hm_show_values", "Show values in cells", value = TRUE),
      sliderInput("hm_cap", "Colour cap (percentile)",
                  min = 0.80, max = 1.00, value = 0.98, step = 0.01),
      helpText("Normalised PCF is right-skewed; capping keeps one extreme",
               "cell from flattening the rest of the map."),

      tags$hr(),
      
      
    ),
    mainPanel(
      tabsetPanel(
        id = "views",
        tabPanel(
          "Interaction violins",
          plotOutput("plot") %>% withSpinner(color="#0dc5c1"),
          actionButton("confirm_pcf", "Confirm")
        ),
        tabPanel(
          "Mean interaction heatmap",
          plotOutput("heatmap", height = "760px") %>% withSpinner(color="#0dc5c1"),
          actionButton("export_heatmap", "Export Heatmap")
        )
      ),

      tags$hr(),

      # Either output is worth downloading on its own, so the button appears
      # once the user has confirmed a violin *or* exported a heatmap.
      conditionalPanel(condition='input.confirm_pcf>0 || input.export_heatmap>0',
                       downloadButton(
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
