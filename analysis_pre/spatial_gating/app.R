library(shiny)
library(DT)
library(plotly)
library(readr)
source('/srv/shiny-server/phenomenalist/utils/provenance.R')

options(shiny.maxRequestSize = 1*1024^3)
# Define UI
ui <- fluidPage(
  titlePanel("Interactive Spatial Gating Tool"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Choose CSV File",
                accept = c(".csv")),
      
      helpText("Upload a CSV file with columns: XMin, XMax, YMin, YMax"),
      
      br(),
      
      conditionalPanel(
        condition = "output.fileUploaded",
        h4("Instructions:"),
        p("1. Use the lasso select tool in the plot toolbar"),
        p("2. Draw around points you want to select"),
        p("3. Download the filtered data using the button below"),
        br(),
        downloadButton("downloadData", "Download Selected Data", 
                       class = "btn-primary")
      )
    ),
    
    mainPanel(
      conditionalPanel(
        condition = "output.fileUploaded",
        tabsetPanel(
          tabPanel("Interactive Plot", 
                   plotlyOutput("scatterPlot", height = "600px")),
          tabPanel("Original Data", 
                   DT::dataTableOutput("originalTable")),
          tabPanel("Selected Data", 
                   DT::dataTableOutput("selectedTable"))
        )
      ),
      
      conditionalPanel(
        condition = "!output.fileUploaded",
        div(style = "text-align: center; margin-top: 100px;",
            h3("Please upload a CSV file to begin"),
            p("The CSV file should contain columns named: XMin, XMax, YMin, YMax"))
      )
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  prov_dir <- file.path(tempdir(), paste0("spatial_gating_legacy_", session$token))
  dir.create(prov_dir, showWarnings = FALSE)
  tracker <- ProvenanceTracker$new("spatial_gating_legacy", session, prov_dir)

  # Reactive value to store the uploaded data
  data <- reactiveVal(NULL)
  processed_data <- reactiveVal(NULL)
  selected_data <- reactiveVal(NULL)
  
  # Check if file is uploaded
  output$fileUploaded <- reactive({
    return(!is.null(data()))
  })
  outputOptions(output, 'fileUploaded', suspendWhenHidden = FALSE)
  
  # Read and process the uploaded file
  observeEvent(input$file, {
    req(input$file)
    
    tryCatch({
      # Read the CSV file
      df <- read_csv(input$file$datapath)
      
      # Check if required columns exist
      required_cols <- c("XMin", "XMax", "YMin", "YMax")
      if (!all(required_cols %in% names(df))) {
        showNotification(
          paste("Error: CSV must contain columns:", paste(required_cols, collapse = ", ")),
          type = "error"
        )
        return()
      }
      
      # Calculate x and y coordinates (averages)
      df$x <- (df$XMin + df$XMax) / 2
      df$y <- (df$YMin + df$YMax) / 2
      
      # Store the data
      data(df)
      processed_data(df)
      
      showNotification("File uploaded and processed successfully!", type = "message")
      
    }, error = function(e) {
      showNotification(paste("Error reading file:", e$message), type = "error")
    })
  })
  
  # Create the interactive scatter plot
  output$scatterPlot <- renderPlotly({
    req(processed_data())
    
    df <- processed_data()
    
    p <- plot_ly(
      data = df,
      x = ~x, 
      y = ~y,
      type = 'scatter',
      mode = 'markers',
      marker = list(size = 8, opacity = 0.7),
      text = ~paste("Row:", 1:nrow(df), "<br>X:", round(x, 3), "<br>Y:", round(y, 3)),
      hovertemplate = "%{text}<extra></extra>",
      source = "scatterplot"
    ) %>%
      layout(
        title = "Interactive Scatter Plot - Use Lasso Select Tool",
        xaxis = list(title = "X (Average of XMin and XMax)"),
        yaxis = list(title = "Y (Average of YMin and YMax)"),
        dragmode = "lasso"
      ) %>%
      config(
        modeBarButtonsToAdd = list("lasso2d"),
        displaylogo = FALSE,
        modeBarButtonsToRemove = c("autoScale2d", "resetScale2d", "toggleHover", 
                                   "toggleSpikelines", "zoom2d", "pan2d", 
                                   "select2d", "zoomIn2d", "zoomOut2d")
      )
    
    p
  })
  
  # Handle lasso selection
  observeEvent(event_data("plotly_selected", source = "scatterplot"), {
    req(processed_data())
    
    selection <- event_data("plotly_selected", source = "scatterplot")
    
    if (!is.null(selection) && nrow(selection) > 0) {
      # Get the indices of selected points
      selected_indices <- selection$pointNumber + 1  # R is 1-indexed
      
      # Filter the original data based on selected points
      df <- processed_data()
      selected_df <- df[selected_indices, ]
      
      selected_data(selected_df)
      
      showNotification(paste("Selected", nrow(selected_df), "data points"), type = "message")
    } else {
      selected_data(NULL)
    }
  })
  
  # Display original data table
  output$originalTable <- DT::renderDataTable({
    req(data())
    DT::datatable(
      data(),
      options = list(
        scrollX = TRUE,
        pageLength = 10,
        dom = 'Bfrtip'
      ),
      caption = paste("Original Data (", nrow(data()), " rows)")
    )
  })
  
  # Display selected data table
  output$selectedTable <- DT::renderDataTable({
    if (is.null(selected_data())) {
      return(DT::datatable(
        data.frame(Message = "No data selected. Use the lasso tool on the plot to select points."),
        options = list(dom = 't'),
        rownames = FALSE,
        colnames = ""
      ))
    }
    
    DT::datatable(
      selected_data(),
      options = list(
        scrollX = TRUE,
        pageLength = 10,
        dom = 'Bfrtip'
      ),
      caption = paste("Selected Data (", nrow(selected_data()), " rows)")
    )
  })
  
  # Download handler for selected data
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("selected_data_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      req(selected_data())
      
      # Remove the calculated x and y columns before downloading
      download_data <- selected_data()
      download_data$x <- NULL
      download_data$y <- NULL
      
      write_csv(download_data, file)
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)
