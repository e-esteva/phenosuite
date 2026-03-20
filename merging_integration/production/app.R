# Spatial Experiment Integration Shiny App
# Load required libraries
library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(shinyjs)
library(DT)
library(plotly)
library(SpatialExperiment)
library(SingleCellExperiment)
library(scater)
library(harmony)
library(uwot)
library(ggplot2)
library(glue)

options(shiny.maxRequestSize=4000*1024^2)

# Define UI
ui <- dashboardPage(
  dashboardHeader(title = "Spatial Experiment Integration Tool"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Upload & Process", tabName = "upload", icon = icon("upload")),
      menuItem("Results", tabName = "results", icon = icon("chart-line"))
    )
  ),
  
  dashboardBody(
    useShinyjs(),  # Initialize shinyjs
    
    # Add custom JavaScript for progress text updates
    tags$script("
      Shiny.addCustomMessageHandler('updateProgressText', function(message) {
        $('#progress_text').text(message);
      });
    "),
    
    tabItems(
      # Upload and Processing Tab
      tabItem(
        tabName = "upload",
        fluidRow(
          box(
            title = "File Upload", status = "primary", solidHeader = TRUE, width = 12,
            fileInput("files", "Upload SpatialExperiment .rds files",
                      multiple = TRUE,
                      accept = c(".rds", ".RDS")),
            br(),
            h4("Processing Options"),
            radioButtons("method", "Choose processing method:",
                         choices = list(
                           "Merge samples" = "merge",
                           "Integration with Harmony batch correction" = "harmony"
                         ),
                         selected = "merge"),
            
            conditionalPanel(
              condition = "input.method == 'harmony'",
              h4("Batch Labels"),
              p("Enter batch labels for each uploaded sample:"),
              uiOutput("batch_inputs")
            ),
            
            br(),
            actionButton("process", "Process Data", 
                         class = "btn-primary btn-lg",
                         disabled = TRUE),
            br(), br(),
            
            # Progress bar and status
            conditionalPanel(
              condition = "input.process > 0",
              div(
                h4("Processing Progress"),
                progressBar(
                  id = "progress_bar",
                  value = 0,
                  total = 100,
                  status = "info",
                  striped = TRUE,
                  
                  title = ""
                ),
                br(),
                div(id = "progress_text", 
                    style = "font-weight: bold; color: #337ab7; text-align: center;",
                    "Ready to process...")
              )
            ),
            br(),
            verbatimTextOutput("status")
          )
        )
      ),
      
      # Results Tab
      tabItem(
        tabName = "results",
        fluidRow(
          box(
            title = "UMAP Visualization", status = "success", solidHeader = TRUE, width = 12,
            plotlyOutput("umap_plot", height = "600px")
          )
        ),
        fluidRow(
          box(
            title = "Sample Summary", status = "info", solidHeader = TRUE, width = 6,
            DT::dataTableOutput("sample_summary")
          ),
          box(
            title = "Download Results", status = "warning", solidHeader = TRUE, width = 6,
            p("Download the processed data and UMAP plot:"),
            br(),
            downloadButton("download", "Download Results (.zip)", 
                           class = "btn-warning btn-lg"),
            br(), br(),
            p("The download includes:"),
            tags$ul(
              tags$li("Processed SpatialExperiment object (.rds)"),
              tags$li("UMAP plot (.png)"),
              tags$li("Sample summary (.csv)")
            )
          )
        )
      )
    )
  )
)

# Define Server
server <- function(input, output, session) {
  # Check if source files exist before sourcing
      source_files <- c(
                        '/srv/shiny-server/Phenoptics-Menu/utils/RunPhenomenalist-shiny/phenomenalist-utils-shiny.R'
      )

      missing_files <- source_files[!file.exists(source_files)]
      if(length(missing_files) > 0) {
        stop(paste("Required source files not found:", paste(missing_files, collapse = ", ")))
      }
      
      # Source files
      for(file in source_files) {
        source(file)
      }
      

  # Reactive values to store data
  values <- reactiveValues(
    uploaded_files = NULL,
    processed_data = NULL,
    umap_plot = NULL,
    sample_info = NULL
  )
  
  # Enable process button when files are uploaded
  observe({
    if (!is.null(input$files)) {
      enable("process")
    } else {
      disable("process")
    }
  })
  
  # Generate batch input fields based on uploaded files
  output$batch_inputs <- renderUI({
    if (is.null(input$files)) return(NULL)
    
    file_names <- input$files$name
    file_names <- gsub("\\.rds$|\\.RDS$", "", file_names)
    
    input_list <- lapply(seq_along(file_names), function(i) {
      textInput(paste0("batch_", i), 
                label = paste("Batch label for", file_names[i], ":"),
                value = paste0("batch_", i),
                placeholder = "e.g., slide1, condition_A, etc.")
    })
    
    do.call(tagList, input_list)
  })
  
  # Process data when button is clicked
  observeEvent(input$process, {
    if (is.null(input$files)) return()
    
    # Initialize progress
    updateProgressBar(session = session, id = "progress_bar", value = 0, title = "Starting...")
    session$sendCustomMessage(type = "updateProgressText", message = "Initializing processing...")
    
    tryCatch({
      # Step 1: Load files (20% progress)
      updateProgressBar(session = session, id = "progress_bar", value = 20, title = "Loading files...")
      output$status <- renderText("Loading RDS files...")
      
      # Load all RDS files
      spe_objects <- list()
      file_names <- gsub("\\.rds$|\\.RDS$", "", input$files$name)
      
      for (i in seq_along(input$files$datapath)) {
        spe_objects[[i]] <- readRDS(input$files$datapath[i])
        colData(spe_objects[[i]])$sample_id <- file_names[i]
        
        # Update progress for each file loaded
        file_progress <- 20 + (i / length(input$files$datapath)) * 10
        updateProgressBar(session = session, id = "progress_bar", value = file_progress, 
                          title = paste("Loaded", i, "of", length(input$files$datapath), "files"))
      }
      
      # Step 2: Harmonize metadata (40% progress)
      updateProgressBar(session = session, id = "progress_bar", value = 40, title = "Harmonizing metadata...")
      output$status <- renderText("Harmonizing colData columns across samples...")
      
      # Harmonize colData before combining
      all_colnames <- unique(unlist(lapply(spe_objects, function(x) colnames(colData(x)))))
      
      for (i in seq_along(spe_objects)) {
        current_cols <- colnames(colData(spe_objects[[i]]))
        missing_cols <- setdiff(all_colnames, current_cols)
        
        if (length(missing_cols) > 0) {
          for (col in missing_cols) {
            colData(spe_objects[[i]])[[col]] <- NA
          }
        }
        
        # Reorder columns
        colData(spe_objects[[i]]) <- colData(spe_objects[[i]])[, all_colnames, drop = FALSE]
      }
      
      if (input$method == "merge") {
        # Step 3: Merge samples (60% progress)
        updateProgressBar(session = session, id = "progress_bar", value = 60, title = "Merging samples...")
        output$status <- renderText("Merging samples...")
        combined_spe <- do.call(cbind, spe_objects)
        
      } else if (input$method == "harmony") {
        # Step 3: Prepare for Harmony integration (50% progress)
        updateProgressBar(session = session, id = "progress_bar", value = 50, title = "Preparing batch correction...")
        output$status <- renderText("Preparing for Harmony integration...")
        
        # Get batch labels from user input
        batch_labels <- character()
        for (i in seq_along(spe_objects)) {
          batch_input <- input[[paste0("batch_", i)]]
          if (is.null(batch_input) || batch_input == "") {
            batch_input <- paste0("batch_", i)
          }
          batch_labels <- c(batch_labels, rep(batch_input, ncol(spe_objects[[i]])))
        }
        
        # Combine objects (colData already harmonized above)
        combined_spe <- do.call(cbind, spe_objects)
        colData(combined_spe)$batch <- batch_labels
        
        # Step 4: Run PCA (65% progress)
        updateProgressBar(session = session, id = "progress_bar", value = 65, title = "Computing PCA...")
        output$status <- renderText("Computing principal components...")
        combined_spe <- runPCA(combined_spe, ncomponents = 50)
        
        # Step 5: Apply Harmony (75% progress)
        updateProgressBar(session = session, id = "progress_bar", value = 75, title = "Running Harmony integration...")
        output$status <- renderText("Applying Harmony batch correction...")
        
        # Use the correct harmony function
        harmony_coords <- harmony::HarmonyMatrix(
          data_mat = reducedDim(combined_spe, "PCA"), 
          meta_data = colData(combined_spe), 
          vars_use = "batch",
          verbose = FALSE
        )
        
        reducedDim(combined_spe, "Harmony") <- harmony_coords
        
        # Step 6: Integration complete (80% progress)
        updateProgressBar(session = session, id = "progress_bar", value = 80, title = "Integration complete!")
        output$status <- renderText("Harmony integration completed successfully!")
      }
      
      # Step 7: Generate UMAP (90% progress)
      updateProgressBar(session = session, id = "progress_bar", value = 90, title = "Generating UMAP...")
      output$status <- renderText("Computing UMAP embedding...")
      
      if (input$method == "harmony" && "Harmony" %in% reducedDimNames(combined_spe)) {
        combined_spe <- runUMAP(combined_spe, dimred = "Harmony", name = "UMAP")
      } else {
        if (!"PCA" %in% reducedDimNames(combined_spe)) {
          combined_spe <- runPCA(combined_spe, ncomponents = 50)
        }
        combined_spe <- runUMAP(combined_spe, dimred = "PCA", name = "UMAP")
      }
      
      updateProgressBar(session = session, id = "progress_bar", value = 85, title = "Initiating clustering...")
      # Add clustering @ low resolution (can sub-cluster later):
      combined_spe = cluster.mod(combined_spe,resolution = 4,max_clust = 300,label='integrated')
      updateProgressBar(session = session, id = "progress_bar", value = 90, title = "Clustering complete!")

      
      # Step 8: Finalize results (100% progress)
      updateProgressBar(session = session, id = "progress_bar", value = 100, title = "Processing complete!")
      
      # Store results
      values$processed_data <- combined_spe
      values$uploaded_files <- input$files
      
      # Create sample summary
      sample_counts <- table(combined_spe$sample_id)
      values$sample_info <- data.frame(
        Sample = names(sample_counts),
        Cell_Count = as.numeric(sample_counts),
        Processing_Method = input$method,
        stringsAsFactors = FALSE
      )
      
      if (input$method == "harmony") {
        batch_counts <- table(colData(combined_spe)$batch)
        values$sample_info$Batch_Label <- sapply(values$sample_info$Sample, function(s) {
          batch_for_sample <- colData(combined_spe)$batch[colData(combined_spe)$sample_id == s][1]
          return(batch_for_sample)
        })
      }
      
      # Final status update
      final_message <- paste("Processing complete!", 
                             ifelse(input$method == "harmony", "Harmony integration", "Sample merging"),
                             "finished successfully. Check the Results tab.")
      output$status <- renderText(final_message)
      
      # Update progress text
      session$sendCustomMessage(type = "updateProgressText", 
                                message = paste("✓ All done!", nrow(combined_spe), "features ×", ncol(combined_spe), "cells"))
      
      # Switch to results tab after a brief delay
      Sys.sleep(1)
      updateTabItems(session, "tabs", "results")
      
    }, error = function(e) {
      updateProgressBar(session = session, id = "progress_bar", value = 0, status = "danger", title = "Error occurred")
      session$sendCustomMessage(type = "updateProgressText", message = "❌ Processing failed")
      output$status <- renderText(paste("Error:", e$message))
    })
  })
  
  # Generate UMAP plot
  output$umap_plot <- renderPlotly({
    if (is.null(values$processed_data)) return(NULL)
    
    spe <- values$processed_data
    umap_coords <- reducedDim(spe, "UMAP")
    
    plot_data <- data.frame(
      UMAP1 = umap_coords[, 1],
      UMAP2 = umap_coords[, 2],
      Sample = spe$sample_id
    )
    
    if ("batch" %in% colnames(colData(spe))) {
      plot_data$Batch <- colData(spe)$batch
    }
    
    p <- ggplot(plot_data, aes(x = UMAP1, y = UMAP2, color = Sample)) +
      geom_point(size = 0.5, alpha = 0.7) +
      theme_minimal() +
      theme(
        legend.position = "right",
        panel.grid = element_blank()
      ) +
      labs(
        title = paste("UMAP -", ifelse(input$method == "harmony", "Harmony Integration", "Merged Samples")),
        x = "UMAP 1",
        y = "UMAP 2"
      )
    
    values$umap_plot <- p
    ggplotly(p, tooltip = c("x", "y", "colour"))
  })
  
  # Sample summary table
  output$sample_summary <- DT::renderDataTable({
    if (is.null(values$sample_info)) return(NULL)
    
    DT::datatable(values$sample_info, 
                  options = list(pageLength = 10, scrollX = TRUE),
                  rownames = FALSE)
  })
  
  # Fixed Download handler
  output$download <- downloadHandler(
    filename = function() {
      paste0("spatial_integration_results_", Sys.Date(), ".zip")
    },
    content = function(file) {
      if (is.null(values$processed_data)) return()
      
      # Create a unique temporary directory
      temp_dir <- file.path(tempdir(), paste0("spatial_results_", Sys.time() %>% 
                                                format("%Y%m%d_%H%M%S")))
      dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)
      
      tryCatch({
        # Save processed data
        rds_file <- file.path(temp_dir, "processed_spatial_experiment.rds")
        saveRDS(values$processed_data, rds_file)
        
        # Save UMAP plot
        if (!is.null(values$umap_plot)) {
          plot_file <- file.path(temp_dir, "umap_plot.png")
          ggsave(plot_file, values$umap_plot, width = 10, height = 8, dpi = 300)
        }
        
        # Save sample summary
        csv_file <- NULL
        if (!is.null(values$sample_info)) {
          csv_file <- file.path(temp_dir, "sample_summary.csv")
          write.csv(values$sample_info, csv_file, row.names = FALSE)
        }
        
        # Get list of files that actually exist
        files_to_zip <- c(rds_file)
        if (!is.null(values$umap_plot) && file.exists(file.path(temp_dir, "umap_plot.png"))) {
          files_to_zip <- c(files_to_zip, file.path(temp_dir, "umap_plot.png"))
        }
        if (!is.null(csv_file) && file.exists(csv_file)) {
          files_to_zip <- c(files_to_zip, csv_file)
        }
        
        # Create zip using base R utils::zip
        old_wd <- getwd()
        setwd(temp_dir)
        
        # Get relative file paths
        rel_files <- basename(files_to_zip)
        
        # Create the zip file
        utils::zip(zipfile = file, files = rel_files)
        
        setwd(old_wd)
        
      }, error = function(e) {
        # Restore working directory on error
        if (exists("old_wd")) setwd(old_wd)
        stop(paste("Download failed:", e$message))
      }, finally = {
        # Clean up temporary directory
        if (dir.exists(temp_dir)) {
          unlink(temp_dir, recursive = TRUE)
        }
      })
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)
