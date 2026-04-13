require(shiny)
require(stringr)
require(glue)
require(shinyjs)
require(shinyFiles)
require(tidyverse)
require(DT)
require(future)
require(promises)
library(stringr)  # for str_detect()
library(tools)    # for file_path_sans_ext()

# Enable parallel processing
plan(multisession)

# analysis engine
options(shiny.maxRequestSize=1000*1024^2) 
jsResetCode <- "shinyjs.resetClick = function() {history.go(0)}" # Define the js method that resets the page
# Define UI for application that accepts analysis inputs
ui <- fluidPage(
  tags$head(tags$style(".shiny-notification {position: fixed; 
                                             opacity: 1 ;
                       top: 35% ;
                       left: 40% ;
                       height: 70px;
                       width: 300px}
                       .sample-config-panel {
                         border: 2px solid #007bff;
                         border-radius: 10px;
                         padding: 20px;
                         margin: 10px 0;
                         background-color: #f8f9fa;
                       }
                       .current-sample {
                         border-color: #28a745 !important;
                         background-color: #e8f5e9 !important;
                       }
                       .completed-sample {
                         border-color: #6c757d !important;
                         background-color: #e9ecef !important;
                       }
                       .global-config-panel {
                         border: 2px solid #17a2b8;
                         border-radius: 10px;
                         padding: 20px;
                         margin: 10px 0;
                         background-color: #e1f7fa;
                       }")),
  
  titlePanel("Multi-Sample Phenomenalist Analysis"),
  
  sidebarLayout( 
    sidebarPanel(
      # Step 1: File Upload
      conditionalPanel(
        condition = "output.current_step == 'upload'",
        h4("Step 1: Upload Files"),
        fileInput("files", "Choose CSV Files",
                  multiple = TRUE,
                  accept = c(
                    "text/csv",
                    "text/comma-separated-values,text/plain",
                    ".csv")
        ),
        tags$hr(),
        textInput('run_label','Results Name Prefix', value = "phenomenalist_run"),
        tags$hr(),
        uiOutput("conditional_button")
      ),
      
      # Step 2: Sample Configuration
      conditionalPanel(
        condition = "output.current_step == 'configure'",
        h4("Step 2: Configure Samples"),
        
        # Configuration mode toggle
        radioButtons("config_mode", "Configuration Mode:",
                     choices = list("Sample-specific configuration" = "individual",
                                    "Apply same configuration to all samples" = "global"),
                     selected = "individual"),
        
        tags$hr(),
        
        # Individual configuration mode
        conditionalPanel(
          condition = "input.config_mode == 'individual'",
          h5(textOutput("current_sample_info")),
          tags$hr(),
          
          # Navigation buttons
          fluidRow(
            column(6,
                   conditionalPanel(
                     condition = "output.show_prev_button",
                     actionButton("prev_sample", "← Previous Sample", class = "btn-secondary")
                   )
            ),
            column(6,
                   conditionalPanel(
                     condition = "output.show_next_button",
                     actionButton("next_sample", "Next Sample →", class = "btn-secondary")
                   )
            )
          ),
          
          tags$hr()
        ),
        
        # Global configuration mode
        conditionalPanel(
          condition = "input.config_mode == 'global'",
          h5("Configure all samples with the same parameters:"),
          tags$hr()
        ),
        
        # Show proceed button
        conditionalPanel(
          condition = "output.show_proceed_button",
          actionButton("proceed_to_analysis", "Proceed to Analysis", class = "btn-success btn-lg")
        )
        
      ),
      
      # Step 3: Analysis
      conditionalPanel(
        condition = "output.current_step == 'analysis'",
        h4("Step 3: Run Analysis"),
        h5("All samples configured. Ready to run analysis."),
        tags$hr(),
        actionButton("run_all", "Run Analysis for All Samples", class = "btn-primary btn-lg"),
        tags$hr(),
        conditionalPanel(
          condition = "input.run_all > 0",
          h5("Processing Status:"),
          verbatimTextOutput("processing_status")
        )
      ),
      
      tags$hr(),
      useShinyjs(),
      extendShinyjs(text = jsResetCode, functions = "resetClick"),
      actionButton("reset_button", "Reset Page", class = "btn-warning")
    ),
    
    mainPanel(
      # Step 1: Welcome/Upload
      conditionalPanel(
        condition = "output.current_step == 'upload'",
        conditionalPanel(
          condition = "!output.files_uploaded",
          h3("Welcome to Multi-Sample Phenomenalist Analysis"),
          p("Please upload one or more CSV files to begin."),
          tags$ul(
            tags$li("Upload multiple segmentation CSV files"),
            tags$li("Configure markers and clustering parameters for each sample individually or globally"),
            tags$li("Run analysis in sequential or parallel mode"),
            tags$li("Download consolidated results")
          )
        ),
        conditionalPanel(
          condition = "output.files_uploaded",
          h3("Files Uploaded Successfully"),
          p("The following files have been uploaded:"),
          tableOutput("uploaded_files_table"),
          p("Click 'Start Sample Configuration' to configure samples.")
        )
      ),
      
      # Step 2: Sample Configuration
      conditionalPanel(
        condition = "output.current_step == 'configure'",
        h3("Sample Configuration"),
        
        # Individual configuration mode
        conditionalPanel(
          condition = "input.config_mode == 'individual'",
          
          # Progress indicator
          div(
            h4("Configuration Progress:"),
            htmlOutput("configuration_progress")
          ),
          
          tags$hr(),
          
          # Current sample configuration
          uiOutput("current_sample_config")
        ),
        
        # Global configuration mode
        conditionalPanel(
          condition = "input.config_mode == 'global'",
          
          # Sample list
          div(
            h4("Samples to be configured:"),
            htmlOutput("samples_list")
          ),
          
          tags$hr(),
          
          # Global configuration panel
          uiOutput("global_config_panel")
        )
      ),
      
      # Step 3: Analysis and Results
      conditionalPanel(
        condition = "output.current_step == 'analysis'",
        h3("Analysis Summary"),
        
        # Configuration summary
        div(
          h4("Configuration Summary:"),
          DTOutput("config_summary_table")
        ),
        
        tags$hr(),
        
        # Results
        conditionalPanel(
          condition = "input.run_all > 0",
          h4("Results"),
          DTOutput("results_table"),
          
          tags$hr(),
          
          # Full error log
          conditionalPanel(
            condition = "output.has_errors",
            h4("Error Details:"),
            verbatimTextOutput("error_log"),
            tags$hr()
          ),
          
          conditionalPanel(
            condition = "output.analysis_complete",
            downloadButton(
              outputId = "download_all_results",
              label = "Download All Results",
              icon = icon("file-download"),
              class = "btn-success btn-lg"
            )
          )
        )
      )
    )
  )
)

# Server logic
server <- function(input, output, session) {
  # Static (non-reactive) copy of the temp dir
  temp_dir_static <- NULL
  
  # Reactive values to store sample data and configurations
  values <- reactiveValues(
    sample_data = list(),
    sample_configs = list(),
    processing_results = list(),
    processing_status = "Ready",
    current_step = "upload",
    current_sample_index = 1,
    temp_dir = NULL,
    global_config = list(), # Store global configuration
    all_markers_union = character(0), # Union of all markers across samples
    all_nuclear_markers_union = character(0), # Union of all nuclear markers
    all_classifier_labels_union = character(0) # Union of all classifier labels
  )
  
  # Helper functions to extract markers
  extract_nuclear_markers <- function(data) {
    cols <- colnames(data)
    # Match columns that contain BOTH nucleus AND a measurement keyword (intensity/mean)
    # This works regardless of word order: "CD3.Nucleus.Mean" or "Nucleus__CD3__Mean"
    has_nucleus <- str_detect(cols, regex('nucleus', ignore_case = TRUE))
    has_measurement <- str_detect(cols, regex('intensity|mean', ignore_case = TRUE))
    nuclear_cols <- cols[has_nucleus & has_measurement]
    message("extract_nuclear_markers found: ", paste(nuclear_cols, collapse = ", "))
    return(nuclear_cols)
  }
  
  extract_all_markers <- function(data) {
    cols <- colnames(data)
    
    # Only consider columns that have BOTH a compartment keyword AND a measurement keyword
    # This prevents metadata columns like "Cell_ID", "Cell_Area", "Nucleus_Count" from sneaking in
    has_compartment <- str_detect(cols, regex('cell|nucleus|cytoplasm|membrane', ignore_case = TRUE))
    has_measurement <- str_detect(cols, regex('intensity|mean', ignore_case = TRUE))
    expression_cols <- cols[has_compartment & has_measurement]
    
    if(length(expression_cols) > 0) {
      # Strip compartment + measurement keywords to extract clean marker names
      markers <- expression_cols
      markers <- gsub("(cell|nucleus|cytoplasm|membrane)", "", markers, ignore.case = TRUE)
      markers <- gsub("(intensity|mean)", "", markers, ignore.case = TRUE)
      # Remove leading/trailing separators and collapse runs of separators
      markers <- gsub("^[\\._ :()]+|[\\._ :()]+$", "", markers)
      markers <- gsub("[\\._ :()]+", "_", markers)
      markers <- unique(markers)
      markers <- markers[markers != "" & markers != "_"]
    } else {
      # Fallback - exclude common non-marker columns (coordinates, morphology, metadata)
      # Exact matches (case-insensitive)
      exclude_exact <- c('x','y','x_min','x_max','y_min','y_max',
                         'xmin','xmax','ymin','ymax',
                         'centroid_x','centroid_y',
                         'centroid_x_um','centroid_y_um',
                         'centroid_x_px','centroid_y_px',
                         'ctrx','ctry',
                         'cell_x_position','cell_y_position',
                         'center_x','center_y',
                         'x_x','y_y',
                         'label','size','cell_id','object_id',
                         'area','perimeter','sample_name',
                         'tissue_category','phenotype','classification',
                         'orig_object','tile_index',
                         'area_convex','min_rot_rect','extent','orientation',
                         'elongation','compactness_circle','compactness_square',
                         'euler_number')
      # Substring matches (case-insensitive) — any column containing these terms
      exclude_substr <- c('centroid','geometry','object','eccentricity',
                          'bbox','convexity','axis','diameter','solidity',
                          'euler')
      cols_lower <- tolower(cols)
      is_exact   <- cols_lower %in% exclude_exact
      is_substr  <- Reduce(`|`, lapply(exclude_substr, function(p) grepl(p, cols_lower, fixed = TRUE)))
      markers <- cols[!(is_exact | is_substr)]
    }
    message("extract_all_markers found: ", paste(markers, collapse = ", "))
    return(markers)
  }
  
  extract_classifier_labels <- function(data) {
    # Case-insensitive search for classifier label column
    cl_col <- colnames(data)[str_detect(colnames(data), regex('^classifier[\\._]?label$', ignore_case = TRUE))]
    if(length(cl_col) > 0) {
      labels <- unique(data[[cl_col[1]]])
      labels <- labels[!is.na(labels)]
      return(labels)
    } else {
      return(character(0))
    }
  }
  
  # Step tracking
  output$current_step <- reactive({
    return(values$current_step)
  })
  outputOptions(output, "current_step", suspendWhenHidden = FALSE)
  
  # Keep your existing reactive
  files_uploaded <- reactive({
    if(!is.null(input$files) && nrow(input$files) > 0){
      message('files uploaded')
      return(TRUE)
    }
    return(FALSE)
  })
  
  
  # Add renderUI
  output$conditional_button <- renderUI({
    if(files_uploaded()) {
      actionButton("start_config", "Start Sample Configuration", class = "btn-primary btn-lg")
    }
  })
  
  # Missing reactive for conditionalPanel elements in UI
  output$files_uploaded <- reactive({
    !is.null(input$files) && nrow(input$files) > 0
  })
  outputOptions(output, "files_uploaded", suspendWhenHidden = FALSE)
  
  # Process uploaded files
  observeEvent(input$files, {
    cat("File upload event triggered\n")
    
    # Clear previous data
    values$sample_data <- list()
    values$sample_configs <- list()
    values$processing_results <- list()
    values$current_sample_index <- 1
    values$global_config <- list()
    values$all_markers_union <- character(0)
    values$all_nuclear_markers_union <- character(0)
    values$all_classifier_labels_union <- character(0)
    
    # Validate that files exist
    if(is.null(input$files) || nrow(input$files) == 0) {
      cat("No files to process\n")
      return()
    }
    
    cat("Processing", nrow(input$files), "files\n")
    
    # Collect all markers across samples for global configuration
    all_markers <- character(0)
    all_nuclear_markers <- character(0)
    all_classifier_labels <- character(0)
    
    # Process each uploaded file
    for(i in 1:nrow(input$files)) {
      file_info <- input$files[i, ]
      cat("Processing file", i, ":", file_info$name, "\n")
      
      # Read the CSV file
      tryCatch({
        # Check if file exists
        if(!file.exists(file_info$datapath)) {
          cat("ERROR: File does not exist at path:", file_info$datapath, "\n")
          showNotification(paste("File not found:", file_info$name), type = "error")
          next
        }
        
        # Try different methods to read CSV
        tbl <- tryCatch({
          read.csv(file_info$datapath, stringsAsFactors = FALSE)
        }, error = function(e) {
          # Try with different separator
          read.csv(file_info$datapath, sep = ";", stringsAsFactors = FALSE)
        })
        
        cat("Successfully read file. Dimensions:", nrow(tbl), "x", ncol(tbl), "\n")
        
        # Validate that we have data
        if(nrow(tbl) == 0) {
          showNotification(paste("Warning: File", file_info$name, "appears to be empty"), 
                           type = "warning")
          next
        }
        
        sample_name <- tools::file_path_sans_ext(file_info$name)
        cat("Sample name:", sample_name, "\n")
        
        # Store sample data
        values$sample_data[[sample_name]] <- list(
          data = tbl,
          file_path = file_info$datapath,
          original_name = file_info$name
        )
        
        # Initialize sample configuration with detected markers
        nuclear_markers <- extract_nuclear_markers(tbl)
        sample_markers <- extract_all_markers(tbl)
        classifier_labels <- extract_classifier_labels(tbl)
        
        # Collect markers for global configuration
        all_markers <- unique(c(all_markers, sample_markers))
        all_nuclear_markers <- unique(c(all_nuclear_markers, nuclear_markers))
        all_classifier_labels <- unique(c(all_classifier_labels, classifier_labels))
        
        cat("Found", length(nuclear_markers), "nuclear markers,", 
            length(sample_markers), "total markers,", 
            length(classifier_labels), "classifier labels\n")
        
        values$sample_configs[[sample_name]] <- list(
          nuclear_markers = nuclear_markers,
          all_markers = sample_markers,
          classifier_labels = classifier_labels,
          selected_nuclear_markers = nuclear_markers, # Pre-select all nuclear markers
          selected_failed_markers = character(0), # No failed markers initially
          selected_classifier_labels = classifier_labels, # Pre-select all classifier labels
          cluster_min = 5,
          cluster_max = 7,
          configured = FALSE
        )
        
        cat("Successfully processed file:", file_info$name, "\n")
        
      }, error = function(e) {
        cat("ERROR processing file", file_info$name, ":", e$message, "\n")
        showNotification(paste("Error reading file", file_info$name, ":", e$message), 
                         type = "error")
      })
    }
    
    # Store union of all markers for global configuration
    values$all_markers_union <- all_markers
    values$all_nuclear_markers_union <- all_nuclear_markers
    values$all_classifier_labels_union <- all_classifier_labels
    
    # Initialize global configuration
    values$global_config <- list(
      selected_nuclear_markers = all_nuclear_markers,
      selected_failed_markers = character(0),
      selected_classifier_labels = all_classifier_labels,
      cluster_min = 5,
      cluster_max = 7
    )
    
    # Check if any files were successfully loaded
    cat("Total samples loaded:", length(values$sample_data), "\n")
    if(length(values$sample_data) == 0) {
      showNotification("No files were successfully loaded. Please check file formats.", 
                       type = "error")
    } else {
      cat("Sample names:", paste(names(values$sample_data), collapse = ", "), "\n")
      showNotification(paste("Successfully loaded", length(values$sample_data), "samples"), 
                       type = "message")
    }
  }, ignoreNULL = TRUE, ignoreInit = TRUE)
  
  # Display uploaded files table
  output$uploaded_files_table <- renderTable({
    req(input$files)
    
    data.frame(
      "File Name" = input$files$name,
      "Size (KB)" = round(input$files$size / 1024, 2),
      "Status" = ifelse(tools::file_path_sans_ext(input$files$name) %in% names(values$sample_data), 
                        "Loaded", "Failed"),
      stringsAsFactors = FALSE
    )
  })
  
  # Start configuration process
  observeEvent(input$start_config, {
    req(values$sample_data)
    if(length(values$sample_data) > 0) {
      values$current_step <- "configure"
      values$current_sample_index <- 1
    } else {
      showNotification("No valid data files found. Please upload valid CSV files.", 
                       type = "error")
    }
  })
  
  # Navigation between samples during configuration (only for individual mode)
  observeEvent(input$next_sample, {
    req(values$sample_data)
    req(input$config_mode == "individual")
    sample_names <- names(values$sample_data)
    
    if(values$current_sample_index < length(values$sample_data)) {
      # Save current sample configuration before moving
      current_sample <- sample_names[values$current_sample_index]
      save_current_sample_config(current_sample)
      
      # Mark current sample as configured
      values$sample_configs[[current_sample]]$configured <- TRUE
      
      values$current_sample_index <- values$current_sample_index + 1
    }
  })
  
  observeEvent(input$prev_sample, {
    req(values$sample_data)
    req(input$config_mode == "individual")
    sample_names <- names(values$sample_data)
    
    if(values$current_sample_index > 1) {
      # Save current sample configuration before moving
      current_sample <- sample_names[values$current_sample_index]
      save_current_sample_config(current_sample)
      
      values$current_sample_index <- values$current_sample_index - 1
    }
  })
  
  # Function to save current sample configuration
  save_current_sample_config <- function(sample_name) {
    req(sample_name)
    
    # Check if sample_name exists in configs
    if(!sample_name %in% names(values$sample_configs)) {
      return()
    }
    
    # Update nuclear markers
    if(!is.null(input$current_nuclear_markers)) {
      values$sample_configs[[sample_name]]$selected_nuclear_markers <- input$current_nuclear_markers
    }
    
    # Update failed markers
    if(!is.null(input$current_failed_markers)) {
      values$sample_configs[[sample_name]]$selected_failed_markers <- input$current_failed_markers
    }
    
    # Update classifier labels
    if(!is.null(input$current_classifier_labels)) {
      values$sample_configs[[sample_name]]$selected_classifier_labels <- input$current_classifier_labels
    }
    
    # Update clustering range
    if(!is.null(input$current_cluster_range)) {
      values$sample_configs[[sample_name]]$cluster_min <- input$current_cluster_range[1]
      values$sample_configs[[sample_name]]$cluster_max <- input$current_cluster_range[2]
    }
  }
  
  # Function to save global configuration
  save_global_config <- function() {
    # Update global nuclear markers
    if(!is.null(input$global_nuclear_markers)) {
      values$global_config$selected_nuclear_markers <- input$global_nuclear_markers
    }
    
    # Update global failed markers
    if(!is.null(input$global_failed_markers)) {
      values$global_config$selected_failed_markers <- input$global_failed_markers
    }
    
    # Update global classifier labels
    if(!is.null(input$global_classifier_labels)) {
      values$global_config$selected_classifier_labels <- input$global_classifier_labels
    }
    
    # Update global clustering range
    if(!is.null(input$global_cluster_range)) {
      values$global_config$cluster_min <- input$global_cluster_range[1]
      values$global_config$cluster_max <- input$global_cluster_range[2]
    }
  }
  
  # Apply global configuration to all samples
  apply_global_config_to_all <- function() {
    save_global_config()
    
    sample_names <- names(values$sample_data)
    for(sample_name in sample_names) {
      if(sample_name %in% names(values$sample_configs)) {
        # Apply global config but only for markers that exist in this sample
        sample_config <- values$sample_configs[[sample_name]]
        
        # Apply nuclear markers (intersection with sample's available markers)
        if(!is.null(values$global_config$selected_nuclear_markers)) {
          available_nuclear <- sample_config$nuclear_markers
          values$sample_configs[[sample_name]]$selected_nuclear_markers <- 
            intersect(values$global_config$selected_nuclear_markers, available_nuclear)
        }
        
        # Apply failed markers (intersection with sample's available markers)  
        if(!is.null(values$global_config$selected_failed_markers)) {
          available_markers <- sample_config$all_markers
          values$sample_configs[[sample_name]]$selected_failed_markers <- 
            intersect(values$global_config$selected_failed_markers, available_markers)
        }
        
        # Apply classifier labels (intersection with sample's available labels)
        if(!is.null(values$global_config$selected_classifier_labels)) {
          available_labels <- sample_config$classifier_labels
          values$sample_configs[[sample_name]]$selected_classifier_labels <- 
            intersect(values$global_config$selected_classifier_labels, available_labels)
        }
        
        # Apply clustering parameters
        values$sample_configs[[sample_name]]$cluster_min <- values$global_config$cluster_min
        values$sample_configs[[sample_name]]$cluster_max <- values$global_config$cluster_max
        
        # Mark as configured
        values$sample_configs[[sample_name]]$configured <- TRUE
      }
    }
  }
  
  # Show/hide navigation buttons (only for individual mode)
  output$show_prev_button <- reactive({
    return(length(values$sample_data) > 0 && values$current_sample_index > 1 && 
             !is.null(input$config_mode) && input$config_mode == "individual")
  })
  outputOptions(output, "show_prev_button", suspendWhenHidden = FALSE)
  
  output$show_next_button <- reactive({
    return(length(values$sample_data) > 0 && values$current_sample_index < length(values$sample_data) &&
             !is.null(input$config_mode) && input$config_mode == "individual")
  })
  outputOptions(output, "show_next_button", suspendWhenHidden = FALSE)
  
  # Show proceed button logic
  output$show_proceed_button <- reactive({
    if(length(values$sample_data) == 0) return(FALSE)
    
    if(is.null(input$config_mode)) return(FALSE)
    
    if(input$config_mode == "global") {
      # For global mode, always show proceed button
      return(TRUE)
    } else {
      # For individual mode, show if we're on the last sample OR if all samples are configured
      on_last_sample <- values$current_sample_index == length(values$sample_data)
      
      # Check if all samples are configured
      all_configured <- TRUE
      if(length(values$sample_configs) > 0) {
        sample_names <- names(values$sample_data)
        all_configured <- all(sapply(sample_names, function(sample_name) {
          if(sample_name %in% names(values$sample_configs)) {
            config <- values$sample_configs[[sample_name]]
            if(!is.null(config) && !is.null(config$configured)) {
              return(isTRUE(config$configured))
            }
          }
          return(FALSE)
        }))
      }
      
      return(on_last_sample || all_configured)
    }
  })
  outputOptions(output, "show_proceed_button", suspendWhenHidden = FALSE)
  
  # Current sample info (only for individual mode)
  output$current_sample_info <- renderText({
    req(values$sample_data)
    req(input$config_mode == "individual")
    sample_names <- names(values$sample_data)
    if(length(sample_names) > 0 && values$current_sample_index <= length(sample_names)) {
      current_sample <- sample_names[values$current_sample_index]
      return(paste("Configuring Sample", values$current_sample_index, "of", length(sample_names), ":", current_sample))
    }
    return("No samples available")
  })
  
  # Configuration progress (only for individual mode)
  output$configuration_progress <- renderUI({
    req(values$sample_data)
    req(values$sample_configs)
    req(input$config_mode == "individual")
    sample_names <- names(values$sample_data)
    
    if(length(sample_names) == 0) {
      return(div("No samples to configure"))
    }
    
    progress_items <- lapply(1:length(sample_names), function(i) {
      sample_name <- sample_names[i]
      is_current <- (i == values$current_sample_index)
      
      # Safely check if sample is configured
      is_configured <- FALSE
      if(sample_name %in% names(values$sample_configs)) {
        config <- values$sample_configs[[sample_name]]
        if(!is.null(config) && !is.null(config$configured)) {
          is_configured <- isTRUE(config$configured)
        }
      }
      
      status_class <- if(is_current) {
        "badge badge-primary"
      } else if(is_configured) {
        "badge badge-success"
      } else {
        "badge badge-secondary"
      }
      
      status_text <- if(is_current) {
        "Current"
      } else if(is_configured) {
        "Configured"
      } else {
        "Pending"
      }
      
      tags$span(
        paste(i, ". ", sample_name, " - "),
        tags$span(status_text, class = status_class),
        style = "margin-right: 15px; display: inline-block; margin-bottom: 5px;"
      )
    })
    
    div(progress_items)
  })
  
  # Samples list (only for global mode)
  output$samples_list <- renderUI({
    req(values$sample_data)
    req(input$config_mode == "global")
    sample_names <- names(values$sample_data)
    
    if(length(sample_names) == 0) {
      return(div("No samples available"))
    }
    
    sample_items <- lapply(1:length(sample_names), function(i) {
      sample_name <- sample_names[i]
      tags$span(
        paste(i, ". ", sample_name),
        style = "margin-right: 15px; display: inline-block; margin-bottom: 5px;"
      )
    })
    
    div(sample_items)
  })
  
  # Current sample configuration UI (for individual mode)
  output$current_sample_config <- renderUI({
    req(values$sample_data)
    req(input$config_mode == "individual")
    sample_names <- names(values$sample_data)
    
    if(length(sample_names) == 0 || values$current_sample_index > length(sample_names)) {
      return(div("No sample data available"))
    }
    
    current_sample <- sample_names[values$current_sample_index]
    
    # Check if config exists
    if(!current_sample %in% names(values$sample_configs)) {
      return(div("Configuration not found for current sample"))
    }
    
    config <- values$sample_configs[[current_sample]]
    
    # Validate config structure
    if(is.null(config)) {
      return(div("Invalid configuration for current sample"))
    }
    
    # Create UI elements with error checking
    tryCatch({
      ui_elements <- list(
        div(
          class = "sample-config-panel current-sample",
          h4(paste("Sample:", current_sample)),
          
          # Nuclear markers (shown when detected)
          if(!is.null(config$nuclear_markers) && length(config$nuclear_markers) > 0) {
            list(
              h5("Nuclear Markers (pre-selected):"),
              checkboxGroupInput(
                inputId = 'current_nuclear_markers',
                label = NULL,
                choices = setNames(config$nuclear_markers, config$nuclear_markers),
                selected = config$selected_nuclear_markers,
                inline = TRUE
              ),
              tags$hr()
            )
          },
          
          # Failed markers
          if(!is.null(config$all_markers) && length(config$all_markers) > 0) {
            list(
              h5("Select Failed Markers:"),
              checkboxGroupInput(
                inputId = 'current_failed_markers',
                label = NULL,
                choices = setNames(config$all_markers, config$all_markers),
                selected = config$selected_failed_markers,
                inline = TRUE
              ),
              tags$hr()
            )
          } else {
            list(
              h5("No markers detected in this sample"),
              tags$hr()
            )
          },
          
          # Classifier labels (shown when detected)
          if(!is.null(config$classifier_labels) && length(config$classifier_labels) > 0) {
            list(
              h5("Classifier Labels (pre-selected):"),
              checkboxGroupInput(
                inputId = 'current_classifier_labels',
                label = NULL,
                choices = setNames(config$classifier_labels, config$classifier_labels),
                selected = config$selected_classifier_labels,
                inline = TRUE
              ),
              tags$hr()
            )
          },
          
          # Clustering parameters
          h5("Clustering Resolution Range:"),
          sliderInput(
            inputId = 'current_cluster_range',
            label = NULL,
            min = 1, max = 10,
            value = c(
              ifelse(is.null(config$cluster_min), 5, config$cluster_min),
              ifelse(is.null(config$cluster_max), 7, config$cluster_max)
            )
          )
        )
      )
      
      return(ui_elements)
    }, error = function(e) {
      return(div(paste("Error creating configuration UI:", e$message)))
    })
  })
  
  # Global configuration panel (for global mode)
  output$global_config_panel <- renderUI({
    req(values$sample_data)
    req(input$config_mode == "global")
    
    if(length(values$sample_data) == 0) {
      return(div("No sample data available"))
    }
    
    # Create UI elements for global configuration
    tryCatch({
      ui_elements <- list(
        div(
          class = "global-config-panel",
          h4("Global Configuration for All Samples"),
          
          # Nuclear markers (shown when detected)
          if(length(values$all_nuclear_markers_union) > 0) {
            list(
              h5("Nuclear Markers (available across all samples):"),
              p("Note: Only markers present in each sample will be applied."),
              checkboxGroupInput(
                inputId = 'global_nuclear_markers',
                label = NULL,
                choices = setNames(values$all_nuclear_markers_union, values$all_nuclear_markers_union),
                selected = values$global_config$selected_nuclear_markers,
                inline = TRUE
              ),
              tags$hr()
            )
          },
          
          # Failed markers
          if(length(values$all_markers_union) > 0) {
            list(
              h5("Select Failed Markers (available across all samples):"),
              p("Note: Only markers present in each sample will be applied."),
              checkboxGroupInput(
                inputId = 'global_failed_markers',
                label = NULL,
                choices = setNames(values$all_markers_union, values$all_markers_union),
                selected = values$global_config$selected_failed_markers,
                inline = TRUE
              ),
              tags$hr()
            )
          } else {
            list(
              h5("No markers detected across samples"),
              tags$hr()
            )
          },
          
          # Classifier labels (shown when detected)
          if(length(values$all_classifier_labels_union) > 0) {
            list(
              h5("Classifier Labels (available across all samples):"),
              p("Note: Only labels present in each sample will be applied."),
              checkboxGroupInput(
                inputId = 'global_classifier_labels',
                label = NULL,
                choices = setNames(values$all_classifier_labels_union, values$all_classifier_labels_union),
                selected = values$global_config$selected_classifier_labels,
                inline = TRUE
              ),
              tags$hr()
            )
          },
          
          # Clustering parameters
          h5("Clustering Resolution Range:"),
          sliderInput(
            inputId = 'global_cluster_range',
            label = NULL,
            min = 1, max = 10,
            value = c(
              ifelse(is.null(values$global_config$cluster_min), 5, values$global_config$cluster_min),
              ifelse(is.null(values$global_config$cluster_max), 7, values$global_config$cluster_max)
            )
          )
        )
      )
      
      return(ui_elements)
    }, error = function(e) {
      return(div(paste("Error creating global configuration UI:", e$message)))
    })
  })
  
  # Proceed to analysis
  observeEvent(input$proceed_to_analysis, {
    req(values$sample_data)
    sample_names <- names(values$sample_data)
    
    if(input$config_mode == "individual") {
      # Save configuration for current sample
      if(length(sample_names) > 0 && values$current_sample_index <= length(sample_names)) {
        current_sample <- sample_names[values$current_sample_index]
        save_current_sample_config(current_sample)
        values$sample_configs[[current_sample]]$configured <- TRUE
      }
    } else if(input$config_mode == "global") {
      # Apply global configuration to all samples
      apply_global_config_to_all()
    }
    
    values$current_step <- "analysis"
    
    # Create temporary directory
    tryCatch({
      temp_path <- file.path(tempdir(), paste0("phenomenalist_", as.integer(Sys.time())))
      values$temp_dir <- temp_path
      temp_dir_static <<- temp_path
      dir.create(values$temp_dir, recursive = TRUE)
    }, error = function(e) {
      showNotification(paste("Error creating temporary directory:", e$message), type = "error")
    })
  })
  
  # Configuration summary table
  output$config_summary_table <- renderDT({
    req(values$sample_configs)
    
    if(length(values$sample_configs) == 0) {
      return(datatable(data.frame(Message = "No configurations available")))
    }
    
    tryCatch({
      summary_data <- do.call(rbind, lapply(names(values$sample_configs), function(sample_name) {
        config <- values$sample_configs[[sample_name]]
        
        data.frame(
          Sample = sample_name,
          `Configuration Mode` = ifelse(exists("config_mode", where = input) && !is.null(input$config_mode), 
                                        ifelse(input$config_mode == "global", "Global", "Individual"), "Individual"),
          `Nuclear Markers` = if(!is.null(config$selected_nuclear_markers) && length(config$selected_nuclear_markers) > 0) 
            paste(config$selected_nuclear_markers, collapse = ", ") else "None",
          `Failed Markers` = if(!is.null(config$selected_failed_markers) && length(config$selected_failed_markers) > 0) 
            paste(config$selected_failed_markers, collapse = ", ") else "None",
          `Classifier Labels` = if(!is.null(config$selected_classifier_labels) && length(config$selected_classifier_labels) > 0) 
            paste(config$selected_classifier_labels, collapse = ", ") else "None",
          `Clustering Range` = paste0(
            ifelse(is.null(config$cluster_min), 5, config$cluster_min), 
            " - ", 
            ifelse(is.null(config$cluster_max), 7, config$cluster_max)
          ),
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
      }))
      
      datatable(summary_data, options = list(pageLength = 10, scrollX = TRUE))
    }, error = function(e) {
      datatable(data.frame(Error = paste("Error creating summary:", e$message)))
    })
  })
  
  # Function to process a single sample
  process_sample <- function(sample_name, sample_data, config, temp_dir, run_label) {
    
    # Use withCallingHandlers to capture the real call stack
    error_info <- NULL
    
    result <- tryCatch(
      withCallingHandlers({
      
      # Validate inputs
      if(is.null(sample_data) || is.null(sample_data$file_path)) {
        stop("Invalid sample data or file path")
      }
      
      if(!file.exists(sample_data$file_path)) {
        stop("Sample data file not found")
      }
      
      # Create sample-specific directory
      sample_dir <- file.path(temp_dir, paste0("sample_", gsub("[^A-Za-z0-9]", "_", sample_name)))
      dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
      
      # Prepare parameters with validation
      nuclear_markers <- if(!is.null(config$selected_nuclear_markers) && length(config$selected_nuclear_markers) > 0) 
        config$selected_nuclear_markers else NULL
      failed_markers <- if(!is.null(config$selected_failed_markers) && length(config$selected_failed_markers) > 0) 
        config$selected_failed_markers else NULL
      classifier_labels <- if(!is.null(config$selected_classifier_labels) && length(config$selected_classifier_labels) > 0) 
        config$selected_classifier_labels else NULL
      
      cluster_min <- ifelse(is.null(config$cluster_min), 5, config$cluster_min)
      cluster_max <- ifelse(is.null(config$cluster_max), 7, config$cluster_max)
      clustering_res <- seq(cluster_min, cluster_max)
      
      sample_label <- paste0(run_label, "_", sample_name)
      
      # ===== DIAGNOSTIC: Log all parameters =====
      message("========== PROCESS_SAMPLE DIAGNOSTICS ==========")
      message("Sample: ", sample_name)
      message("File: ", sample_data$file_path)
      message("Output dir: ", sample_dir)
      message("Label: ", sample_label)
      message("nuclear_markers class: ", class(nuclear_markers), 
              " | length: ", length(nuclear_markers),
              " | is.null: ", is.null(nuclear_markers))
      if(!is.null(nuclear_markers)) message("nuclear_markers: ", paste(nuclear_markers, collapse = " | "))
      message("failed_markers class: ", class(failed_markers), 
              " | length: ", length(failed_markers),
              " | is.null: ", is.null(failed_markers))
      if(!is.null(failed_markers)) message("failed_markers: ", paste(failed_markers, collapse = " | "))
      message("classifier_labels class: ", class(classifier_labels), 
              " | length: ", length(classifier_labels),
              " | is.null: ", is.null(classifier_labels))
      if(!is.null(classifier_labels)) message("classifier_labels: ", paste(classifier_labels, collapse = " | "))
      message("clustering_res: ", paste(clustering_res, collapse = ", "))
      message("=================================================")
      
      # Check if source files exist before sourcing
      source_files <- c(
        '/srv/shiny-server/phenomenalist/analysis_engine_multi_config_streamline/utils/RunPhenomenalist-shiny.R',
        '/srv/shiny-server/phenomenalist/analysis_engine_multi_config_streamline/utils/Phenomenalist-utils-shiny.R'
      )
      
      missing_files <- source_files[!file.exists(source_files)]
      if(length(missing_files) > 0) {
        stop(paste("Required source files not found:", paste(missing_files, collapse = ", ")))
      }
      
      # Source files with individual error catching
      for(f in source_files) {
        message("Sourcing: ", f)
        tryCatch({
          source(f)
          message("  -> Sourced OK")
        }, error = function(e) {
          message("  -> ERROR sourcing ", f, ": ", e$message)
          stop(paste("Failed to source", f, ":", e$message))
        })
      }
      
      # ===== DIAGNOSTIC: Check which functions are available =====
      message("--- Function availability after sourcing ---")
      message("  select_intensity_columns exists: ", exists("select_intensity_columns"))
      message("  RunPhenomenalist.shiny exists: ", exists("RunPhenomenalist.shiny"))
      message("  create_object.mod exists: ", exists("create_object.mod"))
      
      if(!exists("RunPhenomenalist.shiny")) {
        stop("RunPhenomenalist.shiny function not found after sourcing.")
      }
      
      # ===== DIAGNOSTIC: Test select_intensity_columns directly =====
      message("--- Testing select_intensity_columns directly ---")
      tryCatch({
        test_cols <- select_intensity_columns(
          filepath = sample_data$file_path,
          failed_markers = failed_markers,
          nuclear_markers = nuclear_markers,
          prefer_cytoplasm = FALSE
        )
        message("  select_intensity_columns returned ", length(test_cols), " columns")
        if(length(test_cols) > 0) message("  First 5: ", paste(head(test_cols, 5), collapse = " | "))
      }, error = function(e) {
        message("  ERROR in select_intensity_columns: ", e$message)
      })
      
      # ===== Run the analysis =====
      message("--- Calling RunPhenomenalist.shiny ---")
      RunPhenomenalist.shiny(
        segmentation_file = sample_data$file_path,
        label = sample_label,
        failed.markers = failed_markers,
        nuclear.markers = nuclear_markers,
        else_cytoplasm = FALSE,
        clustering_res = clustering_res,
        out_dir = sample_dir,
        classifier_label = classifier_labels,
        min.cells = 10
      )
      
      list(
        sample_name = sample_name,
        status = "Success",
        output_dir = sample_dir,
        message = "Analysis completed successfully"
      )
      
      }, error = function(e) {
        # Capture the REAL call stack before tryCatch unwinds it
        calls <- sys.calls()
        call_text <- paste(lapply(calls, function(c) {
          tryCatch(paste(deparse(c, width.cutoff = 80), collapse = " "), error = function(x) "<unparseable>")
        }), collapse = "\n  -> ")
        error_info <<- paste0("Error: ", e$message, 
                              "\n\nCall stack (most recent last):\n  ", call_text)
        message("========== ERROR CALL STACK ==========")
        message(error_info)
        message("======================================")
      }),
      error = function(e) {
        # Outer tryCatch - shouldn't normally trigger
        error_info <<- paste0("Error: ", e$message)
      }
    )
    
    if(!is.null(error_info)) {
      return(list(
        sample_name = sample_name,
        status = "Error",
        output_dir = NULL,
        message = error_info
      ))
    }
    
    return(result)
  }
  
  # Main processing function
  observeEvent(input$run_all, {
    req(values$sample_data)
    req(values$sample_configs)
    req(input$run_label)
    req(values$temp_dir)
    
    # Validate inputs
    if(length(values$sample_data) == 0) {
      showNotification("No sample data available for processing", type = "error")
      return()
    }
    
    if(is.null(input$run_label) || nchar(trimws(input$run_label)) == 0) {
      showNotification("Please provide a valid run label", type = "error")
      return()
    }
    
    # Clear previous results
    values$processing_results <- list()
    values$processing_status <- "Processing..."
    
    sample_names <- names(values$sample_data)
    
    results <- list()
    
    for(i in seq_along(sample_names)) {
      sample_name <- sample_names[i]
      values$processing_status <- paste("Processing sample", i, "of", length(sample_names), ":", sample_name)
      
      # Force reactive update
      invalidateLater(100, session)
      
      result <- process_sample(
        sample_name = sample_name,
        sample_data = values$sample_data[[sample_name]],
        config = values$sample_configs[[sample_name]],
        temp_dir = values$temp_dir,
        run_label = input$run_label
      )
      
      results[[i]] <- result
      
      # Show intermediate results
      if(result$status == "Error") {
        showNotification(paste("Error processing", sample_name, ":", result$message), type = "error")
      } else {
        showNotification(paste("Successfully processed", sample_name), type = "message")
      }
      
      # Small delay to show progress
      Sys.sleep(0.1)
    }
    
    values$processing_results <- results
    values$processing_status <- "Completed"
    showNotification("All samples processed!", type = "message")
  })
  
  # Display processing status
  output$processing_status <- renderText({
    values$processing_status
  })
  
  # Check if analysis is complete
  output$analysis_complete <- reactive({
    return(length(values$processing_results) > 0 && values$processing_status == "Completed")
  })
  outputOptions(output, "analysis_complete", suspendWhenHidden = FALSE)
  
  # Check if any results have errors
  output$has_errors <- reactive({
    if(length(values$processing_results) == 0) return(FALSE)
    any(sapply(values$processing_results, function(r) identical(r$status, "Error")))
  })
  outputOptions(output, "has_errors", suspendWhenHidden = FALSE)
  
  # Full error log output
  output$error_log <- renderText({
    req(values$processing_results)
    error_results <- Filter(function(r) identical(r$status, "Error"), values$processing_results)
    if(length(error_results) == 0) return("")
    paste(sapply(error_results, function(r) {
      paste0("=== ", r$sample_name, " ===\n", r$message, "\n")
    }), collapse = "\n")
  })
  
  # Display results table
  output$results_table <- renderDT({
    req(values$processing_results)
    
    if(length(values$processing_results) == 0) {
      return(datatable(data.frame(Message = "No results available")))
    }
    
    tryCatch({
      results_df <- do.call(rbind, lapply(values$processing_results, function(result) {
        # Truncate message for table display; full message in error_log
        msg <- ifelse(is.null(result$message), "No message", result$message)
        msg_short <- ifelse(nchar(msg) > 200, paste0(substr(msg, 1, 200), "... [see Error Details below]"), msg)
        data.frame(
          Sample = ifelse(is.null(result$sample_name), "Unknown", result$sample_name),
          Status = ifelse(is.null(result$status), "Unknown", result$status),
          Message = msg_short,
          stringsAsFactors = FALSE
        )
      }))
      
      datatable(results_df, options = list(pageLength = 10, scrollX = TRUE))
    }, error = function(e) {
      datatable(data.frame(Error = paste("Error displaying results:", e$message)))
    })
  })
  
  # Download handler for all results
  output$download_all_results <- downloadHandler(
    filename = function() {
      paste("phenomenalist-multi-sample-", Sys.Date(), ".zip", sep = "")
    },
    content = function(file) {
      req(values$processing_results)
      req(values$temp_dir)
      
      tryCatch({
        # Check if temp directory exists and has content
        if(!dir.exists(values$temp_dir)) {
          stop("Temporary directory not found")
        }
        
        all_files <- list.files(values$temp_dir, recursive = TRUE, full.names = FALSE)
        if(length(all_files) == 0) {
          stop("No output files found")
        }
        
        # Try different zip methods
        zip_success <- FALSE
        
        # Method 1: Try zip package
        if(requireNamespace("zip", quietly = TRUE)) {
          tryCatch({
            zip::zip(
              zipfile = file,
              files = all_files,
              root = values$temp_dir
            )
            zip_success <- TRUE
          }, error = function(e) {
            warning("zip package failed: ", e$message)
          })
        }
        
        # Method 2: Fallback to utils::zip
        if(!zip_success) {
          tryCatch({
            old_wd <- getwd()
            on.exit(setwd(old_wd))
            setwd(values$temp_dir)
            utils::zip(file, all_files)
            zip_success <- TRUE
          }, error = function(e) {
            warning("utils::zip failed: ", e$message)
          })
        }
        
        # Method 3: Try system zip command
        if(!zip_success && Sys.which("zip") != "") {
          tryCatch({
            old_wd <- getwd()
            on.exit(setwd(old_wd))
            setwd(values$temp_dir)
            system(paste("zip -r", shQuote(file), "*"))
            zip_success <- TRUE
          }, error = function(e) {
            warning("system zip failed: ", e$message)
          })
        }
        
        if(!zip_success) {
          stop("All zip methods failed")
        }
        
      }, error = function(e) {
        showNotification(paste("Error creating zip file:", e$message), type = "error")
        # Create a simple text file with error message as fallback
        writeLines(paste("Error creating zip file:", e$message), file)
      })
    },
    contentType = 'application/zip'
  )
  
  # Reset functionality
  observeEvent(input$reset_button, {
    tryCatch({
      js$resetClick()
    }, error = function(e) {
      # Fallback reset method
      session$reload()
    })
  })
  
  # Cleanup on session end
  session$onSessionEnded(function() {
    if(!is.null(temp_dir_static) && dir.exists(temp_dir_static)) {
      tryCatch({
        unlink(temp_dir_static, recursive = TRUE)
      }, error = function(e) {
        # Ignore cleanup errors
      })
    }
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
