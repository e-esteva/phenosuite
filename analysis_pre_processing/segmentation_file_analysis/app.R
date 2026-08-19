library(shiny)
library(DT)
library(ggplot2)
library(dplyr)
library(readr)
library(zip)
library(gridExtra)
library(broom)
source('/srv/shiny-server/phenomenalist/utils/provenance.R')

options(shiny.maxRequestSize = 1*1024^3)
# Define UI
ui <- fluidPage(
  titlePanel("Multi-CSV Data Analysis Tool"),
  
  sidebarLayout(
    sidebarPanel(
      # File upload
      fileInput("files", "Upload CSV files",
                multiple = TRUE,
                accept = c(".csv")),
      
      # Dynamic UI for file labels
      uiOutput("file_labels_ui"),
      
      # Column selection
      uiOutput("column_selection_ui"),
      
      # Action buttons
      br(),
      actionButton("generate_plots", "Generate Plots", 
                   class = "btn-primary"),
      br(), br(),
      
      # Download button
      downloadButton("download_plots", "Download Plots", 
                     class = "btn-success")
    ),
    
    mainPanel(
      # Status and preview
      verbatimTextOutput("status"),
      br(),
      h4("File Preview:"),
      DT::dataTableOutput("file_preview")
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  prov_dir <- file.path(tempdir(), paste0("seg_analysis_", session$token))
  dir.create(prov_dir, showWarnings = FALSE)
  tracker <- ProvenanceTracker$new("segmentation_file_analysis", session, prov_dir)

  # Reactive values to store data
  values <- reactiveValues(
    data_list = NULL,
    all_columns = NULL,
    selected_columns = NULL,
    plots_generated = FALSE,
    plot_files = NULL
  )
  
  # Load and process uploaded files
  observeEvent(input$files, {
    req(input$files)

    tracker$register_input(input$files, input_id = "files")

    # Read all CSV files
    data_list <- list()
    file_paths <- input$files$datapath
    file_names <- input$files$name
    
    for(i in seq_along(file_paths)) {
      tryCatch({
        data_list[[i]] <- read_csv(file_paths[i], show_col_types = FALSE)
        names(data_list)[i] <- file_names[i]
      }, error = function(e) {
        showNotification(paste("Error reading", file_names[i], ":", e$message), 
                         type = "warning")
      })
    }
    
    values$data_list <- data_list
    
    # Get union of all columns (excluding specified ones)
    excluded_cols <- c("XMin", "XMax", "YMin", "YMax", "Analysis.Region", "Analysis Region")
    all_cols <- unique(unlist(lapply(data_list, names)))
    # Also exclude columns containing 'Intensity'
    intensity_cols <- all_cols[grepl("Intensity", all_cols, ignore.case = TRUE)]
    all_excluded <- c(excluded_cols, intensity_cols)
    available_cols <- setdiff(all_cols, all_excluded)
    values$all_columns <- available_cols
    
    showNotification(paste("Loaded", length(data_list), "files successfully"), 
                     type = "default")
  })
  
  # Dynamic UI for file labels
  output$file_labels_ui <- renderUI({
    req(values$data_list)
    
    file_names <- names(values$data_list)
    
    div(
      h4("File Labels:"),
      lapply(seq_along(file_names), function(i) {
        # Get default value from Analysis.Region column if it exists
        default_value <- ""
        data <- values$data_list[[i]]
        
        # Check for Analysis.Region or Analysis Region column
        analysis_col <- NULL
        if("Analysis.Region" %in% names(data)) {
          analysis_col <- "Analysis.Region"
        } else if("Analysis Region" %in% names(data)) {
          analysis_col <- "Analysis Region"
        }
        
        if(!is.null(analysis_col)) {
          unique_regions <- unique(data[[analysis_col]])
          unique_regions <- unique_regions[!is.na(unique_regions)]
          if(length(unique_regions) > 0) {
            default_value <- paste(unique_regions, collapse = ", ")
          }
        }
        
        div(
          strong(file_names[i]),
          textInput(paste0("label_", i), 
                    label = NULL,
                    value = default_value,
                    placeholder = paste("Label for", file_names[i]))
        )
      })
    )
  })
  
  # Dynamic UI for column selection
  output$column_selection_ui <- renderUI({
    req(values$all_columns)
    
    div(
      h4("Select Columns to Analyze:"),
      checkboxGroupInput("selected_columns", 
                         label = NULL,
                         choices = values$all_columns,
                         selected = NULL)
    )
  })
  
  # File preview table
  output$file_preview <- DT::renderDataTable({
    req(values$data_list)
    
    # Create summary of files
    file_summary <- data.frame(
      File = names(values$data_list),
      Rows = sapply(values$data_list, nrow),
      Columns = sapply(values$data_list, ncol),
      stringsAsFactors = FALSE
    )
    
    DT::datatable(file_summary, options = list(pageLength = 10))
  })
  
  # Status output
  output$status <- renderText({
    if(is.null(values$data_list)) {
      return("Please upload CSV files to begin.")
    }
    
    status_text <- paste("Files loaded:", length(values$data_list), "\n")
    status_text <- paste0(status_text, "Available columns:", length(values$all_columns), "\n")
    
    if(!is.null(input$selected_columns)) {
      status_text <- paste0(status_text, "Selected columns:", length(input$selected_columns), "\n")
    }
    
    if(values$plots_generated) {
      status_text <- paste0(status_text, "Plots generated successfully!")
    }
    
    return(status_text)
  })
  
  # Function to determine if column is binary
  is_binary <- function(x) {
    unique_vals <- unique(x[!is.na(x)])
    length(unique_vals) <= 2 && all(unique_vals %in% c(0, 1, TRUE, FALSE, "TRUE", "FALSE", "Yes", "No", "yes", "no"))
  }
  
  # Function to create bar plot for binary data (showing means)
  create_barplot <- function(data_list, column, labels) {
    combined_data <- data.frame()
    
    for(i in seq_along(data_list)) {
      if(column %in% names(data_list[[i]])) {
        # Convert to numeric for mean calculation
        values <- as.numeric(data_list[[i]][[column]])
        values <- values[!is.na(values)]  # Remove NAs
        
        if(length(values) > 0) {
          temp_data <- data.frame(
            Mean_Value = mean(values),
            Group = labels[i],
            stringsAsFactors = FALSE
          )
          combined_data <- rbind(combined_data, temp_data)
        }
      }
    }
    
    if(nrow(combined_data) == 0) return(NULL)
    
    p <- ggplot(combined_data, aes(x = Group, y = Mean_Value, fill = Group)) +
      geom_bar(stat = "identity", alpha = 0.7) +
      geom_text(aes(label = round(Mean_Value, 3)), 
                vjust = -0.5, size = 3) +
      labs(title = paste("Mean Values:", column),
           x = "Group",
           y = paste("Mean", column)) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")
    
    return(p)
  }
  
  # Function to create violin plot for continuous data
  create_violinplot <- function(data_list, column, labels) {
    combined_data <- data.frame()
    
    for(i in seq_along(data_list)) {
      if(column %in% names(data_list[[i]])) {
        temp_data <- data.frame(
          Value = as.numeric(data_list[[i]][[column]]),
          Group = labels[i],
          stringsAsFactors = FALSE
        )
        combined_data <- rbind(combined_data, temp_data)
      }
    }
    
    if(nrow(combined_data) == 0) return(NULL)
    
    # Remove NA values
    combined_data <- combined_data[!is.na(combined_data$Value), ]
    
    if(nrow(combined_data) == 0) return(NULL)
    
    # Perform statistical test (ANOVA if more than 2 groups, t-test if 2 groups)
    stat_result <- NULL
    if(length(unique(combined_data$Group)) == 2) {
      tryCatch({
        test_result <- t.test(Value ~ Group, data = combined_data)
        stat_result <- paste("t-test p-value:", format(test_result$p.value, digits = 4))
      }, error = function(e) {
        stat_result <- "Statistical test failed"
      })
    } else if(length(unique(combined_data$Group)) > 2) {
      tryCatch({
        test_result <- aov(Value ~ Group, data = combined_data)
        p_value <- summary(test_result)[[1]][["Pr(>F)"]][1]
        stat_result <- paste("ANOVA p-value:", format(p_value, digits = 4))
      }, error = function(e) {
        stat_result <- "Statistical test failed"
      })
    }
    
    p <- ggplot(combined_data, aes(x = Group, y = Value, fill = Group)) +
      geom_violin(alpha = 0.7) +
      geom_boxplot(width = 0.1, alpha = 0.8) +
      labs(title = paste("Violin Plot:", column),
           subtitle = stat_result,
           x = "Group",
           y = column) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    return(p)
  }
  
  # Generate plots
  observeEvent(input$generate_plots, {
    req(values$data_list, input$selected_columns)

    tracker$capture_parameters(input)
    tracker$analysis_started()

    # Get labels
    labels <- character(length(values$data_list))
    for(i in seq_along(values$data_list)) {
      label_input <- input[[paste0("label_", i)]]
      labels[i] <- if(is.null(label_input) || label_input == "") {
        names(values$data_list)[i]
      } else {
        label_input
      }
    }
    
    if(any(labels == "")) {
      showNotification("Please provide labels for all files", type = "warning")
      return()
    }
    
    # Create directory for plots inside the provenance dir so the sidecar
    # captures them and the download ZIP ships everything together.
    plot_dir <- file.path(prov_dir, "plots")
    if(dir.exists(plot_dir)) unlink(plot_dir, recursive = TRUE)
    dir.create(plot_dir, recursive = TRUE)
    
    plot_files <- character()
    
    withProgress(message = "Generating plots...", value = 0, {
      
      for(i in seq_along(input$selected_columns)) {
        column <- input$selected_columns[i]
        
        incProgress(1/length(input$selected_columns), 
                    detail = paste("Processing", column))
        
        # Check if column is binary or continuous
        sample_data <- NULL
        for(j in seq_along(values$data_list)) {
          if(column %in% names(values$data_list[[j]])) {
            sample_data <- c(sample_data, values$data_list[[j]][[column]])
          }
        }
        
        if(is.null(sample_data)) next
        
        plot_obj <- NULL
        file_suffix <- ""
        
        if(is_binary(sample_data)) {
          plot_obj <- create_barplot(values$data_list, column, labels)
          file_suffix <- "_barplot.png"
        } else {
          # Try to convert to numeric for continuous data
          numeric_data <- suppressWarnings(as.numeric(sample_data))
          if(!all(is.na(numeric_data))) {
            plot_obj <- create_violinplot(values$data_list, column, labels)
            file_suffix <- "_violinplot.png"
          }
        }
        
        if(!is.null(plot_obj)) {
          file_name <- paste0(gsub("[^A-Za-z0-9]", "_", column), file_suffix)
          file_path <- file.path(plot_dir, file_name)
          
          tryCatch({
            ggsave(file_path, plot_obj, width = 10, height = 6, dpi = 300)
            plot_files <- c(plot_files, file_path)
          }, error = function(e) {
            showNotification(paste("Error saving plot for", column, ":", e$message), 
                             type = "warning")
          })
        }
      }
    })
    
    values$plot_files <- plot_files
    values$plots_generated <- length(plot_files) > 0

    tracker$analysis_completed()

    if(values$plots_generated) {
      showNotification(paste("Generated", length(plot_files), "plots successfully!"),
                       type = "default")
    } else {
      showNotification("No plots could be generated. Please check your data and column selections.",
                       type = "warning")
    }
  })

  # Download handler — bundles plots + provenance sidecar + replay.R
  output$download_plots <- downloadHandler(
    filename = function() {
      paste0("analysis_plots_", Sys.Date(), ".zip")
    },
    content = function(file) {
      req(values$plot_files)

      zip::zip(
        zipfile = file,
        files   = dir(prov_dir, recursive = TRUE),
        root    = prov_dir
      )
    },
    contentType = "application/zip"
  )
  
  # Enable/disable download button
  observe({
    shinyjs::toggleState("download_plots", values$plots_generated)
  })
}

# Run the application
shinyApp(ui = ui, server = server)
