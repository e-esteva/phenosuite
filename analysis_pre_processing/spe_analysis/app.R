# Set maximum upload size to 1GB (adjust as needed)
options(shiny.maxRequestSize = 1*1024^3)

library(shiny)
library(DT)
library(ggplot2)
library(dplyr)
library(zip)
library(gridExtra)
library(broom)
library(SpatialExperiment)
library(SingleCellExperiment)
library(scales)  # For percentage formatting
source('/srv/shiny-server/phenomenalist/utils/provenance.R')

# Define UI
ui <- fluidPage(
  titlePanel("Spatial Experiment Data Analysis Tool"),
  
  sidebarLayout(
    sidebarPanel(
      # File upload
      fileInput("files", "Upload .rds files (Spatial Experiment objects)",
                multiple = TRUE,
                accept = c(".rds")),
      
      # Dynamic UI for file labels
      uiOutput("file_labels_ui"),
      
      # Column selection
      uiOutput("column_selection_ui"),
      
      # Categorical analysis options (shown when categorical columns are selected)
      uiOutput("categorical_options_ui"),
      
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
  prov_dir <- file.path(tempdir(), paste0("spe_analysis_", session$token))
  dir.create(prov_dir, showWarnings = FALSE)
  tracker <- ProvenanceTracker$new("spe_analysis", session, prov_dir)

  # Reactive values to store data
  values <- reactiveValues(
    data_list = NULL,
    all_columns = NULL,
    selected_columns = NULL,
    categorical_columns = NULL,
    plots_generated = FALSE,
    plot_files = NULL
  )
  
  # Load and process uploaded files
  observeEvent(input$files, {
    req(input$files)
    
    # Read all .rds files
    data_list <- list()
    file_paths <- input$files$datapath
    file_names <- input$files$name
    
    for(i in seq_along(file_paths)) {
      tryCatch({
        spe_obj <- readRDS(file_paths[i])
        
        # Check if it's a SpatialExperiment or SingleCellExperiment object
        if(!(is(spe_obj, "SpatialExperiment") || is(spe_obj, "SingleCellExperiment"))) {
          showNotification(paste("Warning:", file_names[i], "is not a Spatial/SingleCell Experiment object"), 
                           type = "warning")
          next
        }
        
        data_list[[i]] <- spe_obj
        names(data_list)[i] <- file_names[i]
      }, error = function(e) {
        showNotification(paste("Error reading", file_names[i], ":", e$message), 
                         type = "warning")
      })
    }
    
    values$data_list <- data_list
    
    # Get union of all colData columns (excluding specified ones)
    excluded_cols <- c("XMin", "XMax", "YMin", "YMax", "Analysis.Region", "Analysis Region")
    all_cols <- unique(unlist(lapply(data_list, function(x) {
      if(is.null(colData(x))) return(character(0))
      return(colnames(colData(x)))
    })))
    
    # Also exclude columns containing 'Intensity'
    intensity_cols <- all_cols[grepl("Intensity", all_cols, ignore.case = TRUE)]
    all_excluded <- c(excluded_cols, intensity_cols)
    available_cols <- setdiff(all_cols, all_excluded)
    values$all_columns <- available_cols
    
    showNotification(paste("Loaded", length(data_list), "spatial experiment objects successfully"), 
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
        spe_obj <- values$data_list[[i]]
        
        # Check for Analysis.Region or Analysis Region column in colData
        analysis_col <- NULL
        col_data <- colData(spe_obj)
        
        if("Analysis.Region" %in% colnames(col_data)) {
          analysis_col <- "Analysis.Region"
        } else if("Analysis Region" %in% colnames(col_data)) {
          analysis_col <- "Analysis Region"
        }
        
        if(!is.null(analysis_col)) {
          unique_regions <- unique(col_data[[analysis_col]])
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
      h4("Select Metadata Columns to Analyze:"),
      checkboxGroupInput("selected_columns", 
                         label = NULL,
                         choices = values$all_columns,
                         selected = NULL)
    )
  })
  
  # Dynamic UI for categorical analysis options
  output$categorical_options_ui <- renderUI({
    req(input$selected_columns)
    
    # Determine which selected columns are categorical
    categorical_cols <- character()
    
    for(col in input$selected_columns) {
      # Check if column is categorical across all datasets
      is_categorical <- TRUE
      for(spe_obj in values$data_list) {
        col_data <- colData(spe_obj)
        if(col %in% colnames(col_data)) {
          col_values <- col_data[[col]]
          # Check if it's categorical (character/factor or limited unique values that are strings)
          if(is.numeric(col_values)) {
            unique_vals <- unique(col_values[!is.na(col_values)])
            if(length(unique_vals) > 10) {  # If more than 10 unique numeric values, treat as continuous
              is_categorical <- FALSE
              break
            }
          }
        }
      }
      if(is_categorical) {
        categorical_cols <- c(categorical_cols, col)
      }
    }
    
    values$categorical_columns <- categorical_cols
    
    if(length(categorical_cols) > 0) {
      div(
        h4("Categorical Analysis Options:"),
        selectInput("categorical_column", 
                    "Select categorical column for proportion analysis:",
                    choices = categorical_cols,
                    selected = categorical_cols[1]),
        
        conditionalPanel(
          condition = "input.categorical_column != null && input.categorical_column != ''",
          uiOutput("celltype_selection_ui"),
          
          radioButtons("comparison_type", 
                       "Comparison type:",
                       choices = list(
                         "Compare one cell type across all samples" = "single_celltype",
                         "Compare all cell types within each sample" = "all_celltypes"
                       ),
                       selected = "single_celltype")
        )
      )
    }
  })
  
  # Dynamic UI for cell type selection
  output$celltype_selection_ui <- renderUI({
    req(input$categorical_column, input$comparison_type)
    
    if(input$comparison_type == "single_celltype") {
      # Get all unique cell types from the selected categorical column
      all_celltypes <- character()
      for(spe_obj in values$data_list) {
        col_data <- colData(spe_obj)
        if(input$categorical_column %in% colnames(col_data)) {
          celltypes <- unique(col_data[[input$categorical_column]])
          celltypes <- celltypes[!is.na(celltypes)]
          all_celltypes <- c(all_celltypes, as.character(celltypes))
        }
      }
      all_celltypes <- unique(all_celltypes)
      
      selectInput("selected_celltype",
                  "Select cell type to compare:",
                  choices = all_celltypes,
                  selected = all_celltypes[1])
    }
  })
  
  # File preview table
  output$file_preview <- DT::renderDataTable({
    req(values$data_list)
    
    # Create summary of spatial experiment objects
    file_summary <- data.frame(
      File = names(values$data_list),
      Cells = sapply(values$data_list, ncol),
      Features = sapply(values$data_list, nrow),
      Metadata_Columns = sapply(values$data_list, function(x) ncol(colData(x))),
      Object_Class = sapply(values$data_list, function(x) class(x)[1]),
      stringsAsFactors = FALSE
    )
    
    DT::datatable(file_summary, options = list(pageLength = 10))
  })
  
  # Status output
  output$status <- renderText({
    if(is.null(values$data_list)) {
      return("Please upload .rds files containing Spatial Experiment objects to begin.")
    }
    
    status_text <- paste("Objects loaded:", length(values$data_list), "\n")
    status_text <- paste0(status_text, "Available metadata columns:", length(values$all_columns), "\n")
    
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
  
  # Function to determine if column is categorical
  is_categorical <- function(x) {
    if(is.character(x) || is.factor(x)) return(TRUE)
    if(is.numeric(x)) {
      unique_vals <- unique(x[!is.na(x)])
      return(length(unique_vals) <= 10)  # Treat as categorical if 10 or fewer unique values
    }
    return(FALSE)
  }
  
  # Function to compute cell type proportions
  compute_proportions <- function(spe_obj, column) {
    col_data <- colData(spe_obj)
    if(!column %in% colnames(col_data)) return(NULL)
    
    cell_types <- col_data[[column]]
    cell_types <- cell_types[!is.na(cell_types)]
    
    if(length(cell_types) == 0) return(NULL)
    
    prop_table <- table(cell_types)
    proportions <- prop_table / sum(prop_table)
    
    return(data.frame(
      CellType = names(proportions),
      Proportion = as.numeric(proportions),
      Count = as.numeric(prop_table),
      stringsAsFactors = FALSE
    ))
  }
  
  # Function to create categorical plots
  create_categorical_plot <- function(data_list, column, labels, comparison_type, selected_celltype = NULL) {
    if(comparison_type == "single_celltype" && is.null(selected_celltype)) {
      return(NULL)
    }
    
    combined_data <- data.frame()
    
    for(i in seq_along(data_list)) {
      spe_obj <- data_list[[i]]
      prop_data <- compute_proportions(spe_obj, column)
      
      if(!is.null(prop_data)) {
        if(comparison_type == "single_celltype") {
          # Filter for selected cell type
          celltype_data <- prop_data[prop_data$CellType == selected_celltype, ]
          if(nrow(celltype_data) > 0) {
            temp_data <- data.frame(
              Sample = labels[i],
              CellType = selected_celltype,
              Proportion = celltype_data$Proportion,
              Count = celltype_data$Count,
              stringsAsFactors = FALSE
            )
            combined_data <- rbind(combined_data, temp_data)
          } else {
            # Cell type not found in this sample, add with 0 proportion
            temp_data <- data.frame(
              Sample = labels[i],
              CellType = selected_celltype,
              Proportion = 0,
              Count = 0,
              stringsAsFactors = FALSE
            )
            combined_data <- rbind(combined_data, temp_data)
          }
        } else {
          # All cell types
          temp_data <- data.frame(
            Sample = labels[i],
            CellType = prop_data$CellType,
            Proportion = prop_data$Proportion,
            Count = prop_data$Count,
            stringsAsFactors = FALSE
          )
          combined_data <- rbind(combined_data, temp_data)
        }
      }
    }
    
    if(nrow(combined_data) == 0) return(NULL)
    
    if(comparison_type == "single_celltype") {
      # Bar plot comparing one cell type across samples
      p <- ggplot(combined_data, aes(x = Sample, y = Proportion, fill = Sample)) +
        geom_bar(stat = "identity", alpha = 0.7) +
        geom_text(aes(label = paste0(round(Proportion * 100, 1), "%")), 
                  vjust = -0.5, size = 3) +
        labs(title = paste("Proportion of", selected_celltype, "across samples"),
             x = "Sample",
             y = "Proportion") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "none") +
        ylim(0, max(combined_data$Proportion) * 1.1)
      
    } else {
      # Stacked bar plot showing all cell types within each sample
      p <- ggplot(combined_data, aes(x = Sample, y = Proportion, fill = CellType)) +
        geom_bar(stat = "identity", position = "stack", alpha = 0.8) +
        labs(title = paste("Cell type proportions by sample -", column),
             x = "Sample",
             y = "Proportion",
             fill = "Cell Type") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        scale_y_continuous(labels = scales::percent_format())
    }
    
    return(p)
  }
  
  # Function to create bar plot for binary data (showing means)
  create_barplot <- function(data_list, column, labels) {
    combined_data <- data.frame()
    
    for(i in seq_along(data_list)) {
      spe_obj <- data_list[[i]]
      col_data <- colData(spe_obj)
      
      if(column %in% colnames(col_data)) {
        # Convert to numeric for mean calculation
        values <- as.numeric(col_data[[column]])
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
      spe_obj <- data_list[[i]]
      col_data <- colData(spe_obj)
      
      if(column %in% colnames(col_data)) {
        temp_data <- data.frame(
          Value = as.numeric(col_data[[column]]),
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
    
    # Create temporary directory for plots
    temp_dir <- tempdir()
    plot_dir <- file.path(temp_dir, "plots")
    if(dir.exists(plot_dir)) unlink(plot_dir, recursive = TRUE)
    dir.create(plot_dir, recursive = TRUE)
    
    plot_files <- character()
    
    withProgress(message = "Generating plots...", value = 0, {
      
      for(i in seq_along(input$selected_columns)) {
        column <- input$selected_columns[i]
        
        incProgress(1/length(input$selected_columns), 
                    detail = paste("Processing", column))
        
        # Check if this column is the selected categorical column
        if(!is.null(input$categorical_column) && column == input$categorical_column) {
          # Handle categorical analysis
          plot_obj <- create_categorical_plot(
            values$data_list, 
            column, 
            labels, 
            input$comparison_type, 
            input$selected_celltype
          )
          
          if(!is.null(plot_obj)) {
            if(input$comparison_type == "single_celltype") {
              file_name <- paste0(gsub("[^A-Za-z0-9]", "_", column), "_", 
                                  gsub("[^A-Za-z0-9]", "_", input$selected_celltype), 
                                  "_proportion.png")
            } else {
              file_name <- paste0(gsub("[^A-Za-z0-9]", "_", column), "_all_celltypes_proportion.png")
            }
            
            file_path <- file.path(plot_dir, file_name)
            
            tryCatch({
              ggsave(file_path, plot_obj, width = 12, height = 8, dpi = 300)
              plot_files <- c(plot_files, file_path)
            }, error = function(e) {
              showNotification(paste("Error saving categorical plot for", column, ":", e$message), 
                               type = "warning")
            })
          }
          
        } else {
          # Handle non-categorical analysis (existing logic)
          # Check if column is binary or continuous
          sample_data <- NULL
          for(j in seq_along(values$data_list)) {
            spe_obj <- values$data_list[[j]]
            col_data <- colData(spe_obj)
            if(column %in% colnames(col_data)) {
              sample_data <- c(sample_data, col_data[[column]])
            }
          }
          
          if(is.null(sample_data)) next
          
          plot_obj <- NULL
          file_suffix <- ""
          
          if(is_binary(sample_data)) {
            plot_obj <- create_barplot(values$data_list, column, labels)
            file_suffix <- "_barplot.png"
          } else if(!is_categorical(sample_data)) {
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
      }
    })
    
    values$plot_files <- plot_files
    values$plots_generated <- length(plot_files) > 0
    
    if(values$plots_generated) {
      showNotification(paste("Generated", length(plot_files), "plots successfully!"), 
                       type = "default")
    } else {
      showNotification("No plots could be generated. Please check your data and column selections.", 
                       type = "warning")
    }
  })
  
  # Download handler
  output$download_plots <- downloadHandler(
    filename = function() {
      paste0("spatial_analysis_plots_", Sys.Date(), ".zip")
    },
    content = function(file) {
      req(values$plot_files)
      
      # Create zip file
      zip::zipr(zipfile = file, 
                files = values$plot_files,
                recurse = FALSE)
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