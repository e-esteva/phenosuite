library(shiny)
library(ggplot2)
library(DT)
library(zip)
library(data.table)
library(ggrastr)

RASTER_THRESHOLD <- 150000
MAX_DISPLAY_POINTS <- 750000

# Load sp package for point.in.polygon or define custom function
if (!requireNamespace("sp", quietly = TRUE)) {
  point.in.polygon <- function(point.x, point.y, pol.x, pol.y) {
    n <- length(pol.x)
    inside <- rep(FALSE, length(point.x))
    
    for (i in seq_along(point.x)) {
      count <- 0
      for (j in 1:n) {
        k <- j %% n + 1
        if (((pol.y[j] <= point.y[i]) && (point.y[i] < pol.y[k])) ||
            ((pol.y[k] <= point.y[i]) && (point.y[i] < pol.y[j]))) {
          x_intersect <- (pol.x[k] - pol.x[j]) * (point.y[i] - pol.y[j]) / 
                         (pol.y[k] - pol.y[j]) + pol.x[j]
          if (point.x[i] < x_intersect) {
            count <- count + 1
          }
        }
      }
      inside[i] <- (count %% 2) == 1
    }
    
    return(as.integer(inside))
  }
} else {
  library(sp)
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

# Helper function to blend colors for overlaps
blend_colors <- function(colors) {
  if (length(colors) == 1) return(colors)
  rgb_vals <- col2rgb(colors)
  avg_rgb <- rowMeans(rgb_vals)
  rgb(avg_rgb[1], avg_rgb[2], avg_rgb[3], maxColorValue = 255)
}

# Helper to find matching gate column in data
find_gate_column <- function(gate_name, col_names) {
  if (gate_name %in% col_names) return(gate_name)
  safe_name <- make.names(gate_name)
  if (safe_name %in% col_names) return(safe_name)
  return(NULL)
}

# Resolve a variable name from JSON against actual column names
# Handles hyphen <-> dot, make.names mangling, etc.
resolve_column_name <- function(var_name, col_names) {
  if (var_name %in% col_names) return(var_name)
  safe <- make.names(var_name)
  if (safe %in% col_names) return(safe)
  for (cn in col_names) {
    if (make.names(cn) == var_name || make.names(cn) == safe) return(cn)
  }
  return(NULL)
}

# Apply a comparison operator to a vector
apply_operator <- function(vals, op, threshold, threshold_max = NULL) {
  result <- switch(op,
    "gt"      = vals > threshold,
    "gte"     = vals >= threshold,
    "lt"      = vals < threshold,
    "lte"     = vals <= threshold,
    "between" = vals >= threshold & vals <= threshold_max
  )
  result[is.na(result)] <- FALSE
  result
}

# Convert operator + threshold(s) to effective min/max for visualization
operator_to_range <- function(op, threshold, threshold_max, data_range) {
  switch(op,
    "gt"      = list(min = threshold, max = data_range[2]),
    "gte"     = list(min = threshold, max = data_range[2]),
    "lt"      = list(min = data_range[1], max = threshold),
    "lte"     = list(min = data_range[1], max = threshold),
    "between" = list(min = threshold, max = threshold_max)
  )
}

shinyServer(function(input, output, session) {
  
  # Reactive values to store data and gates
  rv <- reactiveValues(
    original_data = NULL,
    current_data = NULL,
    full_data_with_gates = NULL,
    gate_number = 0,
    gate_tree = list(),
    current_gate_id = NULL,
    saved_brush = NULL,
    gate_data_cache = list(),
    gate_colors = c(),
    visible_gates = c(),
    x_range = NULL,
    y_range = NULL
  )

  # Debug output to check temp directory settings
  output$debug <- renderText({
    paste(
      "user:", Sys.info()[["user"]],
      "\npid: ", Sys.getpid(),
      "\nTMPDIR:", Sys.getenv("TMPDIR"),
      "\nTMP:", Sys.getenv("TMP"),
      "\nTEMP:", Sys.getenv("TEMP"),
      "\nR_SESSION_TMPDIR:", Sys.getenv("R_SESSION_TMPDIR"),
      "\ntempdir():", tempdir()
    )
  })

  output$tmp_env <- renderPrint({
    list(
      pid = Sys.getpid(),
      user = Sys.info()[["user"]],
      tempdir = tempdir(),
      env = Sys.getenv(c("TMPDIR","TMP","TEMP","R_SESSION_TMPDIR"))
    )
  })

  probe_res <- eventReactive(input$run_tmp_probe, {
    td <- Sys.getenv("TMPDIR")
    p  <- file.path(td, paste0("probe_", Sys.getpid(), "_", as.integer(Sys.time())))
    f  <- file.path(p, "x.txt")
    list(
      pid = Sys.getpid(),
      user = Sys.info()[["user"]],
      TMPDIR = td,
      tempdir = tempdir(),
      mkdir = tryCatch({ dir.create(p, recursive = TRUE); TRUE }, error=function(e) e$message),
      write = tryCatch({ writeLines("ok", f); TRUE }, error=function(e) e$message),
      cleanup = tryCatch({ unlink(p, recursive = TRUE); TRUE }, error=function(e) e$message)
    )
  })

  output$tmp_probe <- renderPrint({
    req(probe_res())
    probe_res()
  })

  ###################
  # Load data automatically when file is selected
  observeEvent(input$csv_file, {
    req(input$csv_file)
    
    tryCatch({
      data <- data.table::fread(input$csv_file$datapath, stringsAsFactors = FALSE)
      data <- as.data.frame(data)
      rownames(data) <- as.character(seq_len(nrow(data)))
      
      rv$original_data <- data
      rv$current_data <- data
      rv$full_data_with_gates <- data
      rv$gate_number <- 0
      rv$gate_tree <- list()
      rv$current_gate_id <- "root"
      rv$saved_brush <- NULL
      rv$gate_data_cache <- list(root = data)
      rv$gate_colors <- c()
      rv$visible_gates <- c()
      
      updateTextInput(session, "gate_name", value = "Gate_1")
      showNotification("Data loaded successfully!", type = "message")
    }, error = function(e) {
      showNotification(paste("Error loading data:", e$message), type = "error")
    })
  })
  
  # Data preview
  output$data_preview_box <- renderUI({
    req(rv$original_data)
    box(title = "Data Preview", width = 12, DTOutput("data_preview"))
  })
  
  output$data_preview <- renderDT({
    req(rv$original_data)
    datatable(head(rv$original_data, 100), options = list(scrollX = TRUE))
  })
  
  # Dynamic UI for axis selection — numeric columns only, gate columns excluded
  output$x_axis_select <- renderUI({
    req(rv$current_data)
    all_cols <- setdiff(names(rv$current_data), 
                        names(rv$current_data)[grepl("[. ]Gate$", names(rv$current_data))])
    original_cols <- all_cols[sapply(all_cols, function(col) is.numeric(rv$current_data[[col]]))]
    req(length(original_cols) > 0)
    selectInput("x_axis", "X-axis Variable:", 
                choices = original_cols,
                selected = original_cols[1])
  })
  
  output$y_axis_select <- renderUI({
    req(rv$current_data)
    all_cols <- setdiff(names(rv$current_data),
                        names(rv$current_data)[grepl("[. ]Gate$", names(rv$current_data))])
    original_cols <- all_cols[sapply(all_cols, function(col) is.numeric(rv$current_data[[col]]))]
    req(length(original_cols) > 0)
    selectInput("y_axis", "Y-axis Variable:", 
                choices = original_cols,
                selected = original_cols[min(2, length(original_cols))])
  })

  # Store brush coordinates
  observe({
    if (!is.null(input$plot_brush)) {
      rv$saved_brush <- list(
        xmin = input$plot_brush$xmin,
        xmax = input$plot_brush$xmax,
        ymin = input$plot_brush$ymin,
        ymax = input$plot_brush$ymax
      )
    }
  })
  
  # Handle gate selection from tree
  observeEvent(input$select_gate, {
    gate_id <- input$select_gate
    if (!is.null(gate_id) && gate_id != rv$current_gate_id && gate_id %in% names(rv$gate_data_cache)) {
      rv$current_gate_id <- gate_id
      rv$current_data <- rv$gate_data_cache[[gate_id]]
      rv$saved_brush <- NULL
      session$resetBrush("plot_brush")
      
      children <- get_children(rv$gate_tree, gate_id)
      next_num <- length(children) + 1
      if (gate_id == "root") {
        updateTextInput(session, "gate_name", value = paste0("Gate_", next_num))
      } else {
        updateTextInput(session, "gate_name", value = paste0(gate_id, ".", next_num))
      }
    }
  }, ignoreInit = TRUE)
  
  ############################################
  # Restore gates from JSON schema
  ############################################
  observeEvent(input$restore_gates, {
    req(rv$original_data, input$gate_metadata_file)
    
    tryCatch({
      metadata <- jsonlite::fromJSON(input$gate_metadata_file$datapath)
      data_cols <- names(rv$original_data)
      
      # --- Phase 1: Validate all variables exist before touching any state ---
      missing_vars <- c()
      col_mapping <- list()
      
      for (gate_name in names(metadata)) {
        gate_info <- metadata[[gate_name]]
        
        for (var_field in c("x_var", "y_var")) {
          json_var <- gate_info[[var_field]]
          if (is.null(json_var)) {
            missing_vars <- c(missing_vars, paste0(gate_name, ": ", var_field, " not specified"))
            next
          }
          if (!json_var %in% names(col_mapping)) {
            resolved <- resolve_column_name(json_var, data_cols)
            if (is.null(resolved)) {
              missing_vars <- c(missing_vars, json_var)
            } else {
              col_mapping[[json_var]] <- resolved
            }
          }
        }
      }
      
      if (length(missing_vars) > 0) {
        unique_missing <- unique(missing_vars)
        output$schema_validation_msg <- renderUI({
          div(
            style = "color: red; padding: 10px; background: #fff0f0; border: 1px solid red; border-radius: 4px;",
            h5("Schema validation failed — variables not found in data:"),
            tags$ul(lapply(unique_missing, tags$li)),
            p(em("Available numeric columns:"), 
              paste(data_cols[sapply(data_cols, function(c) is.numeric(rv$original_data[[c]]))], 
                    collapse = ", "))
          )
        })
        return()
      }
      
      # --- Phase 2: Apply gates hierarchically with transformations ---
      rv$gate_tree <- metadata
      rv$gate_data_cache <- list(root = rv$original_data)
      rv$gate_number <- 0
      rv$gate_colors <- c()
      rv$visible_gates <- c()
      rv$full_data_with_gates <- rv$original_data
      
      color_palette <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
                         "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5",
                         "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F",
                         "#E5C494", "#B3B3B3", "#1B9E77", "#D95F02", "#7570B3")
      
      for (gate_name in names(metadata)) {
        gate_info <- metadata[[gate_name]]
        
        parent_id <- gate_info$parent
        if (!(parent_id %in% names(rv$gate_data_cache))) {
          stop(paste("Parent gate", parent_id, "not found for", gate_name,
                     "— gates in JSON may be out of order"))
        }
        parent_data <- rv$gate_data_cache[[parent_id]]
        
        # Resolve actual column names
        x_col <- col_mapping[[gate_info$x_var]]
        y_col <- col_mapping[[gate_info$y_var]]
        
        # Get raw values and apply stored transformation
        x_raw <- parent_data[[x_col]]
        y_raw <- parent_data[[y_col]]
        
        x_transform <- gate_info$x_transform %||% "none"
        y_transform <- gate_info$y_transform %||% "none"
        
        x_vals <- apply_transformation(x_raw, x_transform)
        y_vals <- apply_transformation(y_raw, y_transform)
        
        # Apply thresholds — use operators if stored, fall back to between
        x_op <- gate_info$x_operator %||% "between"
        y_op <- gate_info$y_operator %||% "between"
        
        x_inside <- apply_operator(x_vals, x_op,
                                   gate_info$x_threshold %||% gate_info$threshold_xmin,
                                   gate_info$x_threshold_max %||% gate_info$threshold_xmax)
        y_inside <- apply_operator(y_vals, y_op,
                                   gate_info$y_threshold %||% gate_info$threshold_ymin,
                                   gate_info$y_threshold_max %||% gate_info$threshold_ymax)
        
        inside <- x_inside & y_inside
        inside[is.na(inside)] <- FALSE
        
        # Cache the gated data
        gated_data <- parent_data[inside, ]
        rv$gate_data_cache[[gate_name]] <- gated_data
        
        # Update cell counts to reflect this dataset
        gate_info$n_before <- nrow(parent_data)
        gate_info$n_after <- sum(inside)
        rv$gate_tree[[gate_name]] <- gate_info
        
        # Update full dataset gate column
        matching_col <- find_gate_column(gate_name, names(rv$full_data_with_gates))
        if (is.null(matching_col)) matching_col <- make.names(gate_name)
        
        if (!(matching_col %in% names(rv$full_data_with_gates))) {
          rv$full_data_with_gates[[matching_col]] <- 0
        }
        rows_that_pass <- rownames(parent_data)[inside]
        rv$full_data_with_gates[rows_that_pass, matching_col] <- 1
        
        # Assign color
        rv$gate_number <- rv$gate_number + 1
        rv$gate_colors[gate_name] <- color_palette[(rv$gate_number - 1) %% length(color_palette) + 1]
        rv$visible_gates <- c(rv$visible_gates, gate_name)
      }
      
      rv$current_gate_id <- "root"
      rv$current_data <- rv$original_data
      
      output$schema_validation_msg <- renderUI({
        div(
          style = "color: green; padding: 10px; background: #f0fff0; border: 1px solid green; border-radius: 4px;",
          h5(sprintf("Successfully applied %d gates!", length(metadata)))
        )
      })
      
      showNotification(sprintf("Applied %d gates from schema!", length(metadata)), 
                       type = "message")
      
    }, error = function(e) {
      showNotification(paste("Error applying schema:", e$message), type = "error")
      output$schema_validation_msg <- renderUI({
        div(style = "color: red;", p(strong("Error:"), e$message))
      })
    })
  })

  ############################################
  # Numeric gate preview
  ############################################
  numeric_gate_mask <- reactive({
    req(input$numeric_gate_yn == "yes")
    req(transformed_current_data(), input$x_operator, input$y_operator)
    req(!is.null(input$x_threshold), !is.null(input$y_threshold))
    
    data <- transformed_current_data()
    
    x_inside <- apply_operator(data$x_transformed, input$x_operator,
                               input$x_threshold, input$x_threshold_max)
    y_inside <- apply_operator(data$y_transformed, input$y_operator,
                               input$y_threshold, input$y_threshold_max)
    
    x_inside & y_inside
  })

  output$numeric_gate_preview_stats <- renderUI({
    req(numeric_gate_mask())
    mask <- numeric_gate_mask()
    n_total <- length(mask)
    n_pass <- sum(mask)
    div(
      style = "padding: 8px; background: #f0f0f0; border-radius: 4px; margin-top: 5px;",
      p(strong("Preview:"), sprintf("%d / %d cells (%.1f%%)", 
                                    n_pass, n_total, 100 * n_pass / n_total))
    )
  })

  ############################################
  # Apply gate — brush or numeric mode
  ############################################
  observeEvent(input$apply_gate, {
    req(rv$current_data, input$x_axis, input$y_axis)
    req(input$gate_name != "")
    req(transformed_current_data())
    
    use_numeric <- isTRUE(input$numeric_gate_yn == "yes")
    
    data <- transformed_current_data()
    x_vals <- data$x_transformed
    y_vals <- data$y_transformed
    
    if (use_numeric) {
      # --- Numeric threshold mode ---
      req(input$x_operator, input$y_operator)
      req(!is.null(input$x_threshold), !is.null(input$y_threshold))
      
      x_inside <- apply_operator(x_vals, input$x_operator,
                                 input$x_threshold, input$x_threshold_max)
      y_inside <- apply_operator(y_vals, input$y_operator,
                                 input$y_threshold, input$y_threshold_max)
      inside <- x_inside & y_inside
      
      x_range <- range(x_vals, na.rm = TRUE)
      y_range <- range(y_vals, na.rm = TRUE)
      x_bounds <- operator_to_range(input$x_operator, input$x_threshold,
                                    input$x_threshold_max, x_range)
      y_bounds <- operator_to_range(input$y_operator, input$y_threshold,
                                    input$y_threshold_max, y_range)
      
      eff_xmin <- x_bounds$min
      eff_xmax <- x_bounds$max
      eff_ymin <- y_bounds$min
      eff_ymax <- y_bounds$max
      
    } else {
      # --- Brush mode ---
      req(rv$saved_brush)
      brush <- rv$saved_brush
      inside <- x_vals >= brush$xmin & x_vals <= brush$xmax &
                y_vals >= brush$ymin & y_vals <= brush$ymax
      
      eff_xmin <- brush$xmin
      eff_xmax <- brush$xmax
      eff_ymin <- brush$ymin
      eff_ymax <- brush$ymax
    }
    
    # --- Common gate creation logic ---
    raw_name <- trimws(input$gate_name)
    if (!grepl("gate", raw_name, ignore.case = TRUE)) {
      gate_name <- paste0(raw_name, " Gate")
    } else {
      gate_name <- raw_name
    }
    gate_name <- make.names(gate_name, unique = FALSE)
    
    if (gate_name %in% names(rv$gate_tree)) {
      showNotification("Gate name already exists! Please choose a different name.", 
                       type = "error")
      return()
    }
    
    rv$gate_number <- rv$gate_number + 1
    
    color_palette <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", 
                       "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5",
                       "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F",
                       "#E5C494", "#B3B3B3", "#1B9E77", "#D95F02", "#7570B3")
    gate_color <- color_palette[(rv$gate_number - 1) %% length(color_palette) + 1]
    rv$gate_colors[gate_name] <- gate_color
    
    gate_info <- list(
      id = gate_name,
      parent = rv$current_gate_id,
      x_var = input$x_axis,
      y_var = input$y_axis,
      x_transform = input$x_scale %||% "none",
      y_transform = input$y_scale %||% "none",
      n_before = nrow(rv$current_data),
      n_after = sum(inside),
      children = list(),
      gate_mode = if (use_numeric) "numeric" else "brush",
      x_operator = if (use_numeric) input$x_operator else "between",
      y_operator = if (use_numeric) input$y_operator else "between",
      x_threshold = if (use_numeric) input$x_threshold else NULL,
      x_threshold_max = if (use_numeric) input$x_threshold_max else NULL,
      y_threshold = if (use_numeric) input$y_threshold else NULL,
      y_threshold_max = if (use_numeric) input$y_threshold_max else NULL,
      threshold_xmin = eff_xmin,
      threshold_xmax = eff_xmax,
      threshold_ymin = eff_ymin,
      threshold_ymax = eff_ymax
    )
    
    rv$gate_tree[[gate_name]] <- gate_info
    
    matching_col <- find_gate_column(gate_name, names(rv$full_data_with_gates))
    if (is.null(matching_col)) matching_col <- make.names(gate_name)
    
    if (!(matching_col %in% names(rv$full_data_with_gates))) {
      rv$full_data_with_gates[[matching_col]] <- 0
    }
    
    rows_that_pass <- rownames(rv$current_data)[inside]
    rv$full_data_with_gates[rows_that_pass, matching_col] <- 1
    
    gated_data <- rv$current_data[inside, ]
    rv$gate_data_cache[[gate_name]] <- gated_data
    
    # Store snapshot plot
    brush_rect <- list(xmin = eff_xmin, xmax = eff_xmax, 
                       ymin = eff_ymin, ymax = eff_ymax)
    temp_plot <- create_gate_plot(data, "x_transformed", "y_transformed", 
                                  brush_rect, gate_name, inside,
                                  get_transform_label(input$x_axis, input$x_scale %||% "none"),
                                  get_transform_label(input$y_axis, input$y_scale %||% "none"))
    gate_info$plot <- temp_plot
    rv$gate_tree[[gate_name]] <- gate_info
    
    rv$current_data <- gated_data
    rv$current_gate_id <- gate_name
    rv$visible_gates <- c(rv$visible_gates, gate_name)
    
    updateTextInput(session, "gate_name", value = paste0(gate_name, ".1"))
    
    showNotification(paste0("Gate '", gate_name, "' applied! Cells remaining: ", sum(inside)), 
                     type = "message")
  })

  # Reset current brush/selection
  observeEvent(input$reset_gate, {
    rv$saved_brush <- NULL
    session$resetBrush("plot_brush")
    showNotification("Selection cleared.", type = "message")
  })
  
  # Reset all gates
  observeEvent(input$reset_all, {
    rv$current_data <- rv$original_data
    rv$full_data_with_gates <- rv$original_data
    rv$gate_number <- 0
    rv$gate_tree <- list()
    rv$current_gate_id <- "root"
    rv$saved_brush <- NULL
    rv$gate_data_cache <- list(root = rv$original_data)
    rv$gate_colors <- c()
    rv$visible_gates <- c()
    session$resetBrush("plot_brush")
    updateTextInput(session, "gate_name", value = "Gate_1")
    showNotification("All gates reset. Starting fresh.", type = "warning")
  })
  
  ############################################
  # Gate tree UI
  ############################################
  output$gate_tree <- renderUI({
    req(rv$original_data)
    gate_count <- length(rv$gate_tree)
    
    if (gate_count == 0) {
      return(div(
        radioButtons("select_gate", NULL, 
                     choices = c("Root (All cells)" = "root"),
                     selected = "root")
      ))
    }
    
    tree_html <- build_gate_tree_ui(rv$gate_tree, rv$current_gate_id)
    div(
      radioButtons("select_gate", NULL, 
                   choices = tree_html$choices,
                   selected = rv$current_gate_id)
    )
  })
  
  build_gate_tree_ui <- function(gate_tree, current_id) {
    choices <- c("Root (All cells)" = "root")
    root_gates <- names(gate_tree)[sapply(gate_tree, function(g) g$parent == "root")]
    
    for (gate_id in root_gates) {
      gate_info <- gate_tree[[gate_id]]
      label <- sprintf("%s (%d\u2192%d, %.1f%%)", 
                       gate_id, gate_info$n_before, gate_info$n_after,
                       100 * gate_info$n_after / gate_info$n_before)
      choices[label] <- gate_id
      
      children <- get_children(gate_tree, gate_id)
      if (length(children) > 0) {
        child_info <- build_children_ui(gate_tree, gate_id, 1)
        choices <- c(choices, child_info$choices)
      }
    }
    list(choices = choices, html = "")
  }
  
  build_children_ui <- function(gate_tree, parent_id, level) {
    children <- get_children(gate_tree, parent_id)
    choices <- c()
    
    for (gate_id in children) {
      gate_info <- gate_tree[[gate_id]]
      indent <- paste(rep("  ", level), collapse = "")
      label <- sprintf("%s\u21b3 %s (%d\u2192%d, %.1f%%)", 
                       indent, gate_id, gate_info$n_before, gate_info$n_after,
                       100 * gate_info$n_after / gate_info$n_before)
      choices[label] <- gate_id
      
      grandchildren <- get_children(gate_tree, gate_id)
      if (length(grandchildren) > 0) {
        child_info <- build_children_ui(gate_tree, gate_id, level + 1)
        choices <- c(choices, child_info$choices)
      }
    }
    list(choices = choices)
  }
  
  get_children <- function(gate_tree, parent_id) {
    names(gate_tree)[sapply(gate_tree, function(g) g$parent == parent_id)]
  }
  
  ############################################
  # Gate display checkboxes — isolate to break reactive loop
  ############################################
  output$gate_display_checkboxes <- renderUI({
    all_children <- get_all_descendants(rv$gate_tree, rv$current_gate_id)
    
    if (length(all_children) == 0) {
      return(p("No child gates to display"))
    }
    
    current_visible <- isolate(rv$visible_gates)
    valid_selected <- intersect(current_visible, all_children)
    
    checkboxGroupInput("visible_gates_check",
                       "Show gates on plot:",
                       choices = setNames(all_children, all_children),
                       selected = valid_selected)
  })
  
  observe({
    if (!is.null(input$visible_gates_check)) {
      rv$visible_gates <- input$visible_gates_check
    }
  })

  output$current_gate_info <- renderUI({
    if (rv$current_gate_id == "root") {
      div(
        style = "padding: 10px; background-color: #f0f0f0; margin-bottom: 10px;",
        h4("Viewing: Root Population (All Cells)"),
        p(sprintf("Total cells: %d", nrow(rv$current_data)))
      )
    } else {
      gate_info <- rv$gate_tree[[rv$current_gate_id]]
      div(
        style = "padding: 10px; background-color: #f0f0f0; margin-bottom: 10px;",
        h4(sprintf("Viewing: %s", rv$current_gate_id)),
        p(sprintf("%s vs %s | %d cells (%.1f%% of parent)",
                  gate_info$x_var, gate_info$y_var,
                  gate_info$n_after,
                  100 * gate_info$n_after / gate_info$n_before))
      )
    }
  })
  
  ############################################
  # Main scatter plot
  ############################################
  output$gate_plot <- renderPlot({
    req(rv$current_data, input$x_axis, input$y_axis)
    req(transformed_current_data())
    
    x_scale <- input$x_scale %||% "none"
    y_scale <- input$y_scale %||% "none"
    x_label <- get_transform_label(input$x_axis, x_scale)
    y_label <- get_transform_label(input$y_axis, y_scale)
    
    plot_data <- transformed_current_data()
    
    # Determine which points are inside saved brush (in transformed space)
    inside_brush <- rep(FALSE, nrow(plot_data))
    if (!is.null(rv$saved_brush) && isTRUE(input$numeric_gate_yn == "no")) {
      brush <- rv$saved_brush
      inside_brush <- plot_data$x_transformed >= brush$xmin & 
                      plot_data$x_transformed <= brush$xmax &
                      plot_data$y_transformed >= brush$ymin & 
                      plot_data$y_transformed <= brush$ymax
    }
    
    # Numeric gate preview highlighting
    if (isTRUE(input$numeric_gate_yn == "yes")) {
      tryCatch({
        mask <- numeric_gate_mask()
        if (!is.null(mask)) inside_brush <- mask
      }, error = function(e) {})
    }
    
    plot_data$selected <- inside_brush
    
    # Gate membership logic — vectorized
    current_row_names <- rownames(plot_data)
    all_children <- get_all_descendants(rv$gate_tree, rv$current_gate_id)
    visible_children <- intersect(all_children, rv$visible_gates)
    
    if (length(visible_children) > 0) {
      gate_member_matrix <- matrix(FALSE, 
                                    nrow = nrow(plot_data), 
                                    ncol = length(visible_children),
                                    dimnames = list(current_row_names, visible_children))
      
      for (gate_id in visible_children) {
        matching_col <- find_gate_column(gate_id, names(rv$full_data_with_gates))
        if (!is.null(matching_col)) {
          passing_rows <- rownames(rv$full_data_with_gates)[rv$full_data_with_gates[[matching_col]] == 1]
          gate_member_matrix[current_row_names %in% passing_rows, gate_id] <- TRUE
        }
      }
      
      gate_depths <- sapply(visible_children, function(g) gate_depth(rv$gate_tree, g))
      gate_parents <- sapply(visible_children, function(g) rv$gate_tree[[g]]$parent)
      
      gate_memberships <- character(nrow(plot_data))
      gate_memberships[] <- "Ungated"
      
      for (i in seq_len(nrow(plot_data))) {
        member_gates <- visible_children[gate_member_matrix[i, ]]
        
        if (length(member_gates) > 0) {
          depths <- gate_depths[member_gates]
          max_depth <- max(depths)
          deepest <- member_gates[depths == max_depth]
          
          if (length(deepest) > 1) {
            parents <- gate_parents[deepest]
            siblings <- deepest[parents == parents[1]]
            if (length(siblings) > 1) {
              gate_memberships[i] <- paste(sort(siblings), collapse = "+")
            } else {
              gate_memberships[i] <- deepest[1]
            }
          } else {
            gate_memberships[i] <- deepest[1]
          }
        }
      }
      
      plot_data$gate_membership <- gate_memberships
    } else {
      plot_data$gate_membership <- "Ungated"
    }
    
    # Build color map
    unique_classes <- unique(plot_data$gate_membership)
    gate_colors_map <- c("Ungated" = "gray80")
    
    for (gate_id in visible_children) {
      if (gate_id %in% names(rv$gate_colors)) {
        gate_colors_map[gate_id] <- rv$gate_colors[gate_id]
      }
    }
    
    overlap_classes <- unique_classes[grepl("\\+", unique_classes)]
    for (ov in overlap_classes) {
      component_gates <- unlist(strsplit(ov, "\\+"))
      parent_colors <- rv$gate_colors[component_gates]
      gate_colors_map[ov] <- blend_colors(parent_colors)
    }
    
    plot_data$display_color <- plot_data$gate_membership
    plot_data$display_color[inside_brush] <- "Selected"
    gate_colors_map["Selected"] <- "red"
    
    class_counts <- table(plot_data$gate_membership)
    
    if (isTRUE(input$show_counts_in_legend)) {
      legend_labels <- names(gate_colors_map)
      legend_labels_with_counts <- paste0(
        legend_labels, " (", class_counts[legend_labels] %||% 0, ")")
      names(gate_colors_map) <- legend_labels_with_counts
      
      plot_data$display_color <- ifelse(
        plot_data$display_color == "Selected",
        "Selected",
        paste0(plot_data$gate_membership, " (", class_counts[plot_data$gate_membership], ")")
      )
      
      plot_data$display_color <- factor(
        plot_data$display_color,
        levels = c(legend_labels_with_counts, "Selected")
      )
    } else {
      plot_data$display_color <- factor(
        plot_data$display_color,
        levels = names(gate_colors_map)
      )
    }
    
    # Downsample for display
    render_data <- plot_data
    if (nrow(render_data) > MAX_DISPLAY_POINTS) {
      set.seed(1)
      render_data <- render_data[sample(nrow(render_data), MAX_DISPLAY_POINTS), ]
    }
    
    use_raster <- nrow(render_data) > RASTER_THRESHOLD
    
    # Draw ungated first (behind), gated on top (visible)
    render_data <- render_data[order(render_data$gate_membership == "Ungated", decreasing = TRUE), ]
    
    # Dynamic point size
    point_size <- max(0.3, min(4, 6 / log10(nrow(render_data) + 1)))
    
    p <- ggplot(render_data, aes(x = x_transformed, y = y_transformed))
    
    if (use_raster) {
      p <- p + geom_point_rast(aes(color = display_color), alpha = 1, size = point_size)
    } else {
      p <- p + geom_point(aes(color = display_color), alpha = 1, size = point_size)
    }
    
    p <- p + scale_color_manual(values = gate_colors_map, name = "Gate", drop = FALSE) +
      labs(x = x_label, y = y_label) +
      theme_bw() +
      theme(text = element_text(size = 14),
            legend.position = "right",
            legend.background = element_rect(fill = "white", color = "black"),
            legend.key.size = unit(0.8, "lines"))
    
    # Apply axis limits from sliders (only in slider mode)
    if (isTRUE(input$numeric_gate_yn == "no")) {
      if (!is.null(input$x_range_slider)) p <- p + xlim(input$x_range_slider)
      if (!is.null(input$y_range_slider)) p <- p + ylim(input$y_range_slider)
    }
    
    # Add gate boundary visualization
    if (isTRUE(input$numeric_gate_yn == "yes") && 
        !is.null(input$x_operator) && !is.null(input$y_operator) &&
        !is.null(input$x_threshold) && !is.null(input$y_threshold)) {
      
      x_range_data <- range(render_data$x_transformed, na.rm = TRUE)
      y_range_data <- range(render_data$y_transformed, na.rm = TRUE)
      x_bounds <- operator_to_range(input$x_operator, input$x_threshold,
                                    input$x_threshold_max, x_range_data)
      y_bounds <- operator_to_range(input$y_operator, input$y_threshold,
                                    input$y_threshold_max, y_range_data)
      
      # Shaded region showing what passes
      p <- p + annotate("rect",
                        xmin = x_bounds$min, xmax = x_bounds$max,
                        ymin = y_bounds$min, ymax = y_bounds$max,
                        fill = "red", alpha = 0.08,
                        color = "red", linewidth = 1.2, linetype = "dashed")
      
      # Threshold lines — solid = inclusive, dotted = strict
      x_op <- input$x_operator
      y_op <- input$y_operator
      
      if (x_op %in% c("gt", "gte")) {
        p <- p + geom_vline(xintercept = input$x_threshold, color = "red", 
                            linewidth = 1, linetype = if (x_op == "gt") "dotted" else "solid")
      } else if (x_op %in% c("lt", "lte")) {
        p <- p + geom_vline(xintercept = input$x_threshold, color = "red",
                            linewidth = 1, linetype = if (x_op == "lt") "dotted" else "solid")
      } else if (x_op == "between") {
        p <- p + geom_vline(xintercept = input$x_threshold, color = "red", linewidth = 1) +
                 geom_vline(xintercept = input$x_threshold_max, color = "red", linewidth = 1)
      }
      
      if (y_op %in% c("gt", "gte")) {
        p <- p + geom_hline(yintercept = input$y_threshold, color = "red",
                            linewidth = 1, linetype = if (y_op == "gt") "dotted" else "solid")
      } else if (y_op %in% c("lt", "lte")) {
        p <- p + geom_hline(yintercept = input$y_threshold, color = "red",
                            linewidth = 1, linetype = if (y_op == "lt") "dotted" else "solid")
      } else if (y_op == "between") {
        p <- p + geom_hline(yintercept = input$y_threshold, color = "red", linewidth = 1) +
                 geom_hline(yintercept = input$y_threshold_max, color = "red", linewidth = 1)
      }
      
    } else if (!is.null(rv$saved_brush) && isTRUE(input$numeric_gate_yn == "no")) {
      brush <- rv$saved_brush
      p <- p + annotate("rect",
                        xmin = brush$xmin, xmax = brush$xmax,
                        ymin = brush$ymin, ymax = brush$ymax,
                        fill = NA, color = "red",
                        linewidth = 1.5, linetype = "dashed")
    }
    
    p
  })

  ############################################
  # Gate statistics
  ############################################
  output$gate_stats <- renderText({
    req(rv$current_data)
    
    if (rv$current_gate_id == "root") {
      sprintf("Viewing root population\nTotal cells: %d", nrow(rv$current_data))
    } else {
      path <- get_gate_path(rv$gate_tree, rv$current_gate_id)
      path_str <- paste(c("Root", path), collapse = " \u2192 ")
      sprintf("Current path: %s\nCells at this gate: %d\nTotal gates in tree: %d",
              path_str, nrow(rv$current_data), length(rv$gate_tree))
    }
  })

  output$download_gate_metadata <- downloadHandler(
    filename = function() paste0("gate_metadata_", Sys.Date(), ".json"),
    content = function(file) {
      clean_tree <- lapply(rv$gate_tree, function(g) { g$plot <- NULL; g })
      writeLines(jsonlite::toJSON(clean_tree, pretty = TRUE, auto_unbox = TRUE), file)
    }
  )

  ############################################
  # Helper functions
  ############################################
  get_gate_path <- function(gate_tree, gate_id) {
    if (gate_id == "root" || !gate_id %in% names(gate_tree)) return(c())
    gate_info <- gate_tree[[gate_id]]
    c(get_gate_path(gate_tree, gate_info$parent), gate_id)
  }
  
  get_all_descendants <- function(gate_tree, gate_id) {
    if (gate_id == "root") return(names(gate_tree))
    children <- get_children(gate_tree, gate_id)
    all_desc <- children
    for (child in children) {
      all_desc <- c(all_desc, get_all_descendants(gate_tree, child))
    }
    all_desc
  }
  
  gate_depth <- function(gate_tree, gate_id) {
    if (gate_id == "root" || gate_id == "Ungated") return(0)
    if (!gate_id %in% names(gate_tree)) return(0)
    1 + gate_depth(gate_tree, gate_tree[[gate_id]]$parent)
  }
  
  ############################################
  # All gates summary
  ############################################
  output$all_gates_summary <- renderText({
    if (length(rv$gate_tree) == 0) return("No gates applied yet.")
    
    summary_lines <- c("=== Complete Gating Strategy ===\n")
    root_gates <- names(rv$gate_tree)[sapply(rv$gate_tree, function(g) g$parent == "root")]
    for (root_gate in root_gates) {
      summary_lines <- c(summary_lines, format_gate_summary(rv$gate_tree, root_gate, 0))
    }
    paste(summary_lines, collapse = "\n")
  })
  
  format_gate_summary <- function(gate_tree, gate_id, level) {
    gate_info <- gate_tree[[gate_id]]
    indent <- paste(rep("  ", level), collapse = "")
    
    path <- get_gate_path(gate_tree, gate_id)
    path_str <- paste(c("Root", path), collapse = " \u2192 ")
    
    total_cells <- nrow(rv$original_data)
    pct_of_total <- 100 * gate_info$n_after / total_cells
    
    lines <- c()
    lines <- c(lines, sprintf("%s%s:", indent, gate_id))
    lines <- c(lines, sprintf("%s  Path: %s", indent, path_str))
    lines <- c(lines, sprintf("%s  Variables: %s vs %s", indent, 
                              get_transform_label(gate_info$x_var, gate_info$x_transform %||% "none"),
                              get_transform_label(gate_info$y_var, gate_info$y_transform %||% "none")))
    lines <- c(lines, sprintf("%s  Thresholds: X[%.2f, %.2f], Y[%.2f, %.2f]",
                              indent, 
                              gate_info$threshold_xmin, gate_info$threshold_xmax,
                              gate_info$threshold_ymin, gate_info$threshold_ymax))
    lines <- c(lines, sprintf("%s  Cells: %d \u2192 %d (%.1f%% of parent, %.1f%% of total)", 
                              indent, gate_info$n_before, gate_info$n_after,
                              100 * gate_info$n_after / gate_info$n_before,
                              pct_of_total))
    
    children <- get_children(gate_tree, gate_id)
    if (length(children) > 0) {
      lines <- c(lines, sprintf("%s  Children: %s", indent, paste(children, collapse = ", ")))
      for (child in children) {
        lines <- c(lines, format_gate_summary(gate_tree, child, level + 1))
      }
    }
    
    c(lines, "")
  }

  ############################################
  # Transformation functions
  ############################################
  apply_transformation <- function(x, transform_type) {
    switch(transform_type,
           "none" = x,
           "log2" = {
             min_val <- min(x, na.rm = TRUE)
             offset <- ifelse(min_val <= 0, abs(min_val) + 1, 0)
             log2(x + offset)
           },
           "ln" = {
             min_val <- min(x, na.rm = TRUE)
             offset <- ifelse(min_val <= 0, abs(min_val) + 1, 0)
             log(x + offset)
           },
           "log10" = {
             min_val <- min(x, na.rm = TRUE)
             offset <- ifelse(min_val <= 0, abs(min_val) + 1, 0)
             log10(x + offset)
           },
           "zscore" = scale(x)[,1],
           x
    )
  }

  get_transform_label <- function(var_name, transform_type) {
    switch(transform_type,
           "none" = var_name,
           "log2" = paste0("log2(", var_name, ")"),
           "ln" = paste0("ln(", var_name, ")"),
           "log10" = paste0("log10(", var_name, ")"),
           "zscore" = paste0("Z-score(", var_name, ")"),
           var_name
    )
  }

  ############################################
  # Reactive: transformed data
  ############################################
  transformed_current_data <- reactive({
    req(rv$current_data, input$x_axis, input$y_axis)
    data <- rv$current_data
    req(is.numeric(data[[input$x_axis]]), is.numeric(data[[input$y_axis]]))
    
    x_scale <- input$x_scale %||% "none"
    y_scale <- input$y_scale %||% "none"
    
    data$x_transformed <- apply_transformation(data[[input$x_axis]], x_scale)
    data$y_transformed <- apply_transformation(data[[input$y_axis]], y_scale)
    data
  })

  x_data_range <- reactive({
    req(transformed_current_data())
    range(transformed_current_data()$x_transformed, na.rm = TRUE)
  })
  
  y_data_range <- reactive({
    req(transformed_current_data())
    range(transformed_current_data()$y_transformed, na.rm = TRUE)
  })

  ############################################
  # Sliders with 3 decimal precision
  ############################################
  output$x_slider_ui <- renderUI({
    req(x_data_range())
    range_vals <- x_data_range()
    range_min <- floor(range_vals[1] * 1000) / 1000
    range_max <- ceiling(range_vals[2] * 1000) / 1000
    step_val <- max(0.001, round((range_max - range_min) / 1000, 3))
    
    sliderInput("x_range_slider", NULL,
                min = range_min, max = range_max,
                value = c(range_min, range_max),
                step = step_val)
  })

  output$y_slider_ui <- renderUI({
    req(y_data_range())
    range_vals <- y_data_range()
    range_min <- floor(range_vals[1] * 1000) / 1000
    range_max <- ceiling(range_vals[2] * 1000) / 1000
    step_val <- max(0.001, round((range_max - range_min) / 1000, 3))
    
    sliderInput("y_range_slider", NULL,
                min = range_min, max = range_max,
                value = c(range_min, range_max),
                step = step_val)
  })

  ############################################
  # Histograms
  ############################################
  output$x_histogram <- renderPlot({
    req(transformed_current_data(), input$x_axis)
    data <- transformed_current_data()
    
    x_slider <- input$x_range_slider
    if (is.null(x_slider)) x_slider <- x_data_range()
    
    p <- ggplot(data, aes(x = x_transformed)) +
      geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
      geom_vline(xintercept = x_slider[1], color = "red", linewidth = 1.5) +
      geom_vline(xintercept = x_slider[2], color = "red", linewidth = 1.5) +
      labs(x = NULL, y = NULL) +
      theme_minimal() +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            panel.grid = element_blank(),
            plot.margin = margin(2, 5, 2, 5))
    print(p)
  }, height = 80)

  output$y_histogram <- renderPlot({
    req(transformed_current_data(), input$y_axis)
    data <- transformed_current_data()
    
    y_slider <- input$y_range_slider
    if (is.null(y_slider)) y_slider <- y_data_range()
    
    p <- ggplot(data, aes(x = y_transformed)) +
      geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
      geom_vline(xintercept = y_slider[1], color = "red", linewidth = 1.5) +
      geom_vline(xintercept = y_slider[2], color = "red", linewidth = 1.5) +
      labs(x = NULL, y = NULL) +
      theme_minimal() +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            panel.grid = element_blank(),
            plot.margin = margin(2, 5, 2, 5))
    print(p)
  }, height = 80)

  ############################################
  # Final data table
  ############################################
  output$final_data <- renderDT({
    req(rv$current_data)
    datatable(rv$current_data, options = list(scrollX = TRUE))
  })

  ############################################
  # Downloads
  ############################################
  output$download_plots <- downloadHandler(
    filename = function() paste0("gate_plots_", Sys.Date(), ".zip"),
    content = function(file) {
      req(length(rv$gate_tree) > 0)
      temp_dir <- tempdir()
      plot_dir <- file.path(temp_dir, "gate_plots")
      dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
      
      for (gate_name in names(rv$gate_tree)) {
        gate_info <- rv$gate_tree[[gate_name]]
        if (!is.null(gate_info$plot)) {
          pdf_file <- file.path(plot_dir, paste0(gate_name, ".pdf"))
          ggsave(pdf_file, plot = gate_info$plot, width = 8, height = 6, device = "pdf")
        }
      }
      zip::zip(file, files = list.files(plot_dir, full.names = TRUE), mode = "cherry-pick")
    }
  )
  
  output$download_csv <- downloadHandler(
    filename = function() paste0("gated_data_final_", Sys.Date(), ".csv"),
    content = function(file) {
      req(rv$full_data_with_gates)
      final_rows <- rownames(rv$current_data)
      output_data <- rv$full_data_with_gates[final_rows, ]
      write.csv(output_data, file, row.names = FALSE)
    }
  )
  
  output$download_all <- downloadHandler(
    filename = function() paste0("gating_complete_package_", Sys.Date(), ".zip"),
    content = function(file) {
      req(length(rv$gate_tree) > 0)
      req(rv$full_data_with_gates)
      
      temp_dir <- tempdir()
      package_dir <- file.path(temp_dir, "gating_package")
      dir.create(package_dir, showWarnings = FALSE, recursive = TRUE)
      
      plot_dir <- file.path(package_dir, "plots")
      dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
      
      for (gate_name in names(rv$gate_tree)) {
        gate_info <- rv$gate_tree[[gate_name]]
        if (!is.null(gate_info$plot)) {
          pdf_file <- file.path(plot_dir, paste0(gate_name, ".pdf"))
          ggsave(pdf_file, plot = gate_info$plot, width = 8, height = 6, device = "pdf")
        }
      }
      
      final_rows <- rownames(rv$current_data)
      final_data <- rv$full_data_with_gates[final_rows, ]
      write.csv(final_data, file.path(package_dir, "gated_data_final.csv"), row.names = FALSE)
      write.csv(rv$full_data_with_gates, file.path(package_dir, "original_data_with_gates.csv"), row.names = FALSE)
      
      # Generate gating strategy summary
      summary_lines <- c("=== Complete Gating Strategy ===", "")
      summary_lines <- c(summary_lines, paste("Generated:", Sys.time()))
      summary_lines <- c(summary_lines, paste("Total gates:", length(rv$gate_tree)))
      summary_lines <- c(summary_lines, paste("Original cells:", nrow(rv$original_data)))
      summary_lines <- c(summary_lines, "")
      
      root_gates <- names(rv$gate_tree)[sapply(rv$gate_tree, function(g) g$parent == "root")]
      for (root_gate in root_gates) {
        summary_lines <- c(summary_lines, format_gate_summary_text(rv$gate_tree, root_gate, 0))
      }
      writeLines(summary_lines, file.path(package_dir, "gating_strategy_summary.txt"))
      
      clean_tree <- lapply(rv$gate_tree, function(g) { g$plot <- NULL; g })
      metadata <- jsonlite::toJSON(clean_tree, pretty = TRUE, auto_unbox = TRUE)
      writeLines(metadata, file.path(package_dir, "gate_metadata.json"))
      
      zip::zip(file, files = list.files(package_dir, recursive = TRUE, full.names = TRUE),
               mode = "cherry-pick")
    }
  )
  
  output$download_current_view <- downloadHandler(
    filename = function() {
      view_name <- if (rv$current_gate_id == "root") "root_view" else rv$current_gate_id
      paste0(view_name, "_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      req(rv$current_data, input$x_axis, input$y_axis)
      
      plot_data <- rv$current_data
      current_row_names <- rownames(plot_data)
      
      all_children <- get_all_descendants(rv$gate_tree, rv$current_gate_id)
      visible_children <- intersect(all_children, rv$visible_gates)
      
      if (length(visible_children) > 0) {
        gate_member_matrix <- matrix(FALSE,
                                      nrow = nrow(plot_data),
                                      ncol = length(visible_children),
                                      dimnames = list(current_row_names, visible_children))
        
        for (gate_id in visible_children) {
          matching_col <- find_gate_column(gate_id, names(rv$full_data_with_gates))
          if (!is.null(matching_col)) {
            passing_rows <- rownames(rv$full_data_with_gates)[rv$full_data_with_gates[[matching_col]] == 1]
            gate_member_matrix[current_row_names %in% passing_rows, gate_id] <- TRUE
          }
        }
        
        gate_depths <- sapply(visible_children, function(g) gate_depth(rv$gate_tree, g))
        gate_parents <- sapply(visible_children, function(g) rv$gate_tree[[g]]$parent)
        
        gate_memberships <- character(nrow(plot_data))
        gate_memberships[] <- "Ungated"
        
        for (i in seq_len(nrow(plot_data))) {
          member_gates <- visible_children[gate_member_matrix[i, ]]
          
          if (length(member_gates) > 0) {
            depths <- gate_depths[member_gates]
            max_depth <- max(depths)
            deepest <- member_gates[depths == max_depth]
            
            if (length(deepest) > 1) {
              parents <- gate_parents[deepest]
              siblings <- deepest[parents == parents[1]]
              if (length(siblings) > 1) {
                gate_memberships[i] <- paste(sort(siblings), collapse = "+")
              } else {
                gate_memberships[i] <- deepest[1]
              }
            } else {
              gate_memberships[i] <- deepest[1]
            }
          }
        }
        
        plot_data$gate_membership <- gate_memberships
      } else {
        plot_data$gate_membership <- "Ungated"
      }
      
      # Build color map
      unique_classes <- unique(plot_data$gate_membership)
      gate_colors_map <- c("Ungated" = "gray80")
      
      for (gate_id in visible_children) {
        if (gate_id %in% names(rv$gate_colors)) {
          gate_colors_map[gate_id] <- rv$gate_colors[gate_id]
        }
      }
      
      overlap_classes <- unique_classes[grepl("\\+", unique_classes)]
      for (ov in overlap_classes) {
        component_gates <- unlist(strsplit(ov, "\\+"))
        parent_colors <- rv$gate_colors[component_gates]
        gate_colors_map[ov] <- blend_colors(parent_colors)
      }
      
      class_counts <- table(plot_data$gate_membership)
      
      if (isTRUE(input$show_counts_in_legend)) {
        legend_labels <- names(gate_colors_map)
        legend_labels_with_counts <- paste0(
          legend_labels, " (", class_counts[legend_labels] %||% 0, ")")
        names(gate_colors_map) <- legend_labels_with_counts
        
        plot_data$display_color <- paste0(
          plot_data$gate_membership, " (", class_counts[plot_data$gate_membership], ")")
        plot_data$display_color <- factor(plot_data$display_color, levels = legend_labels_with_counts)
      } else {
        plot_data$display_color <- factor(plot_data$gate_membership, levels = names(gate_colors_map))
      }
      
      # Draw ungated first (behind), gated on top
      plot_data <- plot_data[order(plot_data$gate_membership == "Ungated", decreasing = TRUE), ]
      
      point_size <- max(0.3, min(4, 6 / log10(nrow(plot_data) + 1)))
      
      p <- ggplot(plot_data, aes(x = .data[[input$x_axis]], y = .data[[input$y_axis]])) +
        geom_point(aes(color = display_color), alpha = 1, size = point_size) +
        scale_color_manual(values = gate_colors_map, name = "Gate", drop = FALSE) +
        labs(title = paste("View:", rv$current_gate_id),
             subtitle = sprintf("%s vs %s", input$x_axis, input$y_axis)) +
        theme_bw() +
        theme(text = element_text(size = 14),
              legend.position = "right",
              legend.background = element_rect(fill = "white", color = "black"),
              legend.key.size = unit(0.8, "lines"))
      
      ggsave(file, plot = p, width = 10, height = 8, device = "pdf")
    }
  )

  ############################################
  # Text file summary helper
  ############################################
  format_gate_summary_text <- function(gate_tree, gate_id, level) {
    gate_info <- gate_tree[[gate_id]]
    indent <- paste(rep("  ", level), collapse = "")
    
    path <- get_gate_path(gate_tree, gate_id)
    path_str <- paste(c("Root", path), collapse = " -> ")
    
    total_cells <- nrow(rv$original_data)
    pct_of_total <- 100 * gate_info$n_after / total_cells
    
    lines <- c()
    lines <- c(lines, sprintf("%s%s:", indent, gate_id))
    lines <- c(lines, sprintf("%s  Path: %s", indent, path_str))
    lines <- c(lines, sprintf("%s  Variables: %s vs %s", indent, 
                              get_transform_label(gate_info$x_var, gate_info$x_transform %||% "none"),
                              get_transform_label(gate_info$y_var, gate_info$y_transform %||% "none")))
    lines <- c(lines, sprintf("%s  Cells: %d -> %d (%.1f%% of parent, %.1f%% of total)", 
                              indent, gate_info$n_before, gate_info$n_after,
                              100 * gate_info$n_after / gate_info$n_before,
                              pct_of_total))
    lines <- c(lines, sprintf("%s  Thresholds: X[%.2f, %.2f], Y[%.2f, %.2f]",
                              indent,
                              gate_info$threshold_xmin, gate_info$threshold_xmax,
                              gate_info$threshold_ymin, gate_info$threshold_ymax))
    
    children <- get_children(gate_tree, gate_id)
    if (length(children) > 0) {
      lines <- c(lines, sprintf("%s  Children: %s", indent, paste(children, collapse = ", ")))
      for (child in children) {
        lines <- c(lines, format_gate_summary_text(gate_tree, child, level + 1))
      }
    }
    c(lines, "")
  }
  
  ############################################
  # Helper: create gate snapshot plot
  ############################################
  create_gate_plot <- function(data, x_var, y_var, brush, gate_name, inside, x_label = x_var, y_label = y_var) {
    data$selected <- inside
    point_size <- max(0.3, min(4, 6 / log10(nrow(data) + 1)))
    
    p <- ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]])) +
      geom_point(aes(color = selected), alpha = 1, size = point_size) +
      scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"),
                         labels = c("Excluded", "Included"),
                         name = "Gate Status") +
      annotate("rect",
               xmin = brush$xmin, xmax = brush$xmax,
               ymin = brush$ymin, ymax = brush$ymax,
               fill = NA, color = "red", linewidth = 1.5, linetype = "dashed") +
      labs(title = gate_name,
           subtitle = sprintf("%s vs %s\n%d cells gated", x_label, y_label, sum(inside)),
           x = x_label, y = y_label) +
      theme_bw() +
      theme(text = element_text(size = 14), legend.position = "bottom")
    
    p
  }
})
