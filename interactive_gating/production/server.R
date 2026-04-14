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
source('/srv/shiny-server/phenomenalist/utils/provenance.R')

`%||%` <- function(a, b) if (!is.null(a)) a else b

# Helper function to blend colors for overlaps
blend_colors <- function(colors) {
  if (length(colors) == 1) return(colors)
  
  # Convert hex to RGB
  rgb_vals <- col2rgb(colors)
  # Average the RGB values
  avg_rgb <- rowMeans(rgb_vals)
  # Convert back to hex
  rgb(avg_rgb[1], avg_rgb[2], avg_rgb[3], maxColorValue = 255)
}
find_gate_column <- function(gate_name, col_names) {
  if (gate_name %in% col_names) return(gate_name)
  safe_name <- make.names(gate_name)
  if (safe_name %in% col_names) return(safe_name)
  return(NULL)
}
shinyServer(function(input, output, session) {
  
  # Reactive values to store data and gates
  rv <- reactiveValues(
    original_data = NULL,
    current_data = NULL,
    full_data_with_gates = NULL,
    gate_number = 0,
    gate_tree = list(),  # Hierarchical gate structure
    current_gate_id = NULL,  # Currently selected gate
    saved_brush = NULL,
    gate_data_cache = list(),  # Store data at each gate for navigation
    gate_colors = c(),  # Store colors for each gate
    visible_gates = c(),  # Gates to display on plot
    x_range = NULL,
    y_range = NULL
  )
  prov_dir <- file.path(tempdir(), paste0("interactive_gating_", session$token))
  dir.create(prov_dir, showWarnings = FALSE)
  tracker <- ProvenanceTracker$new("interactive_gating", session, prov_dir)

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
  output$tmp_probe <- renderPrint({
  	td <- Sys.getenv("TMPDIR")
  	p  <- file.path(td, paste0("probe_", Sys.getpid()))
  	f  <- file.path(p, "x.txt")

  	list(
    	pid = Sys.getpid(),
    	TMPDIR = td,
    	tempdir = tempdir(),
    	mkdir = tryCatch({ dir.create(p); TRUE }, error=function(e) e$message),
    	write = tryCatch({ writeLines("ok", f); TRUE }, error=function(e) e$message),
    	cleanup = tryCatch({ unlink(p, recursive=TRUE); TRUE }, error=function(e) e$message)
  	)
  })
  # shows current env + tempdir live
  output$tmp_env <- renderPrint({
  	list(
    	pid = Sys.getpid(),
    	user = Sys.info()[["user"]],
    	tempdir = tempdir(),
    	env = Sys.getenv(c("TMPDIR","TMP","TEMP","R_SESSION_TMPDIR"))
  	)
  })

 # run probe on demand
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
      # Use fread which handles .csv.gz automatically
      data <- data.table::fread(input$csv_file$datapath, stringsAsFactors = FALSE)
      
      # Convert to regular data frame to ensure proper rowname handling
      data <- as.data.frame(data)
      
      # CRITICAL: Explicitly set rownames to ensure they persist through filtering
      rownames(data) <- as.character(seq_len(nrow(data)))
      
      # Detect gate columns (looking for " Gate" suffix)
      gate_cols <- names(data)[grepl("[. ]Gate$", names(data))]
      
      rv$detected_gates <- gate_cols
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
      
      # Update gate name for first gate
      updateTextInput(session, "gate_name", value = "Gate_1")
      
      showNotification("Data loaded successfully!", type = "message")
    }, error = function(e) {
      showNotification(paste("Error loading data:", e$message), type = "error")
    })
  })
  
  # Data preview (conditional UI)
  output$data_preview_box <- renderUI({
    req(rv$original_data)
    box(
      title = "Data Preview",
      width = 12,
      DTOutput("data_preview")
    )
  })
  
  output$data_preview <- renderDT({
    req(rv$original_data)
    datatable(head(rv$original_data, 100), options = list(scrollX = TRUE))
  })
  
    
    
  # Dynamic UI for axis selection
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
  output$previous_gates_detected <- renderUI({
  	req(rv$detected_gates)
  	if (length(rv$detected_gates) > 0) {
    		box(
      			title = "Previous Gates Detected!",
      			status = "warning",
      			solidHeader = TRUE,
      			width = 12,
      			p(sprintf("Found %d gate columns in your data.", length(rv$detected_gates))),
      			p(strong("To restore your previous gating session with exact thresholds:")),
      			fileInput("gate_metadata_file", 
                	"Upload gate_metadata.json file:",
                	accept = c("application/json", ".json")),
      			actionButton("restore_gates", 
                   	"Restore Gates from Metadata", 
                   	class = "btn-success"),
      			hr(),
      			p(em("Note: Without the metadata file, gate definitions (thresholds, axes) cannot be restored."))
    		)
  	}
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
    
    # Prevent infinite loops - only update if actually changing gates
    if (!is.null(gate_id) && gate_id != rv$current_gate_id && gate_id %in% names(rv$gate_data_cache)) {
      rv$current_gate_id <- gate_id
      rv$current_data <- rv$gate_data_cache[[gate_id]]
      rv$saved_brush <- NULL
      session$resetBrush("plot_brush")
      
      # Update gate name for next child
      children <- get_children(rv$gate_tree, gate_id)
      next_num <- length(children) + 1
      if (gate_id == "root") {
        updateTextInput(session, "gate_name", value = paste0("Gate_", next_num))
      } else {
        updateTextInput(session, "gate_name", value = paste0(gate_id, ".", next_num))
      }
    }
  }, ignoreInit = TRUE)
  
  #observeEvent(input$restore_gates, {
  #	req(rv$original_data, input$gate_metadata_file)
  #
  #	tryCatch({
   # 		# Read JSON metadata
   # 		metadata <- jsonlite::fromJSON(input$gate_metadata_file$datapath)
    
    #		# Restore gate tree directly
    #		rv$gate_tree <- metadata
    #
    #		# Rebuild gate data cache
    #		for (gate_name in names(metadata)) {
     # 			gate_info <- metadata[[gate_name]]
      #
      #			# Get cells that passed this gate
      #			passing_rows <- which(rv$full_data_with_gates[[gate_name]] == 1)
      #			rv$gate_data_cache[[gate_name]] <- rv$original_data[passing_rows, ]
      #
      #			# Restore color
      #			rv$gate_number <- rv$gate_number + 1
      #			color_palette <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")
      #			rv$gate_colors[gate_name] <- color_palette[(rv$gate_number - 1) %% 5 + 1]
      #
      #			rv$visible_gates <- c(rv$visible_gates, gate_name)
    	#	}
    #
    #		showNotification(sprintf("Restored %d gates with exact thresholds!", 
     #                       length(metadata)), type = "message")
    #
  #	}, error = function(e) {
   # 	showNotification(paste("Error restoring gates:", e$message), type = "error")
  #	})
  #})
  observeEvent(input$restore_gates, {
  	req(rv$original_data, input$gate_metadata_file)
  
  	tryCatch({
    	# Read JSON metadata
    	metadata <- jsonlite::fromJSON(input$gate_metadata_file$datapath)
    
    	# Restore gate tree directly
    	rv$gate_tree <- metadata
    
    	# Reset gate data cache with root
    	rv$gate_data_cache <- list(root = rv$original_data)
    
    	# Rebuild gates hierarchically by re-applying thresholds
    	for (gate_name in names(metadata)) {
      		gate_info <- metadata[[gate_name]]
      
      		# Get parent data
      		parent_id <- gate_info$parent
      		if (!(parent_id %in% names(rv$gate_data_cache))) {
        		stop(paste("Parent gate", parent_id, "not found for", gate_name))
      		}
      		parent_data <- rv$gate_data_cache[[parent_id]]
      
      		# RE-APPLY the gate using saved thresholds
      		x_vals <- parent_data[[gate_info$x_var]]
      		y_vals <- parent_data[[gate_info$y_var]]
      
      		inside <- x_vals >= gate_info$threshold_xmin & 
                x_vals <= gate_info$threshold_xmax &
                y_vals >= gate_info$threshold_ymin & 
                y_vals <= gate_info$threshold_ymax
      
      		# Cache the gated data
      		gated_data <- parent_data[inside, ]
      		rv$gate_data_cache[[gate_name]] <- gated_data
      
      		# Update full dataset gate column
      		if (!(gate_name %in% names(rv$full_data_with_gates))) {
        		rv$full_data_with_gates[[gate_name]] <- 0
      		}
      		rows_that_pass <- rownames(parent_data)[inside]
      		rv$full_data_with_gates[rows_that_pass, gate_name] <- 1
      
      		# Restore color
      		rv$gate_number <- rv$gate_number + 1
      		color_palette <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")
      		rv$gate_colors[gate_name] <- color_palette[(rv$gate_number - 1) %% 5 + 1]
      
      		rv$visible_gates <- c(rv$visible_gates, gate_name)
       }
    
       # Set current view to root after restoration
       rv$current_gate_id <- "root"
       rv$current_data <- rv$original_data
    
       showNotification(sprintf("Restored %d gates with exact thresholds!", 
                        length(metadata)), type = "message")
    
  	}, error = function(e) {
    	showNotification(paste("Error restoring gates:", e$message), type = "error")
  	})
   })
  # Apply gate with lasso/brush selection
  # Update the observeEvent(input$apply_gate, {...}) section
  observeEvent(input$apply_gate, {
  	req(rv$current_data, input$x_axis, input$y_axis, rv$saved_brush)
  	req(input$gate_name != "")
  	req(transformed_current_data())
  
  	brush <- rv$saved_brush
  
  	# Get TRANSFORMED coordinates
  	data <- transformed_current_data()
  	x_vals <- data$x_transformed
  	y_vals <- data$y_transformed
  
  	# Apply selection in transformed space
  	inside <- x_vals >= brush$xmin & x_vals <= brush$xmax &
            y_vals >= brush$ymin & y_vals <= brush$ymax
  
  	# Rest of the gate application logic remains the same...
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
  
  	# Store transformation info in gate metadata
  	gate_info <- list(
    	id = gate_name,
    	parent = rv$current_gate_id,
    	x_var = input$x_axis,
    	y_var = input$y_axis,
    	x_transform = input$x_scale %||% "none",
    	y_transform = input$y_scale %||% "none",
    	n_before = nrow(rv$current_data),
    	n_after = sum(inside),
    	brush = brush,
    	children = list(),
    	threshold_xmin = brush$xmin,
    	threshold_xmax = brush$xmax,
    	threshold_ymin = brush$ymin,
    	threshold_ymax = brush$ymax
  	)
  
  	rv$gate_tree[[gate_name]] <- gate_info
  
  	matching_col <- find_gate_column(gate_name, names(rv$full_data_with_gates))
  	if (is.null(matching_col)) {
    		matching_col <- make.names(gate_name)
  	}
  
  	if (!(matching_col %in% names(rv$full_data_with_gates))) {
    		rv$full_data_with_gates[[matching_col]] <- 0
  	}
  
  	rows_that_pass <- rownames(rv$current_data)[inside]
  	rv$full_data_with_gates[rows_that_pass, matching_col] <- 1
  
  	gated_data <- rv$current_data[inside, ]
  	rv$gate_data_cache[[gate_name]] <- gated_data
  
  	temp_plot <- create_gate_plot(data, "x_transformed", "y_transformed", 
                                 brush, gate_name, inside,
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
  # Reset current brush/lasso
  observeEvent(input$reset_gate, {
    rv$saved_brush <- NULL
    session$resetBrush("plot_brush")
    showNotification("Lasso selection cleared.", type = "message")
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
  
  # Render gate tree with radio buttons
  output$gate_tree <- renderUI({
    req(rv$original_data)
    
    # Force a dependency on gate_tree to update when gates change
    gate_count <- length(rv$gate_tree)
    
    if (gate_count == 0) {
      return(div(
        radioButtons("select_gate", NULL, 
                     choices = c("Root (All cells)" = "root"),
                     selected = "root")
      ))
    }
    
    # Build hierarchical display
    tree_html <- build_gate_tree_ui(rv$gate_tree, rv$current_gate_id)
    
    div(
      radioButtons("select_gate", NULL, 
                   choices = tree_html$choices,
                   selected = rv$current_gate_id)
    )
  })
  
  # Helper function to build tree UI
  build_gate_tree_ui <- function(gate_tree, current_id) {
    choices <- c("Root (All cells)" = "root")
    html_parts <- c()
    
    # Find root gates (gates with parent = "root")
    root_gates <- names(gate_tree)[sapply(gate_tree, function(g) g$parent == "root")]
    
    for (gate_id in root_gates) {
      gate_info <- gate_tree[[gate_id]]
      label <- sprintf("%s (%d→%d, %.1f%%)", 
                       gate_id, 
                       gate_info$n_before, 
                       gate_info$n_after,
                       100 * gate_info$n_after / gate_info$n_before)
      choices[label] <- gate_id
      
      # Add children recursively
      children <- get_children(gate_tree, gate_id)
      if (length(children) > 0) {
        child_info <- build_children_ui(gate_tree, gate_id, 1)
        choices <- c(choices, child_info$choices)
      }
    }
    
    list(choices = choices, html = "")
  }
  
  # Helper to build children display
  build_children_ui <- function(gate_tree, parent_id, level) {
    children <- get_children(gate_tree, parent_id)
    choices <- c()
    
    for (gate_id in children) {
      gate_info <- gate_tree[[gate_id]]
      indent <- paste(rep("  ", level), collapse = "")
      label <- sprintf("%s↳ %s (%d→%d, %.1f%%)", 
                       indent,
                       gate_id, 
                       gate_info$n_before, 
                       gate_info$n_after,
                       100 * gate_info$n_after / gate_info$n_before)
      choices[label] <- gate_id
      
      # Recursively add grandchildren
      grandchildren <- get_children(gate_tree, gate_id)
      if (length(grandchildren) > 0) {
        child_info <- build_children_ui(gate_tree, gate_id, level + 1)
        choices <- c(choices, child_info$choices)
      }
    }
    
    list(choices = choices)
  }
  
  # Helper to get children of a gate
  get_children <- function(gate_tree, parent_id) {
    children <- names(gate_tree)[sapply(gate_tree, function(g) g$parent == parent_id)]
    return(children)
  }
  
  # Render gate display checkboxes
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
  
  # Update visible gates based on checkbox selection
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
  
  # Create plot
  # Replace the existing output$gate_plot with this updated version
  output$gate_plot <- renderPlot({
  	req(rv$current_data, input$x_axis, input$y_axis)
  	req(transformed_current_data())
  
  	# Get transformation settings
  	x_scale <- input$x_scale %||% "none"
  	y_scale <- input$y_scale %||% "none"
  
  	# Get axis labels with transformation
  	x_label <- get_transform_label(input$x_axis, x_scale)
  	y_label <- get_transform_label(input$y_axis, y_scale)
  
  	# Work on transformed data
  	plot_data <- transformed_current_data()
  
  	# Determine which points are inside saved brush (in TRANSFORMED space)
  	inside_brush <- rep(FALSE, nrow(plot_data))
  	if (!is.null(rv$saved_brush)) {
    		brush <- rv$saved_brush
    		inside_brush <- plot_data$x_transformed >= brush$xmin & 
                    plot_data$x_transformed <= brush$xmax &
                    plot_data$y_transformed >= brush$ymin & 
                    plot_data$y_transformed <= brush$ymax
  	}
  
  	plot_data$selected <- inside_brush
  
  	# Gate membership logic (same as before)
  	current_row_names <- rownames(plot_data)
  	all_children <- get_all_descendants(rv$gate_tree, rv$current_gate_id)
  	visible_children <- intersect(all_children, rv$visible_gates)
  
  	if (length(visible_children) > 0) {
    		
   		# Pre-compute gate membership (once per gate, not per cell)
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
        				}else {
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
      			legend_labels,
      			" (",
      			class_counts[legend_labels] %||% 0,
      			")"
    			)
    		names(gate_colors_map) <- legend_labels_with_counts
    
    		plot_data$display_color <- ifelse(
      			plot_data$display_color == "Selected",
      			"Selected",
      			paste0(
        			plot_data$gate_membership,
        			" (",
        			class_counts[plot_data$gate_membership],
        			")"
      			)
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
  
  	# Performance: downsample for display
  	render_data <- plot_data
  	if (nrow(render_data) > MAX_DISPLAY_POINTS) {
    		set.seed(1)
    		render_data <- render_data[sample(nrow(render_data), MAX_DISPLAY_POINTS),]
  	}
  
  	use_raster <- nrow(render_data) > RASTER_THRESHOLD
        
	# Right before the ggplot call, after downsampling:
  	render_data <- render_data[order(render_data$gate_membership == "Ungated", decreasing = TRUE), ]
  	
	# Create plot with transformed axes
  	p <- ggplot(render_data,
              aes(x = x_transformed,
                  y = y_transformed))
  
 	# Dynamic point size based on cell count
  	point_size <- max(0.3, min(4, 6 / log10(nrow(render_data) + 1)))

  	if (use_raster) {
    		p <- p + geom_point_rast(aes(color = display_color),
                             alpha = 1,
                             size = point_size)
  	} else {
    		p <- p + geom_point(aes(color = display_color),
                        alpha = 1,
                        size = point_size)
  	}	 
  	p <- p + scale_color_manual(values = gate_colors_map,
                              name = "Gate",
                              drop = FALSE) +
    	labs(x = x_label, y = y_label) +
    	theme_bw() +
    	theme(text = element_text(size = 14),
          legend.position = "right",
          legend.background = element_rect(fill = "white", color = "black"),
          legend.key.size = unit(0.8, "lines"))
  
  	# Apply axis limits from sliders
  	if (!is.null(input$x_range_slider)) {
    		p <- p + xlim(input$x_range_slider)
  	}
  
  	if (!is.null(input$y_range_slider)) {
    		p <- p + ylim(input$y_range_slider)
  	}
  
  	# Add brush rectangle if present
  	if (!is.null(rv$saved_brush)) {
    		brush <- rv$saved_brush
    		p <- p + annotate("rect",
                      xmin = brush$xmin, xmax = brush$xmax,
                      ymin = brush$ymin, ymax = brush$ymax,
                      fill = NA, color = "red",
                      linewidth = 1.5, linetype = "dashed")
  	}
  
  	p
  }) 
  # Gate statistics
  output$gate_stats <- renderText({
    req(rv$current_data)
    
    if (rv$current_gate_id == "root") {
      sprintf("Viewing root population\nTotal cells: %d", nrow(rv$current_data))
    } else {
      # Build path from root to current gate
      path <- get_gate_path(rv$gate_tree, rv$current_gate_id)
      path_str <- paste(c("Root", path), collapse = " → ")
      
      sprintf("Current path: %s\nCells at this gate: %d\nTotal gates in tree: %d",
              path_str, nrow(rv$current_data), length(rv$gate_tree))
    }
  })
  output$download_gate_metadata <- downloadHandler(
  	filename = function() paste0("gate_metadata_", Sys.Date(), ".json"),
  	content = function(file) {
    		writeLines(jsonlite::toJSON(rv$gate_tree, pretty = TRUE), file)
  	}
  )
  # Helper to get path from root to gate
  get_gate_path <- function(gate_tree, gate_id) {
    if (gate_id == "root" || !gate_id %in% names(gate_tree)) {
      return(c())
    }
    
    gate_info <- gate_tree[[gate_id]]
    parent_path <- get_gate_path(gate_tree, gate_info$parent)
    return(c(parent_path, gate_id))
  }
  
  # Helper to get all descendants of a gate
  get_all_descendants <- function(gate_tree, gate_id) {
    if (gate_id == "root") {
      # Return all gates in order of creation
      return(names(gate_tree))
    }
    
    children <- get_children(gate_tree, gate_id)
    all_desc <- children
    
    for (child in children) {
      grandchildren <- get_all_descendants(gate_tree, child)
      all_desc <- c(all_desc, grandchildren)
    }
    
    return(all_desc)
  }
  
  # Helper to get gate depth
  gate_depth <- function(gate_tree, gate_id) {
    if (gate_id == "root" || gate_id == "Ungated") return(0)
    if (!gate_id %in% names(gate_tree)) return(0)
    
    gate_info <- gate_tree[[gate_id]]
    return(1 + gate_depth(gate_tree, gate_info$parent))
  }
  
  # All gates summary
  output$all_gates_summary <- renderText({
    if (length(rv$gate_tree) == 0) {
      return("No gates applied yet.")
    }
    
    summary_lines <- c("=== Complete Gating Strategy ===\n")
    
    # Get all gates in hierarchical order
    root_gates <- names(rv$gate_tree)[sapply(rv$gate_tree, function(g) g$parent == "root")]
    
    for (root_gate in root_gates) {
      summary_lines <- c(summary_lines, format_gate_summary(rv$gate_tree, root_gate, 0))
    }
    
    paste(summary_lines, collapse = "\n")
  })
  
  # Helper to format gate summary recursively
  format_gate_summary <- function(gate_tree, gate_id, level) {
    gate_info <- gate_tree[[gate_id]]
    indent <- paste(rep("  ", level), collapse = "")
    
    # Build path
    path <- get_gate_path(gate_tree, gate_id)
    path_str <- paste(c("Root", path), collapse = " → ")
    
    # Calculate percentage of total
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
    lines <- c(lines, sprintf("%s  Cells: %d → %d (%.1f%% of parent, %.1f%% of total)", 
                              indent, gate_info$n_before, gate_info$n_after,
                              100 * gate_info$n_after / gate_info$n_before,
                              pct_of_total))
    
    # Add children
    children <- get_children(gate_tree, gate_id)
    if (length(children) > 0) {
      lines <- c(lines, sprintf("%s  Children: %s", indent, paste(children, collapse = ", ")))
      for (child in children) {
        lines <- c(lines, format_gate_summary(gate_tree, child, level + 1))
      }
    }
    
    lines <- c(lines, "")  # Empty line between gates
    return(lines)
  }
  # Transformation functions
  apply_transformation <- function(x, transform_type) {
  	switch(transform_type,
         "none" = x,
         "log2" = {
           # Handle negative/zero values
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
         x  # default
  	)
  }

  # Get transformation label for axis
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
  # Final data table
  output$final_data <- renderDT({
    req(rv$current_data)
    datatable(rv$current_data, options = list(scrollX = TRUE))
  })
  # Reactive: Get transformed current data
  transformed_current_data <- reactive({
  	req(rv$current_data, input$x_axis, input$y_axis)
  
  	data <- rv$current_data
	req(is.numeric(data[[input$x_axis]]), is.numeric(data[[input$y_axis]]))

  	x_scale <- input$x_scale %||% "none"
  	y_scale <- input$y_scale %||% "none"
  
  	# Apply transformations
  	data$x_transformed <- apply_transformation(data[[input$x_axis]], x_scale)
  	data$y_transformed <- apply_transformation(data[[input$y_axis]], y_scale)
  
  	return(data)
   })

  # Reactive: Get X range based on transformation
  x_data_range <- reactive({
  	req(transformed_current_data())
  	data <- transformed_current_data()
  	range(data$x_transformed, na.rm = TRUE)
  })
  
  # Reactive: Get Y range based on transformation
  y_data_range <- reactive({
  	req(transformed_current_data())
  	data <- transformed_current_data()
  	range(data$y_transformed, na.rm = TRUE)
  })
  # Dynamic X slider
  output$x_slider_ui <- renderUI({
  	req(x_data_range())
  	range_vals <- x_data_range()
	range_min <- floor(range_vals[1] * 1000) / 1000
    	range_max <- ceiling(range_vals[2] * 1000) / 1000
    	step_val <- max(0.001, round((range_max - range_min) / 1000, 3))

  	sliderInput("x_range_slider",
              NULL,
              min = range_min,
              max = range_max,
              value = c(range_min,range_max),
              step = step_val)
   })

   # Dynamic Y slider
   output$y_slider_ui <- renderUI({
  	req(y_data_range())
  	
	range_vals <- y_data_range()
	range_min <- floor(range_vals[1] * 1000) / 1000
    	range_max <- ceiling(range_vals[2] * 1000) / 1000
    	step_val <- max(0.001, round((range_max - range_min) / 1000, 3))

	
  	sliderInput("y_range_slider",
              NULL,
              min = range_min,
              max = range_max,
              value = c(range_min,range_max),
              step = step_val)
   })
   # X-axis histogram
   output$x_histogram <- renderPlot({
  	req(transformed_current_data(), input$x_axis)
  
  	data <- transformed_current_data()
  	x_scale <- input$x_scale %||% "none"
  	x_label <- get_transform_label(input$x_axis, x_scale)
  
  	# Get current slider values
  	x_slider <- input$x_range_slider
  	if (is.null(x_slider)) {
    		x_slider <- x_data_range()
  	}
  
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

   # Y-axis histogram
   output$y_histogram <- renderPlot({
  	req(transformed_current_data(), input$y_axis)
  
  	data <- transformed_current_data()
  	y_scale <- input$y_scale %||% "none"
  	y_label <- get_transform_label(input$y_axis, y_scale)
  
  	# Get current slider values
  	y_slider <- input$y_range_slider
  	if (is.null(y_slider)) {
    		y_slider <- y_data_range()
  	}
  
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

    # Download plots as ZIP
    output$download_plots <- downloadHandler(
    	filename = function() {
      		paste0("gate_plots_", Sys.Date(), ".zip")
    	},
    	content = function(file) {
      		req(length(rv$gate_tree) > 0)
      
      		temp_dir <- tempdir()
      		plot_dir <- file.path(temp_dir, "gate_plots")
      		dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
      
      		for (gate_name in names(rv$gate_tree)) {
        		gate_info <- rv$gate_tree[[gate_name]]
        		if (!is.null(gate_info$plot)) {
          			pdf_file <- file.path(plot_dir, paste0(gate_name, ".pdf"))
          			ggsave(pdf_file, plot = gate_info$plot, 
                 			width = 8, height = 6, device = "pdf")
        		}
      		}
      
      		zip::zip(file, files = list.files(plot_dir, full.names = TRUE),
               		mode = "cherry-pick")
    	}
     )
  
  # Download final CSV
  output$download_csv <- downloadHandler(
    filename = function() {
      paste0("gated_data_final_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(rv$full_data_with_gates)
      
      final_rows <- rownames(rv$current_data)
      output_data <- rv$full_data_with_gates[final_rows, ]
      
      write.csv(output_data, file, row.names = FALSE)
    }
  )
  
  # Download complete package
  output$download_all <- downloadHandler(
    filename = function() {
      paste0("gating_complete_package_", Sys.Date(), ".zip")
    },
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
          ggsave(pdf_file, plot = gate_info$plot, 
                 width = 8, height = 6, device = "pdf")
        }
      }
      
      final_rows <- rownames(rv$current_data)
      final_data <- rv$full_data_with_gates[final_rows, ]
      write.csv(final_data, 
                file.path(package_dir, "gated_data_final.csv"), 
                row.names = FALSE)
      
      write.csv(rv$full_data_with_gates, 
                file.path(package_dir, "original_data_with_gates.csv"), 
                row.names = FALSE)
      
      # Generate gating strategy summary text file
      summary_lines <- c("=== Complete Gating Strategy ===", "")
      summary_lines <- c(summary_lines, paste("Generated:", Sys.time()))
      summary_lines <- c(summary_lines, paste("Total gates:", length(rv$gate_tree)))
      summary_lines <- c(summary_lines, paste("Original cells:", nrow(rv$original_data)))
      summary_lines <- c(summary_lines, "")
      
      # Get all gates in hierarchical order
      root_gates <- names(rv$gate_tree)[sapply(rv$gate_tree, function(g) g$parent == "root")]
      
      for (root_gate in root_gates) {
        summary_lines <- c(summary_lines, format_gate_summary_text(rv$gate_tree, root_gate, 0))
      }
      
      writeLines(summary_lines, file.path(package_dir, "gating_strategy_summary.txt"))
      clean_tree <- lapply(rv$gate_tree, function(g) {
      		g$plot <- NULL
  		g
      })
      metadata <- jsonlite::toJSON(clean_tree, pretty = TRUE, auto_unbox = TRUE)
      writeLines(metadata, file.path(package_dir, "gate_metadata.json"))   
      
      zip::zip(file, 
               files = list.files(package_dir, recursive = TRUE, full.names = TRUE),
               mode = "cherry-pick")
    }
  )
  
  output$download_current_view <- downloadHandler(
  	filename = function() {
    		view_name <- if(rv$current_gate_id == "root") "root_view" else rv$current_gate_id
    		paste0(view_name, "_", Sys.Date(), ".pdf")
  	},
  	content = function(file) {
    	req(rv$current_data, input$x_axis, input$y_axis)
    
    	# Create the current plot with SAME logic as main plot
    	plot_data <- rv$current_data
    	current_row_names <- rownames(plot_data)
    
    	all_children <- get_all_descendants(rv$gate_tree, rv$current_gate_id)
    	visible_children <- intersect(all_children, rv$visible_gates)
    
    	# REPLICATE THE OVERLAP LOGIC FROM output$gate_plot
    	if (length(visible_children) > 0) {
                # Pre-compute gate membership (once per gate, not per cell)
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
                                        }else { 
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
    
    	# Build color map with overlaps
    	unique_classes <- unique(plot_data$gate_membership)
    	gate_colors_map <- c("Ungated" = "gray80")
    
    	for (gate_id in visible_children) {
      		if (gate_id %in% names(rv$gate_colors)) {
        		gate_colors_map[gate_id] <- rv$gate_colors[gate_id]
      		}
    	}
    
    	# Add colors for overlap classes
    	overlap_classes <- unique_classes[grepl("\\+", unique_classes)]
    	for (ov in overlap_classes) {
      		component_gates <- unlist(strsplit(ov, "\\+"))
      		parent_colors <- rv$gate_colors[component_gates]
      		gate_colors_map[ov] <- blend_colors(parent_colors)
    	}
    
    	# Add counts to legend if enabled
    	class_counts <- table(plot_data$gate_membership)
    
    	if (isTRUE(input$show_counts_in_legend)) {
      		legend_labels <- names(gate_colors_map)
      		legend_labels_with_counts <- paste0(
        	legend_labels,
        		" (",
        	class_counts[legend_labels] %||% 0,
        	")"
      		)
      		names(gate_colors_map) <- legend_labels_with_counts
      
      		plot_data$display_color <- paste0(
        	plot_data$gate_membership,
        	" (",
        	class_counts[plot_data$gate_membership],
        	")"
      		)
      
      		plot_data$display_color <- factor(
        		plot_data$display_color,
        		levels = legend_labels_with_counts
      			)
    	} else {
      		plot_data$display_color <- factor(
        		plot_data$gate_membership,
        		levels = names(gate_colors_map)
      		)
    	}
        
	# Draw ungated first (behind), gated on top (visible)
  	plot_data <- plot_data[order(plot_data$gate_membership == "Ungated", decreasing = TRUE), ]

    	# Create the plot
    	p <- ggplot(plot_data, aes(x = .data[[input$x_axis]], y = .data[[input$y_axis]])) +
      	geom_point(aes(color = display_color), alpha = 1,
             size = max(0.3, min(4, 6 / log10(nrow(plot_data) + 1)))) +
		scale_color_manual(values = gate_colors_map, 
                         name = "Gate",
                         drop = FALSE) +
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
  # Helper to format gate summary for text file
  format_gate_summary_text <- function(gate_tree, gate_id, level) {
    gate_info <- gate_tree[[gate_id]]
    indent <- paste(rep("  ", level), collapse = "")
    
    # Build path
    path <- get_gate_path(gate_tree, gate_id)
    path_str <- paste(c("Root", path), collapse = " -> ")
    
    # Calculate percentage of total
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
    lines <- c(lines,sprintf("%s  Thresholds: X[%.2f, %.2f], Y[%.2f, %.2f]",
                            indent,
                            gate_info$threshold_xmin, gate_info$threshold_xmax,
                            gate_info$threshold_ymin, gate_info$threshold_ymax))
    
    # Add children
    children <- get_children(gate_tree, gate_id)
    if (length(children) > 0) {
      lines <- c(lines, sprintf("%s  Children: %s", indent, paste(children, collapse = ", ")))
      for (child in children) {
        lines <- c(lines, format_gate_summary_text(gate_tree, child, level + 1))
      }
    }
    
    lines <- c(lines, "")  # Empty line between gates
    return(lines)
  }
  
  # Helper function to create gate plot
  # Update helper function
  create_gate_plot <- function(data, x_var, y_var, brush, gate_name, inside, x_label = x_var, y_label = y_var) {
  	data$selected <- inside

  	p <- ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]])) +
    	geom_point(aes(color = selected), alpha = 1,
             size = max(0.3, min(4, 6 / log10(nrow(data) + 1)))) +
		scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"),
                       labels = c("Excluded", "Included"),
                       name = "Gate Status") +
    	annotate("rect",
             xmin = brush$xmin, xmax = brush$xmax,
             ymin = brush$ymin, ymax = brush$ymax,
             fill = NA, color = "red", linewidth = 1.5, linetype = "dashed") +
    	labs(title = gate_name,
         subtitle = sprintf("%s vs %s\n%d cells gated", x_label, y_label, sum(inside)),
         x = x_label,
         y = y_label) +
    	theme_bw() +
    	theme(text = element_text(size = 14),
          legend.position = "bottom")

  	return(p)
  }
})
