library(shiny)
library(SpatialExperiment)
library(SingleCellExperiment)
library(ggplot2)
library(colourpicker)
library(plotly)
source('/srv/shiny-server/phenomenalist/utils/provenance.R')

options(shiny.maxRequestSize = 4*500 * 1024^2)

# ─────────────────────────────────────────────────────────────────────────────
# UI
# ─────────────────────────────────────────────────────────────────────────────
ui <- fluidPage(
  tags$head(
    tags$title("SpatialExploreR"),
    tags$style(HTML("
    body { background: #f7f8fa; font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; }
    .sidebar-panel  { background: #fff; border-radius: 8px; padding: 20px;
                      box-shadow: 0 1px 4px rgba(0,0,0,.08); }
    .main-panel     { background: #fff; border-radius: 8px; padding: 20px;
                      box-shadow: 0 1px 4px rgba(0,0,0,.08); }
    .section-title  { font-weight: 600; margin-bottom: 10px; color: #333; }
    .mutate-box     { background: #f0f4ff; border: 1px solid #c5d2f0;
                      border-radius: 6px; padding: 15px; margin-top: 12px; }
    .col-entry      { background: #fff; border: 1px solid #dce3f0;
                      border-radius: 4px; padding: 8px 10px; margin-bottom: 6px; }
    hr.sep          { border-top: 1px solid #e0e0e0; }
    .btn-primary    { background: #4a6cf7; border: none; }
    .btn-primary:hover { background: #3958d9; }
    .btn-success    { background: #34b370; border: none; }
    .btn-danger     { background: #e55353; border: none; }
    #plot_area, #plot_static, #plot_interactive { min-height: 500px; }
  "))),
  
  titlePanel(
    div(
      tags$span("SpatialExploreR",
                style = "font-weight:700; font-size:1.5em;"),
      tags$span(" — upload, visualise & mutate",
                style = "color:#888; font-size:0.9em;")
    )
  ),
  
  sidebarLayout(
    
    # ── Sidebar ──────────────────────────────────────────────────────────────
    sidebarPanel(
      class = "sidebar-panel", width = 4,
      
      # Upload
      h5(class = "section-title", "1. Upload SpatialExperiment (.rds)"),
      fileInput("rds_file", NULL, accept = ".rds",
                placeholder = "Choose .rds file"),
      uiOutput("upload_status"),
      hr(class = "sep"),
      
      # Plot type + DR + colour
      h5(class = "section-title", "2. Plot Settings"),
      uiOutput("plot_type_ui"),
      uiOutput("dr_ui"),
      uiOutput("colour_ui"),
      uiOutput("palette_ui"),
      uiOutput("feature_toggle_ui"),
      uiOutput("assay_select_ui"),
      uiOutput("feature_select_ui"),
      hr(class = "sep"),
      
      # ── Mutate ─────────────────────────────────────────────────────────────
      h5(class = "section-title", "3. Create New Column"),
      div(
        class = "mutate-box",
        textInput("new_col_name", "New column name", placeholder = "e.g. my_score"),
        uiOutput("mutate_entries_ui"),
        fluidRow(
          column(6, actionButton("add_col_btn", "Add Column",
                                 icon = icon("plus"), class = "btn-primary btn-sm btn-block")),
          column(6, actionButton("create_col_btn", "Create New Column",
                                 icon = icon("wand-magic-sparkles"), class = "btn-success btn-sm btn-block"))
        ),
        uiOutput("mutate_status")
      ),
      
      hr(class = "sep"),
      
      # Downloads
      h5(class = "section-title", "4. Downloads"),
      fluidRow(
        column(6, downloadButton("dl_pdf",  "Download PDF",  class = "btn-primary btn-sm btn-block")),
        column(6, uiOutput("dl_rds_ui"))
      )
    ),
    
    # ── Main panel ───────────────────────────────────────────────────────────
    mainPanel(
      class = "main-panel", width = 8,
      uiOutput("plot_container"),
      verbatimTextOutput("obj_summary")
    )
  )
)

# ─────────────────────────────────────────────────────────────────────────────
# Server
# ─────────────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {
  prov_dir <- file.path(tempdir(), paste0("spatialExploreR_", session$token))
  dir.create(prov_dir, showWarnings = FALSE)
  tracker <- ProvenanceTracker$new("spatialExploreR", session, prov_dir)

  # Reactive: loaded SPE object (mutable copy) --------------------------------
  spe <- reactiveVal(NULL)
  mutated <- reactiveVal(FALSE)      # has the user created any new columns?
  n_entries <- reactiveVal(1)        # number of column entries in mutate UI
  
  # ── Upload ─────────────────────────────────────────────────────────────────
  observeEvent(input$rds_file, {
    req(input$rds_file)
    tryCatch({
      obj <- readRDS(input$rds_file$datapath)
      stopifnot(is(obj, "SpatialExperiment") || is(obj, "SingleCellExperiment"))
      spe(obj)
      mutated(FALSE)
      n_entries(1)
    }, error = function(e) {
      showNotification(paste("Error loading file:", e$message),
                       type = "error", duration = 8)
    })
  })
  
  output$upload_status <- renderUI({
    req(spe())
    tags$span(style = "color:#34b370; font-weight:600;",
              icon("check-circle"),
              sprintf("Loaded: %d cells × %d features",
                      ncol(spe()), nrow(spe())))
  })
  
  # ── Helpers ────────────────────────────────────────────────────────────────
  dr_names <- reactive({
    req(spe())
    reducedDimNames(spe())
  })
  
  meta_cols <- reactive({
    req(spe())
    names(colData(spe()))
  })
  
  numeric_cols <- reactive({
    req(spe())
    cd <- colData(spe())
    nms <- names(cd)
    nms[vapply(nms, function(n) is.numeric(cd[[n]]), logical(1))]
  })
  
  # ── Plot type selector ──────────────────────────────────────────────────
  has_spatial <- reactive({
    req(spe())
    is(spe(), "SpatialExperiment") &&
      !is.null(tryCatch(spatialCoords(spe()), error = function(e) NULL)) &&
      nrow(spatialCoords(spe())) > 0
  })
  
  output$plot_type_ui <- renderUI({
    req(spe())
    choices <- c("Reductions")
    if (has_spatial()) choices <- c("Spatial", choices)
    selectInput("plot_type", "Plot type", choices = choices, selected = choices[1])
  })
  
  # ── DR selector (only shown for Reductions) ────────────────────────────────
  output$dr_ui <- renderUI({
    req(dr_names(), input$plot_type)
    if (input$plot_type != "Reductions") return(NULL)
    selectInput("dr_choice", "Dimensionality Reduction",
                choices = dr_names(), selected = dr_names()[1])
  })
  
  # ── Colour-by selector ────────────────────────────────────────────────────
  output$colour_ui <- renderUI({
    req(meta_cols())
    selectInput("colour_col", "Colour by (metadata column)",
                choices = meta_cols(), selected = meta_cols()[1])
  })
  
  # Show palette selector for continuous vs discrete
  output$palette_ui <- renderUI({
    req(spe())
    # Show palette when colouring by numeric metadata OR by a feature
    show_palette <- FALSE
    if (isTRUE(input$plot_features == "Yes")) {
      show_palette <- TRUE
    } else if (!is.null(input$colour_col)) {
      vals <- colData(spe())[[input$colour_col]]
      if (is.numeric(vals)) show_palette <- TRUE
    }
    if (show_palette) {
      selectInput("cont_palette", "Continuous palette",
                  choices = c("viridis", "inferno", "plasma", "magma",
                              "cividis", "RdYlBu", "Spectral"),
                  selected = "viridis")
    } else {
      NULL
    }
  })
  
  # ── Feature plotting toggle ────────────────────────────────────────────────
  output$feature_toggle_ui <- renderUI({
    req(spe())
    radioButtons("plot_features", "Plot features?",
                 choices = c("No", "Yes"), selected = "No", inline = TRUE)
  })
  
  output$assay_select_ui <- renderUI({
    req(spe(), input$plot_features)
    if (input$plot_features != "Yes") return(NULL)
    assay_names <- assayNames(spe())
    selectInput("assay_choice", "Assay",
                choices = assay_names, selected = assay_names[1])
  })
  
  output$feature_select_ui <- renderUI({
    req(spe(), input$plot_features, input$assay_choice)
    if (input$plot_features != "Yes") return(NULL)
    selectizeInput("feature_name", "Feature name",
                   choices = NULL, selected = NULL,
                   options = list(placeholder = "Type to search…"))
  })
  
  # Server-side autocomplete for feature names (re-fires when assay changes)
  observeEvent(list(input$plot_features, input$assay_choice), {
    req(input$plot_features == "Yes", spe(), input$assay_choice)
    features <- rownames(spe())
    updateSelectizeInput(session, "feature_name",
                         choices = features, selected = NULL,
                         server = TRUE)
  })
  
  # ── Is the current plot discrete? ──────────────────────────────────────────
  is_discrete <- reactive({
    req(spe())
    # Feature plotting is always continuous
    if (isTRUE(input$plot_features == "Yes") && !is.null(input$feature_name) &&
        input$feature_name != "") {
      return(FALSE)
    }
    req(input$colour_col)
    vals <- colData(spe())[[input$colour_col]]
    is.factor(vals) || is.character(vals) || is.logical(vals)
  })
  
  # ── Dynamic plot container ─────────────────────────────────────────────────
  output$plot_container <- renderUI({
    req(spe())
    if (is_discrete()) {
      plotlyOutput("plot_interactive", height = "600px")
    } else {
      plotOutput("plot_static", height = "600px")
    }
  })
  
  # ── Main ggplot (used by both renderers & PDF download) ────────────────────
  make_plot <- reactive({
    req(spe(), input$plot_type)
    obj <- spe()
    
    # Determine colour source: feature expression or metadata column
    use_feature <- isTRUE(input$plot_features == "Yes") &&
      !is.null(input$assay_choice) &&
      !is.null(input$feature_name) && input$feature_name != ""
    
    if (use_feature) {
      feat <- input$feature_name
      if (!(feat %in% rownames(obj))) {
        showNotification(paste0("Feature '", feat, "' not found in assay."),
                         type = "warning")
        return(NULL)
      }
      # Pull expression from the selected assay
      colour_vec <- as.numeric(assay(obj, input$assay_choice)[feat, ])
      colour_label <- feat
    } else {
      req(input$colour_col)
      cc <- input$colour_col
      # Pull directly from S4 colData — preserves hyphens unlike as.data.frame()
      colour_vec <- colData(obj)[[cc]]
      colour_label <- cc
    }
    
    if (input$plot_type == "Spatial") {
      # ── Spatial coordinates plot ──────────────────────────────────────────
      req(has_spatial())
      sc <- as.data.frame(spatialCoords(obj))
      colnames(sc)[1:2] <- c("X", "Y")
      df <- data.frame(X = sc[, 1], Y = sc[, 2], .colour_by. = colour_vec,
                       check.names = FALSE)
      
      p <- ggplot(df, aes(x = X, y = Y, colour = .colour_by.)) +
        geom_point(size = 0.6, alpha = 0.75) +
        labs(x = "Spatial X", y = "Spatial Y",
             colour = colour_label,
             title = paste0("Spatial — coloured by ", colour_label)) +
        coord_fixed() +
        scale_y_reverse() +
        theme_minimal(base_size = 14) +
        theme(
          plot.title       = element_text(face = "bold", size = 15),
          legend.position  = "right",
          panel.grid.minor = element_blank()
        )
      
    } else {
      # ── Dimensionality reduction plot ─────────────────────────────────────
      req(input$dr_choice)
      dr <- input$dr_choice
      
      coords <- as.data.frame(reducedDim(obj, dr))
      if (ncol(coords) < 2) return(NULL)
      colnames(coords)[1:2] <- c("Dim1", "Dim2")
      df <- data.frame(Dim1 = coords[, 1], Dim2 = coords[, 2],
                       .colour_by. = colour_vec, check.names = FALSE)
      
      p <- ggplot(df, aes(x = Dim1, y = Dim2, colour = .colour_by.)) +
        geom_point(size = 0.6, alpha = 0.75) +
        labs(x = paste0(dr, " 1"), y = paste0(dr, " 2"),
             colour = colour_label,
             title = paste0(dr, " — coloured by ", colour_label)) +
        theme_minimal(base_size = 14) +
        theme(
          plot.title       = element_text(face = "bold", size = 15),
          legend.position  = "right",
          panel.grid.minor = element_blank()
        )
    }
    
    # ── Colour scale ────────────────────────────────────────────────────────
    is_num <- is.numeric(df[[".colour_by."]])
    if (is_num) {
      pal <- if (!is.null(input$cont_palette)) input$cont_palette else "viridis"
      if (pal %in% c("viridis", "inferno", "plasma", "magma", "cividis")) {
        p <- p + scale_colour_viridis_c(option = pal)
      } else {
        p <- p + scale_colour_distiller(palette = pal)
      }
    }
    
    p
  })
  
  # ── Plotly renderer (discrete / factor columns) ────────────────────────────
  output$plot_interactive <- renderPlotly({
    p <- make_plot()
    req(p)
    
    obj <- spe()
    cc  <- input$colour_col
    cell_ids <- colnames(obj)
    
    # Pull colour values from the ggplot's own data (already safe)
    plot_data <- ggplot_build(p)$plot$data
    colour_vals <- plot_data[[".colour_by."]]
    hover_text <- paste0("cell: ", cell_ids, "\n", cc, ": ", colour_vals)
    
    # Convert ggplot → plotly with custom tooltip
    pp <- ggplotly(p, tooltip = "none") %>%
      style(
        text = hover_text,
        hoverinfo = "text"
      ) %>%
      layout(
        hoverlabel = list(bgcolor = "white", font = list(size = 12)),
        dragmode   = "zoom",
        legend     = list(itemsizing = "constant")
      ) %>%
      toWebGL()   # GPU-accelerated for large datasets
    
    pp
  })
  
  # ── Static renderer (numeric / continuous columns) ─────────────────────────
  output$plot_static <- renderPlot({ make_plot() })
  
  output$obj_summary <- renderPrint({
    req(spe())
    obj <- spe()
    cat("Object class :", class(obj), "\n")
    cat("Dimensions   :", nrow(obj), "features ×", ncol(obj), "cells\n")
    cat("Reduced dims :", paste(reducedDimNames(obj), collapse = ", "), "\n")
    cat("Metadata cols:", paste(names(colData(obj)), collapse = ", "), "\n")
  })
  
  # ── Mutate: dynamic column entries ─────────────────────────────────────────
  # Each entry = (column selector) + (operation to apply BEFORE the next col)
  # The chain is evaluated left-to-right: col1 OP1 col2 OP2 col3 ...
  
  output$mutate_entries_ui <- renderUI({
    req(meta_cols())
    n <- n_entries()
    num_cols <- numeric_cols()
    all_cols <- meta_cols()
    
    entries <- lapply(seq_len(n), function(i) {
      col_id <- paste0("mut_col_", i)
      op_id  <- paste0("mut_op_",  i)
      
      col_row <- div(
        class = "col-entry",
        fluidRow(
          column(
            if (i < n) 6 else 12,
            selectInput(col_id,
                        if (i == 1) "Column" else paste0("Column ", i),
                        choices = all_cols, selected = NULL)
          ),
          if (i < n) {
            column(6,
                   selectInput(op_id, "then…",
                               choices = c("+" = "+", "−" = "-",
                                           "×" = "*", "÷" = "/",
                                           "AND (&)" = "&", "OR (|)" = "|",
                                           "==" = "==", "!=" = "!=",
                                           ">" = ">", "<" = "<",
                                           ">=" = ">=", "<=" = "<=",
                                           "paste" = "paste",
                                           "pmax" = "pmax", "pmin" = "pmin"),
                               selected = "+")
            )
          }
        )
      )
      col_row
    })
    tagList(entries)
  })
  
  observeEvent(input$add_col_btn, {
    n_entries(n_entries() + 1)
  })
  
  # ── Mutate: create the new column ──────────────────────────────────────────
  observeEvent(input$create_col_btn, {
    req(spe())
    new_name <- trimws(input$new_col_name)
    if (nchar(new_name) == 0) {
      showNotification("Please provide a name for the new column.",
                       type = "warning"); return()
    }
    
    n   <- n_entries()
    obj <- spe()
    cd  <- colData(obj)  # S4 DataFrame — preserves hyphenated names
    
    # Collect column names & operators
    cols <- character(n)
    ops  <- character(n - 1)
    for (i in seq_len(n)) {
      cols[i] <- input[[paste0("mut_col_", i)]]
      if (i < n) ops[i] <- input[[paste0("mut_op_", i)]]
    }
    
    if (any(is.null(cols)) || any(cols == "")) {
      showNotification("Select a column for every entry.", type = "warning")
      return()
    }
    
    tryCatch({
      # Start with the first column's values
      result <- cd[[cols[1]]]
      
      for (i in seq_along(ops)) {
        rhs <- cd[[cols[i + 1]]]
        op  <- ops[i]
        
        result <- switch(op,
                         "+"     = as.numeric(result) + as.numeric(rhs),
                         "-"     = as.numeric(result) - as.numeric(rhs),
                         "*"     = as.numeric(result) * as.numeric(rhs),
                         "/"     = as.numeric(result) / as.numeric(rhs),
                         "&"     = as.logical(result) & as.logical(rhs),
                         "|"     = as.logical(result) | as.logical(rhs),
                         "=="    = result == rhs,
                         "!="    = result != rhs,
                         ">"     = as.numeric(result) >  as.numeric(rhs),
                         "<"     = as.numeric(result) <  as.numeric(rhs),
                         ">="    = as.numeric(result) >= as.numeric(rhs),
                         "<="    = as.numeric(result) <= as.numeric(rhs),
                         "paste" = paste(result, rhs, sep = "_"),
                         "pmax"  = pmax(as.numeric(result), as.numeric(rhs), na.rm = TRUE),
                         "pmin"  = pmin(as.numeric(result), as.numeric(rhs), na.rm = TRUE),
                         stop("Unknown operation: ", op)
        )
      }
      
      colData(obj)[[new_name]] <- result
      spe(obj)
      mutated(TRUE)
      
      showNotification(
        paste0("Column '", new_name, "' created successfully (",
               class(result)[1], ", length ", length(result), ")"),
        type = "message", duration = 5)
      
    }, error = function(e) {
      showNotification(paste("Column creation failed:", e$message),
                       type = "error", duration = 8)
    })
  })
  
  output$mutate_status <- renderUI({
    req(mutated())
    tags$span(style = "color:#34b370; font-size:0.85em; margin-top:6px; display:block;",
              icon("circle-check"), " Object has been modified — download available.")
  })
  
  # ── Downloads ──────────────────────────────────────────────────────────────
  output$dl_pdf <- downloadHandler(
    filename = function() {
      plot_label <- if (input$plot_type == "Spatial") "Spatial" else input$dr_choice
      use_feature <- isTRUE(input$plot_features == "Yes") &&
        !is.null(input$feature_name) && input$feature_name != ""
      colour_label <- if (use_feature) input$feature_name else input$colour_col
      paste0(plot_label, "_", colour_label, "_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      p <- make_plot()
      ggsave(file, plot = p, device = "pdf", width = 10, height = 7)
    }
  )
  
  output$dl_rds_ui <- renderUI({
    if (isTRUE(mutated())) {
      downloadButton("dl_rds", "Download .rds", class = "btn-success btn-sm btn-block")
    } else {
      tags$button("Download .rds", class = "btn btn-secondary btn-sm btn-block",
                  disabled = "disabled", title = "Mutate the object first")
    }
  })
  
  output$dl_rds <- downloadHandler(
    filename = function() {
      paste0("spe_modified_", Sys.Date(), ".rds")
    },
    content = function(file) {
      saveRDS(spe(), file)
    }
  )
}

# ─────────────────────────────────────────────────────────────────────────────
shinyApp(ui, server)
