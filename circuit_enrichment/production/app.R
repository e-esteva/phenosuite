# Circuit Enrichment — PhenoSuite Shiny App
# N-way spatial co-localization analysis for multiplexed tissue imaging.
#
# Architecture: Python engine (circuit_engine.py) via reticulate,
# matching PhenoSuite's existing pattern (masquerade, image processing).
# Single source of truth for the analytical core — same Python code
# runs under both the CLI (sbatch) and the GUI (Shiny).
#
# Conventions:
#   - dplyr:: namespace prefixes throughout
#   - as.data.frame() for any colData access; S4 [[ for hyphenated columns
#   - data.table::fwrite for large outputs
#   - Full self-contained script

suppressPackageStartupMessages({
  library(shiny)
  library(SpatialExperiment)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(S4Vectors)
  library(ggplot2)
  library(patchwork)
  library(reticulate)
  library(jsonlite)
  library(data.table)
  library(scales)
  library(shinycssloaders)
})

source('/srv/shiny-server/phenomenalist/utils/provenance.R')

# ============================================================================
# RETICULATE SETUP  (matches masquerade / pcf-v2 convention)
# ============================================================================

use_virtualenv("r-reticulate", required = FALSE)
source_python("/srv/shiny-server/phenomenalist/circuit_enrichment/production/circuit_engine.py")

# ============================================================================
# BRIDGE HELPERS
# ============================================================================

to_py_array <- function(mat) reticulate::np_array(mat)
to_py_str   <- function(v)   as.character(v)
from_py     <- function(x)   reticulate::py_to_r(x)

# Converts a Python list-of-dicts into an R data.frame
from_py_lod <- function(x) {
  do.call(rbind, lapply(reticulate::py_to_r(x),
    as.data.frame, stringsAsFactors = FALSE))
}

# ============================================================================
# UI
# ============================================================================

ui <- fluidPage(
  tags$head(tags$style(HTML("
    body { font-family: 'Helvetica Neue', Helvetica, sans-serif; }
    .sidebar { background: #f8f9fa; padding: 15px; border-radius: 4px; }
    .circuit-tag {
      display: inline-block; padding: 3px 10px; margin: 2px;
      background: #2c3e50; color: white; border-radius: 12px; font-size: 12px;
    }
    .stat-box { background: #ecf0f1; padding: 12px; border-radius: 6px; margin-bottom: 10px; }
    .stat-val   { font-size: 22px; font-weight: bold; color: #2c3e50; }
    .stat-label { font-size: 11px; color: #7f8c8d; text-transform: uppercase; }
    .engine-badge {
      display: inline-block; padding: 2px 8px; background: #3498db;
      color: white; border-radius: 4px; font-size: 10px;
    }
  "))),

  titlePanel("Circuit Enrichment"),
  tags$span(class = "engine-badge", "Python engine via reticulate"),

  sidebarLayout(
    sidebarPanel(width = 4, class = "sidebar",
      fileInput("spe_file", "Upload SPE (.rds)", accept = ".rds"),
      uiOutput("celltype_col_ui"),

      h4("Define circuit"),
      uiOutput("circuit_ui"),
      uiOutput("circuit_tags"),

      hr(), h4("Neighborhood"),
      numericInput("radius", "Radius (coord units)", 50, min = 1, step = 5),
      selectInput("score_method", "Scoring",
        choices = c("Min fraction"    = "min_fraction",
                    "Geometric mean"  = "geometric_mean")),
      numericInput("n_perm", "Permutations", 200, min = 50, max = 2000, step = 50),

      hr(), h4("Threshold sweep"),
      sliderInput("thresh_range", "Range", 0, 0.5, c(0.01, 0.3), step = 0.005),
      numericInput("n_thresholds", "Steps", 30, min = 5, max = 100),

      hr(),
      actionButton("run", "Run", class = "btn-success btn-block"),

      hr(),
      uiOutput("thresh_slider_ui"),

      hr(),
      downloadButton("dl_zip", "Download (.zip)")
    ),

    mainPanel(width = 8,
      tabsetPanel(id = "tabs",
        tabPanel("Data summary",
          fluidRow(
            column(6, plotOutput("spatial_celltypes", height = "500px") %>% withSpinner()),
            column(6, plotOutput("celltype_bar",      height = "500px") %>% withSpinner())
          )
        ),
        tabPanel("Circuit scores",
          fluidRow(
            column(4, uiOutput("stat_boxes")),
            column(8, plotOutput("score_histogram", height = "350px") %>% withSpinner())
          ),
          plotOutput("spatial_scores", height = "550px") %>% withSpinner()
        ),
        tabPanel("Threshold optimization",
          plotOutput("sweep_plot",       height = "350px") %>% withSpinner(),
          plotOutput("spatial_threshold", height = "550px") %>% withSpinner()
        ),
        tabPanel("Circuit domains",
          plotOutput("domain_plot", height = "600px") %>% withSpinner(),
          verbatimTextOutput("domain_summary")
        )
      )
    )
  )
)

# ============================================================================
# SERVER
# ============================================================================

server <- function(input, output, session) {

  rv <- reactiveValues(
    spe       = NULL,
    celltypes = NULL,
    xy        = NULL,
    comp      = NULL,
    scores    = NULL,
    ztest     = NULL,
    sweep     = NULL
  )

  # --- Provenance tracker -----------------------------------------------
  tempdir0 <- file.path(tempdir(), as.integer(Sys.time()))
  dir.create(tempdir0, showWarnings = FALSE)
  tracker <- ProvenanceTracker$new("circuit_enrichment", session, tempdir0)

  # --- Load SPE ---------------------------------------------------------
  observeEvent(input$spe_file, {
    req(input$spe_file)
    rv$spe <- readRDS(input$spe_file$datapath)
    tracker$register_input(input$spe_file, input_id = "spe_file")
  })

  output$celltype_col_ui <- renderUI({
    req(rv$spe)
    cd <- as.data.frame(colData(rv$spe))
    cc <- names(cd)[vapply(cd, function(x) is.character(x) || is.factor(x), logical(1))]
    selectInput("ct_col", "Cell-type column", choices = cc)
  })

  output$circuit_ui <- renderUI({
    req(rv$spe, input$ct_col)
    cd    <- as.data.frame(colData(rv$spe))
    types <- sort(unique(as.character(cd[[input$ct_col]])))
    selectizeInput("circuit_members", "Circuit members (3+)",
      choices  = types,
      multiple = TRUE,
      options  = list(placeholder = "Pick 3+ cell types"))
  })

  output$circuit_tags <- renderUI({
    req(input$circuit_members)
    tags$div(lapply(input$circuit_members, function(ct)
      tags$span(class = "circuit-tag", ct)))
  })

  # --- Summary plots ----------------------------------------------------
  output$spatial_celltypes <- renderPlot({
    req(rv$spe, input$ct_col)
    cd <- as.data.frame(colData(rv$spe))
    xy <- as.data.frame(spatialCoords(rv$spe)); names(xy) <- c("x", "y")
    xy$celltype <- as.character(cd[[input$ct_col]])
    ggplot(xy, aes(x, y, color = celltype)) +
      geom_point(size = 0.3, alpha = 0.6) +
      coord_equal() + theme_minimal() +
      theme(legend.position = "bottom", legend.text = element_text(size = 7)) +
      guides(color = guide_legend(override.aes = list(size = 2), ncol = 3)) +
      labs(title = "Cell types in tissue space", color = NULL)
  })

  output$celltype_bar <- renderPlot({
    req(rv$spe, input$ct_col)
    cd <- as.data.frame(colData(rv$spe))
    ct <- as.character(cd[[input$ct_col]])
    df <- as.data.frame(table(celltype = ct)); names(df) <- c("celltype", "count")
    df <- df[order(-df$count), ]
    df$celltype   <- factor(df$celltype, levels = df$celltype)
    df$in_circuit <- df$celltype %in% input$circuit_members
    ggplot(df, aes(celltype, count, fill = in_circuit)) +
      geom_col() +
      scale_fill_manual(values = c("FALSE" = "#bdc3c7", "TRUE" = "#2c3e50"),
                        labels = c("Other", "In circuit")) +
      coord_flip() + theme_minimal() +
      labs(title = "Cell-type composition", x = NULL, y = "Count", fill = NULL)
  })

  # --- Main analysis: Python engine via reticulate ----------------------
  observeEvent(input$run, {
    req(rv$spe, input$ct_col, input$circuit_members,
        length(input$circuit_members) >= 3)

    tracker$capture_parameters(input)
    tracker$analysis_started()

    seed_val <- 42L

    cd        <- as.data.frame(colData(rv$spe))
    xy        <- as.matrix(spatialCoords(rv$spe))
    celltypes <- as.character(cd[[input$ct_col]])
    rv$xy        <- xy
    rv$celltypes <- celltypes

    withProgress(message = "Python: neighborhoods...", value = 0.1, {

      comp <- from_py(compute_neighborhood_composition(
        to_py_array(xy),
        to_py_str(celltypes),
        as.numeric(input$radius),
        as.list(input$circuit_members)))
      rv$comp <- comp
      setProgress(0.4, message = "Python: scoring...")

      scores <- as.numeric(from_py(circuit_score(
        to_py_array(comp),
        method = input$score_method)))
      rv$scores <- scores
      setProgress(0.6, message = "Python: permutation test...")

      rv$ztest <- from_py(circuit_zscore(
        to_py_array(scores),
        to_py_array(comp),
        as.integer(input$n_perm),
        as.integer(seed_val)))
      setProgress(0.8, message = "Python: threshold sweep...")

      thresholds <- seq(input$thresh_range[1], input$thresh_range[2],
                        length.out = input$n_thresholds)
      rv$sweep <- from_py_lod(threshold_sweep(
        to_py_array(scores),
        as.list(thresholds)))
      setProgress(1.0)
    })

    # analysis_completed() (with provenance write) is deferred to the
    # download handler so output-file hashes are included in provenance.json.
    updateTabsetPanel(session, "tabs", selected = "Circuit scores")
  })

  # --- Stat boxes -------------------------------------------------------
  output$stat_boxes <- renderUI({
    req(rv$ztest)
    z <- rv$ztest
    tags$div(
      tags$div(class = "stat-box",
        tags$div(class = "stat-val",   sprintf("%.2f",  z$z)),
        tags$div(class = "stat-label", "Z-score")),
      tags$div(class = "stat-box",
        tags$div(class = "stat-val",   sprintf("%.4f",  z$p_value)),
        tags$div(class = "stat-label", "P-value")),
      tags$div(class = "stat-box",
        tags$div(class = "stat-val",   sprintf("%.4f",  z$obs_mean)),
        tags$div(class = "stat-label", "Observed mean")),
      tags$div(class = "stat-box",
        tags$div(class = "stat-val",
          sprintf("%.4f \u00b1 %.4f", z$null_mean, z$null_sd)),
        tags$div(class = "stat-label", "Null mean \u00b1 SD")),
      tags$div(class = "stat-box",
        tags$div(class = "stat-val",   sum(rv$scores > 0)),
        tags$div(class = "stat-label",
          paste0("Active (",
            sprintf("%.1f%%", 100 * mean(rv$scores > 0)), ")")))
    )
  })

  output$score_histogram <- renderPlot({
    req(rv$scores)
    ggplot(data.frame(s = rv$scores[rv$scores > 0]), aes(s)) +
      geom_histogram(bins = 60, fill = "#2c3e50", color = "white", linewidth = 0.2) +
      theme_minimal() +
      labs(title = "Score distribution (non-zero)", x = "Score", y = "Count")
  })

  output$spatial_scores <- renderPlot({
    req(rv$scores, rv$xy)
    df <- data.frame(x = rv$xy[,1], y = rv$xy[,2], score = rv$scores)
    ggplot(df, aes(x, y, color = score)) +
      geom_point(size = 0.3) +
      scale_color_viridis_c(option = "inferno") +
      coord_equal() + theme_minimal() +
      labs(title = paste("Circuit:",
             paste(input$circuit_members, collapse = " + ")),
           color = "Score")
  })

  # --- Threshold sweep --------------------------------------------------
  output$sweep_plot <- renderPlot({
    req(rv$sweep)
    df <- rv$sweep
    p1 <- ggplot(df, aes(threshold, n_positive)) +
      geom_line(linewidth = 1, color = "#2c3e50") +
      geom_point(size = 2, color = "#2c3e50") +
      theme_minimal() +
      labs(title = "Threshold sweep", x = "Threshold", y = "N positive")
    p2 <- ggplot(df, aes(threshold, frac_positive)) +
      geom_line(linewidth = 1, color = "#e74c3c") +
      geom_point(size = 2, color = "#e74c3c") +
      theme_minimal() +
      labs(x = "Threshold", y = "Fraction positive")
    p1 | p2
  })

  output$thresh_slider_ui <- renderUI({
    req(rv$sweep)
    rng <- range(rv$sweep$threshold)
    sliderInput("active_thresh", "Explore threshold",
      min = rng[1], max = rng[2],
      value = median(rv$sweep$threshold), step = 0.005)
  })

  output$spatial_threshold <- renderPlot({
    req(rv$scores, rv$xy, input$active_thresh)
    df           <- data.frame(x = rv$xy[,1], y = rv$xy[,2], score = rv$scores)
    df$positive  <- df$score >= input$active_thresh
    n_pos        <- sum(df$positive)
    frac         <- sprintf("%.1f%%", 100 * mean(df$positive))
    ggplot(df, aes(x, y, color = positive)) +
      geom_point(size = 0.3) +
      scale_color_manual(values = c("FALSE" = "#ecf0f1", "TRUE" = "#e74c3c"),
                         labels = c("Below", "Circuit-positive")) +
      coord_equal() + theme_minimal() +
      theme(legend.position = "bottom") +
      labs(title = paste0("Threshold = ", input$active_thresh,
             "  |  ", n_pos, " cells  (", frac, ")"),
           color = NULL)
  })

  # --- Circuit domains --------------------------------------------------
  output$domain_plot <- renderPlot({
    req(rv$scores, rv$xy, rv$celltypes, input$active_thresh, input$circuit_members)
    df           <- data.frame(x = rv$xy[,1], y = rv$xy[,2],
                               celltype = rv$celltypes, score = rv$scores)
    df$positive  <- df$score >= input$active_thresh
    df$label     <- ifelse(df$positive & df$celltype %in% input$circuit_members,
                           df$celltype, "other")
    df$label     <- factor(df$label, levels = c(input$circuit_members, "other"))
    n_ct  <- length(input$circuit_members)
    pal   <- c(setNames(scales::hue_pal()(n_ct), input$circuit_members),
               "other" = "#ecf0f1")
    ggplot(df, aes(x, y, color = label)) +
      geom_point(size = 0.3) +
      scale_color_manual(values = pal) +
      coord_equal() + theme_minimal() +
      theme(legend.position = "bottom") +
      guides(color = guide_legend(override.aes = list(size = 3))) +
      labs(title    = paste0("Circuit domains (threshold = ", input$active_thresh, ")"),
           subtitle = paste("Circuit:", paste(input$circuit_members, collapse = " + ")),
           color    = NULL)
  })

  output$domain_summary <- renderPrint({
    req(rv$scores, rv$ztest, input$active_thresh, input$circuit_members)
    pos <- rv$scores >= input$active_thresh
    cat("Engine: circuit_engine.py via reticulate\n")
    cat("Circuit:", paste(input$circuit_members, collapse = " + "), "\n")
    cat("Scoring:", input$score_method, "\n")
    cat("Radius:", input$radius, "\n")
    cat("Threshold:", input$active_thresh, "\n")
    cat("Positive:", sum(pos), "/", length(pos),
        sprintf("(%.1f%%)", 100 * mean(pos)), "\n")
    cat("Z-score:", sprintf("%.3f", rv$ztest$z), "\n")
    cat("P-value:", sprintf("%.4f", rv$ztest$p_value), "\n")
  })

  # --- Download ---------------------------------------------------------
  output$dl_zip <- downloadHandler(
    filename = function()
      paste0("circuit_enrichment_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".zip"),
    content = function(file) {
      dir.create(tempdir0, showWarnings = FALSE, recursive = TRUE)

      sp <- file.path(tempdir0, "circuit_scores.csv")
      data.table::fwrite(data.frame(
        x             = rv$xy[,1],
        y             = rv$xy[,2],
        celltype      = rv$celltypes,
        circuit_score = rv$scores,
        positive      = rv$scores >= input$active_thresh), sp)

      cp <- file.path(tempdir0, "neighborhood_composition.csv")
      cd <- as.data.frame(rv$comp)
      colnames(cd) <- input$circuit_members
      cd$x <- rv$xy[,1]; cd$y <- rv$xy[,2]
      data.table::fwrite(cd, cp)

      tp <- file.path(tempdir0, "threshold_sweep.csv")
      data.table::fwrite(rv$sweep, tp)

      jp <- file.path(tempdir0, "enrichment_summary.json")
      write(jsonlite::toJSON(list(
        circuit       = input$circuit_members,
        radius        = input$radius,
        score_method  = input$score_method,
        n_perm        = input$n_perm,
        threshold     = input$active_thresh,
        n_cells       = length(rv$scores),
        n_positive    = sum(rv$scores >= input$active_thresh),
        frac_positive = mean(rv$scores >= input$active_thresh),
        z             = rv$ztest$z,
        p             = rv$ztest$p_value,
        engine        = "circuit_engine.py via reticulate"),
        pretty = TRUE, auto_unbox = TRUE), jp)

      # Provenance sidecar: set output_dir and close the analysis;
      # analysis_completed() calls write_provenance() which hashes every
      # file now present in tempdir0 and writes provenance.json + replay.R.
      tracker$output_dir <- tempdir0
      tracker$analysis_completed()

      zip::zip(zipfile = file, files = dir(tempdir0), root = tempdir0)
    },
    contentType = "application/zip"
  )
}

shinyApp(ui, server)
