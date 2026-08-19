require(shiny)
require(shinyjs)
require(glue)
require(SpatialExperiment)
require(SingleCellExperiment)
require(jsonlite)
require(DT)
require(dplyr)
require(zip)
require(ggplot2)
require(cowplot)
source('/srv/shiny-server/phenomenalist/utils/provenance.R')
source('/srv/shiny-server/phenomenalist/utils/anndata_export.R')
source('/srv/shiny-server/phenomenalist/utils/vectra_export.R')
source('/srv/shiny-server/phenomenalist/gemma_phenotyper/production/gemma_utils.R')
source('/srv/shiny-server/phenomenalist/utils/RunPhenomenalist-shiny/RunPhenomenalist-shiny.R')
source('/srv/shiny-server/phenomenalist/utils/RunPhenomenalist-shiny/phenomenalist-utils-shiny.R')
source('/srv/shiny-server/phenomenalist/utils/plot-scatter.R')

options(shiny.maxRequestSize = 4000 * 1024^2)   # 4 GB — for large merged checkpoints
jsResetCode <- "shinyjs.resetClick = function() { history.go(0) }"

INFER_SCRIPT <- '/srv/shiny-server/phenomenalist/gemma_phenotyper/production/infer.py'
GEMMA_PYTHON <- Sys.getenv("PHENOSUITE_GEMMA_PYTHON", "/opt/venv/bin/python3")

# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------
ui <- fluidPage(
  useShinyjs(),
  extendShinyjs(text = jsResetCode, functions = "resetClick"),

  tags$head(tags$style(HTML("
    @import url('https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;600&family=IBM+Plex+Sans:wght@300;400;600&display=swap');

    * { box-sizing: border-box; }

    body {
      font-family: 'IBM Plex Sans', sans-serif;
      font-weight: 400;
      background: #f4f3ef;
      color: #1a1a1a;
      font-size: 14px;
    }

    .app-header {
      background: #1a1a1a;
      color: #f4f3ef;
      padding: 18px 28px 14px;
      margin: -15px -15px 0 -15px;
      display: flex;
      align-items: baseline;
      gap: 14px;
      border-bottom: 3px solid #c8f55a;
    }
    .app-header h1 {
      font-family: 'IBM Plex Mono', monospace;
      font-size: 18px;
      font-weight: 600;
      letter-spacing: -0.3px;
      margin: 0;
      color: #f4f3ef;
    }
    .app-header .subtitle {
      font-size: 12px;
      color: #888;
      letter-spacing: 0.5px;
      text-transform: uppercase;
    }

    .well {
      background: #ffffff !important;
      border: 1px solid #e0ddd6 !important;
      border-radius: 4px !important;
      box-shadow: none !important;
      padding: 0 !important;
    }

    .step-block {
      padding: 14px 16px;
      border-bottom: 1px solid #f0ede6;
    }
    .step-block:last-child { border-bottom: none; }

    .step-label {
      font-family: 'IBM Plex Mono', monospace;
      font-size: 10px;
      font-weight: 600;
      letter-spacing: 1.2px;
      text-transform: uppercase;
      color: #888;
      margin-bottom: 10px;
      display: flex;
      align-items: center;
      gap: 8px;
    }
    .step-num {
      background: #1a1a1a;
      color: #c8f55a;
      font-size: 9px;
      font-weight: 600;
      width: 16px;
      height: 16px;
      border-radius: 50%;
      display: inline-flex;
      align-items: center;
      justify-content: center;
      flex-shrink: 0;
    }

    .form-control {
      border: 1px solid #d4d0c8;
      border-radius: 3px;
      background: #fafaf8;
      font-family: 'IBM Plex Sans', sans-serif;
      font-size: 13px;
      color: #1a1a1a;
      height: 32px;
      padding: 4px 10px;
    }
    .form-control:focus {
      border-color: #1a1a1a;
      box-shadow: none;
      background: #fff;
    }

    textarea.form-control { height: auto; }

    .control-label {
      font-size: 12px;
      font-weight: 600;
      color: #444;
      margin-bottom: 4px;
    }

    .key-status {
      font-family: 'IBM Plex Mono', monospace;
      font-size: 11px;
      padding: 3px 8px;
      border-radius: 3px;
      display: inline-block;
      margin-top: 6px;
    }
    .key-ok   { background: #e8f5e9; color: #2e7d32; }
    .key-err  { background: #ffebee; color: #c62828; }
    .key-idle { background: #f0ede6; color: #888; }

    .btn {
      font-family: 'IBM Plex Sans', sans-serif;
      font-size: 12px;
      font-weight: 600;
      letter-spacing: 0.3px;
      border-radius: 3px;
      border: none;
      cursor: pointer;
      height: 32px;
      padding: 0 14px;
      transition: background 0.15s, opacity 0.15s;
    }
    .btn-load {
      background: #1a1a1a;
      color: #c8f55a;
      width: 100%;
    }
    .btn-load:hover { background: #333; color: #c8f55a; }

    .btn-run {
      background: #c8f55a;
      color: #1a1a1a;
      width: 100%;
      height: 36px;
      font-size: 13px;
    }
    .btn-run:hover { background: #b8e040; }

    .btn-dl {
      background: #f0ede6;
      color: #1a1a1a;
      width: 100%;
      border: 1px solid #d4d0c8;
    }
    .btn-dl:hover { background: #e8e4da; }

    .btn-retry {
      background: #fff3e0;
      color: #b25000;
      width: 100%;
      border: 1px solid #ffcc80;
    }
    .btn-retry:hover { background: #ffe8cc; }
    .failure-note {
      font-family: 'IBM Plex Mono', monospace;
      font-size: 11px;
      color: #b25000;
      background: #fff3e0;
      border-radius: 3px;
      padding: 4px 8px;
      display: block;
      margin-bottom: 6px;
    }

    .btn-reset {
      background: transparent;
      color: #888;
      width: 100%;
      border: 1px solid #d4d0c8;
      font-size: 11px;
    }
    .btn-reset:hover { background: #f0ede6; color: #555; }

    .tab-content { background: #fff; border: 1px solid #e0ddd6; border-top: none; padding: 20px; border-radius: 0 0 4px 4px; }
    .nav-tabs > li > a {
      font-family: 'IBM Plex Mono', monospace;
      font-size: 11px;
      letter-spacing: 0.5px;
      text-transform: uppercase;
      color: #888;
      border-radius: 3px 3px 0 0 !important;
      border: 1px solid #e0ddd6 !important;
      background: #f4f3ef;
      margin-right: 3px;
    }
    .nav-tabs > li.active > a,
    .nav-tabs > li.active > a:focus,
    .nav-tabs > li.active > a:hover {
      color: #1a1a1a !important;
      background: #fff !important;
      border-bottom-color: #fff !important;
      font-weight: 600;
    }
    .nav-tabs { border-bottom: 1px solid #e0ddd6; margin-bottom: 0; }

    .marker-scroll { max-height: 180px; overflow-y: auto; border: 1px solid #e0ddd6; border-radius: 3px; padding: 6px 8px; background: #fafaf8; }
    .checkbox-inline { font-size: 12px; margin-right: 10px; }
    .shiny-input-checkboxgroup > .checkbox-inline { font-size: 12px; }

    /* ── Temperature slider (Ollama) ─────────────────────────────── */
    .temp-slider .irs--shiny .irs-line {
      background: #fafaf8; border: 1px solid #d4d0c8; border-radius: 3px; height: 6px; top: 33px;
    }
    .temp-slider .irs--shiny .irs-bar {
      background: #c8f55a; border-top: 1px solid #b8e040; border-bottom: 1px solid #b8e040; height: 6px; top: 33px;
    }
    .temp-slider .irs--shiny .irs-handle {
      background: #1a1a1a; border: 2px solid #c8f55a; border-radius: 50%;
      width: 16px; height: 16px; top: 27px; box-shadow: none; cursor: pointer;
    }
    .temp-slider .irs--shiny .irs-handle:hover,
    .temp-slider .irs--shiny .irs-handle.state_hover { background: #333; }
    .temp-slider .irs--shiny .irs-handle > i:first-child { display: none; }
    .temp-slider .irs--shiny .irs-single {
      background: #1a1a1a; color: #c8f55a; font-family: 'IBM Plex Mono', monospace;
      font-size: 10px; font-weight: 600; border-radius: 3px; padding: 2px 6px;
    }
    .temp-slider .irs--shiny .irs-single:before { border-top-color: #1a1a1a; }
    .temp-slider .irs--shiny .irs-min, .temp-slider .irs--shiny .irs-max {
      font-family: 'IBM Plex Mono', monospace; font-size: 10px; color: #aaa; background: transparent; top: 33px;
    }
    .temp-slider .irs--shiny .irs-grid { display: none; }
    .temp-scale {
      display: flex; justify-content: space-between; font-family: 'IBM Plex Mono', monospace;
      font-size: 9px; letter-spacing: 0.5px; text-transform: uppercase; color: #aaa; margin-top: 2px;
    }

    .btn-file {
      background: #1a1a1a;
      color: #f4f3ef;
      font-size: 11px;
      font-weight: 600;
      border-radius: 2px;
      padding: 3px 10px;
    }

    .progress-bar { background-color: #c8f55a; }
    .progress { border-radius: 2px; }

    .shiny-notification {
      font-family: 'IBM Plex Mono', monospace;
      font-size: 12px;
      background: #1a1a1a;
      color: #f4f3ef;
      border: none;
      border-left: 3px solid #c8f55a;
      border-radius: 3px;
    }
  "))),

  # ── Header ──────────────────────────────────────────────────────────────
  div(class = "app-header",
    h1("phenomenalist"),
    span(class = "subtitle", "mIF Cell Phenotyper (Gemma)")
  ),

  br(),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      # ── Step 1: Data ────────────────────────────────────────────────
      div(class = "step-block",
        div(class = "step-label", span(class = "step-num", "1"), "Load Data"),
        fileInput("rds_file", NULL, accept = ".rds",
                  placeholder = "Choose .rds file…",
                  buttonLabel = tags$span(icon("folder-open"), "Browse")),
        conditionalPanel(
          condition = "output.has_data",
          selectInput("cluster_col", "Cluster column:", choices = ""),
          div(class = "marker-scroll",
            checkboxGroupInput("marker_cols", NULL, choices = "", inline = FALSE)
          ),
          tags$small(style = "color:#999;", "Select markers to include in prompts")
        )
      ),

      # ── Step 2: Model Weights ────────────────────────────────────────
      div(class = "step-block",
        div(class = "step-label", span(class = "step-num", "2"), "Model Weights"),
        radioButtons("model_source", NULL,
          choices  = c("Upload model zip"       = "upload",
                       "Server directory path"  = "server_path",
                       "Ollama server (zero-shot)" = "ollama"),
          selected = "server_path"
        ),
        conditionalPanel(
          condition = "input.model_source == 'upload'",
          fileInput("model_zip", NULL,
                    accept = c(".zip", "application/zip"),
                    buttonLabel = tags$span(icon("file-archive"), "Choose zip"))
        ),
        conditionalPanel(
          condition = "input.model_source == 'server_path'",
          textInput("model_server_path", NULL,
                    placeholder = "/srv/models/mif-phenotyper-v1.0-merged/")
        ),
        conditionalPanel(
          condition = "input.model_source == 'ollama'",
          textInput("ollama_host", "Ollama host:",
                    value = Sys.getenv("OLLAMA_HOST", "http://host.docker.internal:11434")),
          textInput("ollama_model", "Model tag:",
                    value = "", placeholder = "gemma3:12b-it-qat"),
          tags$small(style = "color:#999;",
            "No fine-tuning needed — the model classifies zero-shot from the tissue type and distinctive markers directly.")
        ),
        actionButton("validate_model", "Validate Model",
                     class = "btn btn-load", icon = icon("check-circle")),
        uiOutput("model_status"),
        conditionalPanel(
          condition = "output.model_is_adapter",
          br(),
          textInput("base_model_path", "Base model (HF ID or local path):",
                    value = "google/gemma-3-1b-it")
        )
      ),

      # ── Step 3: Inference Config ─────────────────────────────────────
      div(class = "step-block",
        div(class = "step-label", span(class = "step-num", "3"), "Inference Config"),
        radioButtons("mode", "Inference mode:",
          choices  = c("Cluster (fast — recommended)" = "cluster",
                       "Per-cell (slow — full resolution)"  = "cell"),
          selected = "cluster"
        ),
        conditionalPanel(
          condition = "input.model_source != 'ollama'",
          checkboxInput("load_in_8bit", "8-bit quantization (halves memory)", value = FALSE),
          br(),
          textAreaInput("ontology", "Cell type ontology:",
                        value = DEFAULT_GEMMA_ONTOLOGY, rows = 6),
          br()
        ),
        conditionalPanel(
          condition = "input.model_source == 'ollama'",
          div(class = "temp-slider",
            sliderInput("ollama_temperature", "Temperature:",
                        min = 0, max = 1, value = 0, step = 0.05, width = "100%")
          ),
          div(class = "temp-scale", span("Deterministic"), span("Varied")),
          tags$small(style = "color:#999; display:block; margin-top:-4px;",
            "0 = greedy decoding: identical input always yields the same label, which is what a reproducible annotation run needs. Raise it only to sample one cluster repeatedly and use the agreement rate across runs as an empirical confidence."),
          br(),
          numericInput("min_z", "Min. elevation to count a marker:",
                       value = 0.5, min = 0, max = 5, step = 0.1),
          tags$small(style = "color:#999; display:block; margin-top:-6px;",
            "A marker must exceed this to be listed as elevated. Clusters with nothing above it are labelled Unclassified rather than guessed."),
          tags$small(style = "color:#999; display:block; margin-top:8px;",
            "Zero-shot: classifies straight from the tissue + distinctive markers below — no ontology table needed."),
          br()
        ),
        textInput("tissue", "Tissue type (optional):", placeholder = "e.g. spleen")
      ),

      # ── Step 4: Execute ──────────────────────────────────────────────
      div(class = "step-block",
        div(class = "step-label", span(class = "step-num", "4"), "Execute"),
        conditionalPanel(
          condition = "input.model_source == 'ollama'",
          checkboxInput("vote_unclassified",
                        "Re-sample Unclassified clusters (10 votes)", value = FALSE),
          tags$small(style = "color:#999; display:block; margin-top:-8px; margin-bottom:8px;",
            "Samples each Unclassified cluster 10x at temperature 0.7 and reports the agreement rate as an empirical confidence. A stable majority replaces the label — except where the min-elevation floor fired, which is only diagnosed, never overridden."),
        ),
        checkboxInput("export_h5ad", "Also export .h5ad (AnnData) file", value = FALSE),
        actionButton("run", "Run Phenotyping",
                     class = "btn btn-run", icon = icon("play")),
        conditionalPanel(
          condition = "output.has_failures",
          br(),
          uiOutput("failure_note"),
          actionButton("retry_failed", "Re-run Failed Clusters",
                       class = "btn btn-retry", icon = icon("rotate-right"))
        ),
        br(), br(),
        downloadButton("download_results", "Download Results",
                       class = "btn btn-dl", icon = icon("file-download")),
        br(), br(),
        actionButton("reset_button", "Reset Page",
                     class = "btn btn-reset", icon = icon("redo"))
      )
    ),

    mainPanel(
      width = 9,
      tabsetPanel(type = "tabs",
        tabPanel("Input Heatmap",
          plotOutput("plots_input", height = "520px")
        ),
        tabPanel("Annotated Heatmap",
          plotOutput("plots_annotated", height = "520px")
        ),
        tabPanel("UMAP & Spatial",
          div(style = "padding:12px 0 6px;",
            selectInput("dr_view", NULL,
              choices  = c("UMAP" = "umap", "Spatial" = "spatial"),
              selected = "umap", width = "200px")
          ),
          plotOutput("plots_dr", height = "540px")
        ),
        tabPanel("Predictions",
          div(style = "padding-top:12px;",
            DT::dataTableOutput("predictions_table")
          )
        )
      )
    )
  )
)

# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------
server <- function(input, output, session) {

  # ── Reactive values ────────────────────────────────────────────────────
  rv <- reactiveValues(
    model_type    = NULL,    # "merged" | "adapter"
    model_dir     = NULL,    # resolved directory path
    model_status  = "idle",  # "idle" | "ok" | "err"
    model_msg     = "",
    model_ontology_manifest = NULL,  # trained markers + ontology, if the checkpoint ships one
    run_complete  = FALSE,
    annotated_spe = NULL,
    # -- retry support: remember what the last run used so failed clusters can
    #    be re-requested without redoing the whole run.
    last_prompts  = NULL,   # path to the prompts JSONL
    last_preds    = NULL,   # path to the predictions JSONL
    last_backend  = NULL,   # "ollama" | "merged" | "adapter"
    n_failed      = 0L
  )

  # ── Session temp dir ───────────────────────────────────────────────────
  tempdir0 <- local({
    d <- glue('{Sys.getenv("PHENOSUITE_TMPDIR", tempdir())}/{session$token}')
    dir.create(d, showWarnings = FALSE, recursive = TRUE)
    as.character(d)
  })
  tracker <- ProvenanceTracker$new("gemma_phenotyper", session, tempdir0)

  session$onSessionEnded(function() {
    unlink(tempdir0, recursive = TRUE)
  })

  # ── Load SPE ───────────────────────────────────────────────────────────
  spe_raw <- reactive({
    req(input$rds_file)
    tracker$register_input(input$rds_file, input_id = "rds_file")
    obj <- tryCatch(readRDS(input$rds_file$datapath), error = function(e) NULL)
    shiny::validate(shiny::need(
      inherits(obj, "SpatialExperiment") || inherits(obj, "SingleCellExperiment"),
      "Uploaded file must be a SpatialExperiment or SingleCellExperiment .rds object."
    ))
    obj
  })

  # ── Expose flags to UI conditionalPanels ──────────────────────────────
  output$has_data <- reactive({ !is.null(input$rds_file) })
  outputOptions(output, "has_data", suspendWhenHidden = FALSE)

  output$model_is_adapter <- reactive({ isTRUE(rv$model_type == "adapter") })
  outputOptions(output, "model_is_adapter", suspendWhenHidden = FALSE)

  output$run_complete <- reactive({ isTRUE(rv$run_complete) })
  outputOptions(output, "run_complete", suspendWhenHidden = FALSE)

  # ── Populate cluster column and marker choices from SPE ───────────────
  observe({
    spe      <- spe_raw()
    manifest <- rv$model_ontology_manifest   # read first so this re-fires once model validation completes
    req(spe)
    cd_cols     <- names(colData(spe))
    cluster_cols <- cd_cols[grepl("cluster", cd_cols, ignore.case = TRUE)]
    if (length(cluster_cols) == 0) cluster_cols <- cd_cols
    markers <- rownames(spe)
    updateSelectInput(session, "cluster_col",
                      choices = cluster_cols, selected = cluster_cols[1])
    updateCheckboxGroupInput(session, "marker_cols",
                             label   = paste0("Markers (", length(markers), "):"),
                             choices = setNames(markers, clean_marker_name(markers)),
                             selected = markers,
                             inline  = FALSE)
    updateTextAreaInput(session, "ontology",
                        value = build_ontology_text(clean_marker_name(markers), manifest))
  })

  # ── Model validation ───────────────────────────────────────────────────
  observeEvent(input$validate_model, {
    src <- input$model_source

    withProgress(message = "Validating model…", value = 0.3, {
      dir_to_check <- NULL
      rv$model_ontology_manifest <- NULL   # reset; repopulated below only on success

      if (src == "ollama") {
        incProgress(0.3, detail = "Pinging Ollama server…")
        res <- validate_ollama(input$ollama_host, input$ollama_model)
        incProgress(0.4, detail = "Done")
        if (!res$ok) {
          rv$model_status <- "err"
          rv$model_msg    <- res$msg
          return()
        }
        rv$model_dir    <- NULL
        rv$model_type   <- "ollama"
        rv$model_status <- "ok"
        rv$model_msg    <- res$msg
        return()
      }

      if (src == "upload") {
        req(input$model_zip)
        zip_path <- input$model_zip$datapath
        out_dir  <- file.path(tempdir0, "model_weights")
        dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

        incProgress(0.3, detail = "Extracting zip…")
        tryCatch(
          utils::unzip(zipfile = zip_path, exdir = out_dir),
          error = function(e) {
            rv$model_status <- "err"
            rv$model_msg    <- paste("Zip extraction failed:", e$message)
            return()
          }
        )
        dir_to_check <- find_model_dir(out_dir)

      } else {
        req(nchar(trimws(input$model_server_path)) > 0)
        dir_to_check <- trimws(input$model_server_path)
      }

      incProgress(0.3, detail = "Checking model type…")

      mtype <- detect_model_type(dir_to_check)

      if (mtype == "not_found") {
        rv$model_status <- "err"
        rv$model_msg    <- paste("Directory not found:", dir_to_check)
        return()
      }
      if (mtype == "unknown") {
        rv$model_status <- "err"
        rv$model_msg    <- "Could not identify model type — expected config.json or adapter_config.json"
        return()
      }

      # Approximate size
      all_f       <- list.files(dir_to_check, full.names = TRUE, recursive = TRUE)
      total_bytes <- sum(file.info(all_f)$size, na.rm = TRUE)
      size_str    <- if (total_bytes > 1e9) sprintf("%.1f GB", total_bytes / 1e9) else
                                             sprintf("%.0f MB", total_bytes / 1e6)

      rv$model_dir    <- dir_to_check
      rv$model_type   <- mtype
      rv$model_status <- "ok"
      rv$model_ontology_manifest <- read_ontology_manifest(dir_to_check)
      rv$model_msg    <- sprintf("%s · %s%s",
        if (mtype == "merged") "Merged checkpoint" else "LoRA adapter", size_str,
        if (is.null(rv$model_ontology_manifest)) "" else " · trained ontology loaded")

      incProgress(0.4, detail = "Done")
    })
  })

  # ── Model status badge ─────────────────────────────────────────────────
  output$model_status <- renderUI({
    cls <- switch(rv$model_status,
      ok  = "key-status key-ok",
      err = "key-status key-err",
      "key-status key-idle"
    )
    lbl <- switch(rv$model_status,
      ok  = paste0("✓ ", rv$model_msg),
      err = paste0("✗ ", rv$model_msg),
      "Click Validate Model after selecting weights"
    )
    tags$span(class = cls, lbl)
  })

  # ── Reset ──────────────────────────────────────────────────────────────
  observeEvent(input$reset_button, { js$resetClick() })

  # ── Input heatmap (available as soon as SPE + cluster col are set) ─────
  output$plots_input <- renderPlot({
    spe <- spe_raw()
    req(spe, input$cluster_col, input$cluster_col %in% names(colData(spe)))
    plot_heatmap.mod(x = spe, group_by = input$cluster_col,
                     out_dir = NULL, size.row = 8, size.col = 8)
  })

  # ── Post-inference finalisation ────────────────────────────────────────
  # Harmonise -> inject -> plot -> write. Factored out of the run handler so
  # the "re-run failed clusters" button can redo exactly these steps after a
  # retry, instead of duplicating them.
  finalize_run <- function(spe, preds, backend_is_ollama, marker_cols, prog = function(...) NULL) {
    # Step 3b – harmonise near-duplicate labels (Ollama only; the
    # fine-tuned checkpoints already emit a controlled vocabulary)
    if (backend_is_ollama && nrow(preds) > 0) {
      prog(0.05, detail = "Harmonising annotation labels…")
      harmonise_map <- tryCatch(
        harmonize_ollama_labels(input$ollama_host, input$ollama_model,
                                preds$cell_type, temperature = input$ollama_temperature),
        error = function(e) { message("Harmonisation failed: ", e$message); NULL }
      )
      if (!is.null(harmonise_map)) {
        hit <- preds$cell_type %in% names(harmonise_map)
        preds$cell_type[hit] <- unname(harmonise_map[preds$cell_type[hit]])
      } else {
        message("Harmonisation step failed or skipped; using raw labels.")
      }
    }

    # Step 4 – inject into SPE
    prog(0.80, detail = "Injecting annotations…")
    if (input$mode == "cluster") {
      spe <- broadcast_cluster_labels(spe, preds, input$cluster_col)
    } else {
      spe[["gemma_cell_type"]]  <- as.character(preds$cell_type)
      spe[["gemma_confidence"]] <- as.numeric(preds$confidence)
    }

    # Step 5 – generate plots + save outputs
    prog(0.85, detail = "Generating plots…")
    generate_colors <- function(n) {
      if (n <= 102) rainbow(n)
      else hsv(seq(0, 1, length.out = n + 1)[1:n], s = 0.8, v = 0.8)
    }

    n_types <- length(unique(spe[["gemma_cell_type"]]))

    plot_heatmap.mod(x = spe, group_by = "gemma_cell_type",
                     out_dir = tempdir0, size.row = 8, size.col = 8)

    tryCatch(
      plot_dr.mod(spe, dr = "UMAP", color_by = "gemma_cell_type",
                  out_dir = tempdir0, h = 20, w = 20),
      error = function(e) message("UMAP plot skipped: ", e$message)
    )
    tryCatch(
      plot_spatial.mod(spe, color_by = "gemma_cell_type",
                       out_dir = tempdir0, h = 20, w = 20,
                       colors  = generate_colors(n_types)),
      error = function(e) message("Spatial plot skipped: ", e$message)
    )

    saveRDS(spe, file.path(tempdir0, "spe_gemma_annotated.rds"))
    write.csv(as.data.frame(colData(spe)),
              file.path(tempdir0, "gemma_colData.csv"))
    # Attach per-cluster signal QC so implausibly-bright clusters are visible in
    # the output rather than only on manual inspection. Reported, not corrected.
    if (input$mode == "cluster" && !is.null(input$cluster_col)) {
      qc <- tryCatch(cluster_signal_qc(spe, marker_cols, input$cluster_col),
                     error = function(e) NULL)
      if (!is.null(qc))
        preds <- merge(preds, qc, by = "cluster_id", all.x = TRUE, sort = FALSE)
    }
    write.csv(preds,
              file.path(tempdir0, "gemma_predictions_raw.csv"), row.names = FALSE)

    # ── Vectra-format CSV export (for PCF-toolkit compatibility) ──
    write_vectra_csv(spe, phenotype_col = "gemma_cell_type",
                      tissue_fallback = input$tissue, out_dir = tempdir0,
                      file_prefix = "vectra_gemma")

    # ── Optional .h5ad (AnnData) export ────────────────────────────
    if (isTRUE(input$export_h5ad)) {
      # Whole object, not a hand-picked subset: AnnData carries colData,
      # spatialCoords (-> obsm["spatial"]) and every assay natively, so there is
      # no reason to flatten it down to a few column attributes as loom required.
      write_h5ad(spe, file.path(tempdir0, "gemma_annotated.h5ad"))
    }


    rv$annotated_spe <- spe
    rv$run_complete  <- TRUE
    invisible(TRUE)
  }

  # ── Main inference pipeline ────────────────────────────────────────────
  observeEvent(input$run, {

    # Guards
    spe <- spe_raw()
    if (is.null(spe)) {
      showNotification("Upload an .rds file first.", type = "error"); return()
    }
    if (length(input$marker_cols) == 0) {
      showNotification("Select at least one marker.", type = "error"); return()
    }
    if (!isTRUE(rv$model_status == "ok") ||
        (rv$model_type != "ollama" && is.null(rv$model_dir))) {
      showNotification("Validate model weights first (Step 2).", type = "error"); return()
    }

    model_dir    <- rv$model_dir
    model_type   <- rv$model_type
    base_model   <- trimws(input$base_model_path %||% "google/gemma-3-1b-it")
    cluster_col  <- if (input$mode == "cluster") input$cluster_col else NULL
    marker_cols  <- input$marker_cols
    ontology     <- trimws(input$ontology)
    load_in_8bit <- isTRUE(input$load_in_8bit)

    rv$run_complete  <- FALSE
    rv$annotated_spe <- NULL

    tracker$capture_parameters(input)
    tracker$analysis_started()

    tmp_in  <- tempfile(fileext = ".jsonl", tmpdir = tempdir0)
    tmp_out <- tempfile(fileext = "_preds.jsonl", tmpdir = tempdir0)

    withProgress(message = "Gemma Phenotyping", value = 0, {

      # Step 1 – serialise prompts
      incProgress(0.05, detail = "Building inference prompts…")
      tryCatch(
        if (model_type == "ollama") {
          prepare_ollama_prompts(
            spe         = spe,
            marker_cols = marker_cols,
            cluster_col = cluster_col,
            tissue      = input$tissue,
            min_z       = if (is.null(input$min_z) || is.na(input$min_z)) 0.5 else input$min_z,
            out_path    = tmp_in
          )
        } else {
          prepare_inference_input(
            spe          = spe,
            marker_cols  = marker_cols,
            ontology_text = ontology,
            cluster_col  = cluster_col,
            out_path     = tmp_in
          )
        },
        error = function(e) {
          showNotification(paste("Prompt build failed:", e$message), type = "error")
          stop(e)
        }
      )

      # Step 2 – run inference
      if (model_type == "ollama") {
        incProgress(0.10, detail = sprintf("Calling %s on %s…", input$ollama_model, input$ollama_host))

        ok <- tryCatch({
          run_ollama_inference(
            tmp_in, tmp_out, input$ollama_host, input$ollama_model,
            temperature = input$ollama_temperature,
            on_progress = function(i, n) {
              step_every <- max(1L, n %/% 20L)
              if (i %% step_every == 0 || i == n)
                incProgress(0.6 * step_every / n, detail = sprintf("Predicting %d / %d…", i, n))
            }
          )
          TRUE
        }, error = function(e) {
          showNotification(paste("Ollama inference failed:", e$message), type = "error", duration = 15)
          FALSE
        })
        if (!isTRUE(ok)) return()

      } else {
        incProgress(0.10, detail = "Loading Gemma model (larger checkpoints take longer)…")

        infer_args <- c(
          INFER_SCRIPT,
          if (model_type == "merged") c("--model_path", model_dir) else
            c("--adapter_path", model_dir, "--base_model", base_model),
          "--input_jsonl",  tmp_in,
          "--output_jsonl", tmp_out
        )
        if (load_in_8bit) infer_args <- c(infer_args, "--load_in_8bit")

        ret    <- system2(GEMMA_PYTHON, args = infer_args,
                          stdout = TRUE, stderr = TRUE, wait = TRUE)
        status <- attr(ret, "status")
        if (!is.null(status) && status != 0L) {
          err_lines <- paste(tail(ret, 30), collapse = "\n")
          showNotification(paste("Inference failed (exit", status, "):\n", err_lines),
                           type = "error", duration = 15)
          return()
        }
        incProgress(0.6, detail = "Inference complete…")
      }

      # Step 3 – parse predictions
      incProgress(0.70, detail = "Parsing predictions…")
      preds <- tryCatch(
        jsonlite::stream_in(file(tmp_out), verbose = FALSE),
        error = function(e) {
          showNotification(paste("Failed to parse predictions:", e$message), type = "error")
          NULL
        }
      )
      req(!is.null(preds))

      # Optional 10-sample vote on Unclassified clusters (Ollama only).
      if (model_type == "ollama" && isTRUE(input$vote_unclassified)) {
        incProgress(0.03, detail = "Re-sampling Unclassified clusters…")
        vres <- tryCatch(
          vote_unclassified_ollama(tmp_in, tmp_out, input$ollama_host, input$ollama_model,
                                   n_samples = 10L, vote_temperature = 0.7),
          error = function(e) { message("Voting failed: ", e$message); NULL })
        if (!is.null(vres) && vres$n > 0) {
          preds <- vres$preds
          showNotification(sprintf("Voting: %d Unclassified cluster(s) re-sampled, %d resolved.",
                                   vres$n, vres$resolved),
                           type = "message", duration = 10)
        }
      }

      # Record state so "re-run failed clusters" can target just the failures.
      rv$last_prompts <- tmp_in
      rv$last_preds   <- tmp_out
      rv$last_backend <- model_type
      rv$n_failed     <- length(failed_prediction_idx(preds))
      if (rv$n_failed > 0)
        showNotification(sprintf("%d cluster(s) failed and were labelled 'unknown'. Use 'Re-run Failed' in Step 4.",
                                 rv$n_failed), type = "warning", duration = 12)

      # finalize_run() stores the *annotated* object in rv itself. Don't assign
      # rv$annotated_spe here: `spe` in this scope is still the un-annotated
      # original (R copies on modify, so the annotation inside finalize_run
      # never touched this binding), and overwriting would strip
      # gemma_cell_type back out from under the plot/table outputs.
      finalize_run(spe, preds, model_type == "ollama", marker_cols,
                   prog = function(...) incProgress(...))

      incProgress(0.15, detail = "Done ✓")
    })

    tracker$analysis_completed()
    showNotification("Phenotyping complete.", type = "message", duration = 5)
  })

  # ── Failure flag + note for the retry control ─────────────────────────
  output$has_failures <- reactive({ isTRUE(rv$n_failed > 0) })
  outputOptions(output, "has_failures", suspendWhenHidden = FALSE)

  output$failure_note <- renderUI({
    req(rv$n_failed > 0)
    tags$span(class = "failure-note",
              sprintf("%d cluster(s) returned no usable label", rv$n_failed))
  })

  # ── Re-run only the clusters that failed ──────────────────────────────
  # Re-requests just those prompts, merges any successes back into the
  # predictions file, then re-runs harmonisation/plots/outputs so every
  # downstream artefact reflects the recovered labels.
  observeEvent(input$retry_failed, {
    if (is.null(rv$last_preds) || !file.exists(rv$last_preds)) {
      showNotification("No previous run to retry.", type = "error"); return()
    }
    if (!identical(rv$last_backend, "ollama")) {
      showNotification(
        "Re-run is only available for the Ollama backend; the fine-tuned backends already retry in-process.",
        type = "warning", duration = 10); return()
    }

    spe <- spe_raw()
    req(spe)

    withProgress(message = "Re-running failed clusters", value = 0.1, {
      res <- tryCatch(
        retry_failed_ollama(
          rv$last_prompts, rv$last_preds,
          input$ollama_host, input$ollama_model,
          temperature = input$ollama_temperature,
          max_attempts = 1L,
          on_progress = function(i, n, attempt)
            incProgress(0.6 / max(n, 1), detail = sprintf("Retrying %d / %d…", i, n))
        ),
        error = function(e) {
          showNotification(paste("Retry failed:", e$message), type = "error", duration = 12)
          NULL
        })
      req(!is.null(res))

      preds <- jsonlite::stream_in(file(rv$last_preds), verbose = FALSE)
      rv$n_failed <- length(failed_prediction_idx(preds))

      incProgress(0.2, detail = "Re-generating outputs…")
      finalize_run(spe, preds, TRUE, input$marker_cols,
                   prog = function(...) incProgress(...))

      showNotification(
        sprintf("Recovered %d of %d failed cluster(s); %d still failing.",
                res$recovered, res$before, res$after),
        type = if (res$after == 0) "message" else "warning", duration = 10)
    })
  })

  # ── Annotated heatmap ──────────────────────────────────────────────────
  output$plots_annotated <- renderPlot({
    req(rv$run_complete, !is.null(rv$annotated_spe))
    plot_heatmap.mod(x = rv$annotated_spe, group_by = "gemma_cell_type",
                     out_dir = NULL, size.row = 8, size.col = 8)
  })

  # ── UMAP / Spatial ─────────────────────────────────────────────────────
  output$plots_dr <- renderPlot({
    req(rv$run_complete, !is.null(rv$annotated_spe))
    spe    <- rv$annotated_spe
    n_types <- length(unique(spe[["gemma_cell_type"]]))
    generate_colors <- function(n) {
      if (n <= 102) rainbow(n)
      else hsv(seq(0, 1, length.out = n + 1)[1:n], s = 0.8, v = 0.8)
    }
    if (is.null(input$dr_view) || input$dr_view == "umap") {
      tryCatch(
        plot_dr.mod(spe, dr = "UMAP", color_by = "gemma_cell_type",
                    out_dir = NULL, h = 20, w = 20),
        error = function(e)
          plot(1, type = "n", main = paste("UMAP not available:", e$message))
      )
    } else {
      tryCatch(
        plot_spatial.mod(spe, color_by = "gemma_cell_type",
                         out_dir = NULL, h = 20, w = 20,
                         colors  = generate_colors(n_types)),
        error = function(e)
          plot(1, type = "n", main = paste("Spatial not available:", e$message))
      )
    }
  })

  # ── Predictions table ──────────────────────────────────────────────────
  output$predictions_table <- DT::renderDataTable({
    req(rv$run_complete, !is.null(rv$annotated_spe))
    cd   <- as.data.frame(colData(rv$annotated_spe))
    keep <- intersect(c(input$cluster_col, "gemma_cluster_id", "gemma_cell_type",
                        "gemma_cell_type_marked", "gemma_confidence", "sample_id"),
                      names(cd))
    DT::datatable(
      cd[, keep, drop = FALSE],
      options = list(pageLength = 25, scrollX = TRUE),
      filter  = "top",
      rownames = FALSE
    )
  })

  # ── Download (zip excludes potentially huge model_weights/ dir) ────────
  output$download_results <- downloadHandler(
    filename    = function() glue("gemma-phenotyping-{Sys.Date()}.zip"),
    content     = function(file) {
      all_files    <- list.files(tempdir0, recursive = TRUE, full.names = FALSE)
      result_files <- all_files[!grepl("^model_weights[/\\\\]", all_files)]
      if (length(result_files) == 0) {
        showNotification("No results to download yet.", type = "warning"); return()
      }
      zip::zip(zipfile = file, files = result_files, root = tempdir0)
    },
    contentType = "application/zip"
  )
}

shinyApp(ui = ui, server = server)
