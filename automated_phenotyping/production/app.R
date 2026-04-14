require(shiny)
require(stringr)
require(glue)
require(shinyjs)
require(phenomenalist)
require(openai)
require(cowplot)
require(dplyr)
require(lubridate)
require(jsonlite)
source('/srv/shiny-server/phenomenalist/utils/provenance.R')

options(shiny.maxRequestSize = 1000 * 1024^2)
jsResetCode <- "shinyjs.resetClick = function() { history.go(0) }"

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

    /* ── Page header ─────────────────────────────────────────────── */
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

    /* ── Sidebar ─────────────────────────────────────────────────── */
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

    /* ── Form controls ───────────────────────────────────────────── */
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

    .control-label {
      font-size: 12px;
      font-weight: 600;
      color: #444;
      margin-bottom: 4px;
    }

    /* ── Key status badge ────────────────────────────────────────── */
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

    /* ── Buttons ─────────────────────────────────────────────────── */
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

    .btn-reset {
      background: transparent;
      color: #888;
      width: 100%;
      border: 1px solid #d4d0c8;
      font-size: 11px;
    }
    .btn-reset:hover { background: #f0ede6; color: #555; }

    /* ── Model panel (hidden until loaded) ──────────────────────── */
    #model-panel {
      transition: opacity 0.3s;
    }

    /* ── Main panel ──────────────────────────────────────────────── */
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

    /* ── Checkbox group ──────────────────────────────────────────── */
    .checkbox-inline { font-size: 12px; margin-right: 10px; }
    .shiny-input-checkboxgroup > .checkbox-inline { font-size: 12px; }

    /* ── File input ──────────────────────────────────────────────── */
    .btn-file {
      background: #1a1a1a;
      color: #f4f3ef;
      font-size: 11px;
      font-weight: 600;
      border-radius: 2px;
      padding: 3px 10px;
    }

    /* Progress bar ──────────────────────────────────────────────── */
    .progress-bar { background-color: #c8f55a; }
    .progress { border-radius: 2px; }

    /* Shiny notification ────────────────────────────────────────── */
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
    span(class = "subtitle", "Automated Cell Phenotyping")
  ),

  br(),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      # ── Step 1: Data ────────────────────────────────────────────
      div(class = "step-block",
        div(class = "step-label", span(class = "step-num", "1"), "Load Data"),
        fileInput("file1", NULL, accept = ".rds", placeholder = "Choose .rds file…",
                  buttonLabel = tags$span(icon("folder-open"), "Browse")),
        conditionalPanel(
          condition = "output.has_data",
          checkboxGroupInput("rb0", "Clustering resolution:",
                             choices = "", inline = TRUE, selected = NULL)
        )
      ),

      # ── Step 2: API Key ─────────────────────────────────────────
      div(class = "step-block",
        div(class = "step-label", span(class = "step-num", "2"), "OpenAI API Key"),
        passwordInput("API_KEY", NULL, placeholder = "sk-…"),
        actionButton("load_models", "Load Models",
                     class = "btn btn-load", icon = icon("robot")),
        uiOutput("key_status")
      ),

      # ── Step 3: Model & Prompt (shown after key validated) ──────
      div(class = "step-block", id = "model-panel",
        div(class = "step-label", span(class = "step-num", "3"), "Model & Prompt"),
        uiOutput("model_ui"),
        uiOutput("prompt_ui")
      ),

      # ── Step 4: Run ─────────────────────────────────────────────
      div(class = "step-block",
        div(class = "step-label", span(class = "step-num", "4"), "Execute"),
        actionButton("run", "Run Phenotyping",
                     class = "btn btn-run", icon = icon("play")),
        br(), br(),
        downloadButton("phenomenalist_download", "Download Results",
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
          plotOutput("plots_ai", height = "520px")
        ),
        tabPanel("Annotated Heatmap",
          plotOutput("plots_ai3", height = "520px")
        ),
        tabPanel("UMAP & Spatial",
          div(style = "padding: 12px 0 6px;",
            selectInput("dr_view", NULL,
              choices  = c("UMAP" = "umap", "Spatial" = "spatial"),
              selected = "umap", width = "200px")
          ),
          plotOutput("plots_ai2", height = "540px")
        )
      )
    )
  )
)

# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------
server <- shinyServer(function(input, output, session) {

  source('/srv/shiny-server/phenomenalist/utils/RunPhenomenalist-shiny/RunPhenomenalist-shiny.R')
  source('/srv/shiny-server/phenomenalist/utils/RunPhenomenalist-shiny/phenomenalist-utils-shiny.R')
  source('/srv/shiny-server/phenomenalist/utils/plot-scatter.R')

  # ── Reactive values ────────────────────────────────────────────────────
  rv <- reactiveValues(
    models         = NULL,
    key_status     = "idle",   # idle | ok | err
    key_msg        = "",
    run_complete   = FALSE,    # flips TRUE after observeEvent finishes
    annotated_spe  = NULL      # holds the annotated SPE for plots
  )

  # ── Temp dir (per session) ────────────────────────────────────────────
  tempdir0 <- local({
    d <- glue('/apps/home/rtmp/{session$token}')
    dir.create(d, showWarnings = FALSE, recursive = TRUE)
    as.character(d)
  })
  tracker <- ProvenanceTracker$new("automated_phenotyping", session, tempdir0)

  # ── Load SPE ──────────────────────────────────────────────────────────
  mydata0 <- reactive({
    req(input$file1)
    readRDS(input$file1$datapath)
  })

  # ── Flag for conditional panel: data loaded ───────────────────────────
  output$has_data <- reactive({ !is.null(input$file1) })
  outputOptions(output, "has_data", suspendWhenHidden = FALSE)

  # ── Populate clustering resolutions ───────────────────────────────────
  observe({
    spe  <- mydata0()
    req(spe)
    cols <- str_subset(names(colData(spe)), "cluster")
    updateCheckboxGroupInput(session, "rb0", choices = cols, selected = "")
  })

  # ── Load model list when button clicked ───────────────────────────────
  observeEvent(input$load_models, {
    key <- trimws(input$API_KEY)
    if (nchar(key) < 10) {
      rv$key_status <- "err"
      rv$key_msg    <- "Key too short"
      rv$models     <- NULL
      return()
    }

    withProgress(message = "Fetching model list…", value = 0.5, {
      result <- tryCatch({
        cache_path <- '/srv/shiny-server/phenomenalist/utils/gpts/gpt_models-updated.csv'
        use_cache  <- FALSE
        if (file.exists(cache_path)) {
          cached    <- read.csv(cache_path)
          cache_age <- as.numeric(Sys.Date() - max(as.Date(cached$created), na.rm = TRUE))
          use_cache <- cache_age <= 90
        }
        if (use_cache) {
          models_vec <- cached$id
        } else {
          raw        <- openai::list_models(openai_api_key = key)
          models_df  <- raw$data %>%
            filter(owned_by == "system") %>%
            mutate(created = lubridate::as_datetime(created)) %>%
            arrange(desc(created))
          models_vec <- str_subset(models_df$id, "gpt")
          write.csv(models_df[match(models_vec, models_df$id), ],
                    cache_path, row.names = FALSE)
        }
        list(ok = TRUE, models = models_vec)
      }, error = function(e) list(ok = FALSE, msg = conditionMessage(e)))
    })

    if (result$ok && length(result$models) > 0) {
      rv$models     <- result$models
      rv$key_status <- "ok"
      rv$key_msg    <- glue("{length(result$models)} models available")
    } else {
      rv$models     <- NULL
      rv$key_status <- "err"
      rv$key_msg    <- if (!result$ok) result$msg else "No models returned"
    }
  })

  # ── Key status badge ──────────────────────────────────────────────────
  output$key_status <- renderUI({
    cls <- switch(rv$key_status,
      ok  = "key-status key-ok",
      err = "key-status key-err",
      "key-status key-idle"
    )
    lbl <- switch(rv$key_status,
      ok  = paste0("\u2713 ", rv$key_msg),
      err = paste0("\u2717 ", rv$key_msg),
      "Enter key and click Load Models"
    )
    tags$span(class = cls, lbl)
  })

  # ── Model selectInput (only when loaded) ─────────────────────────────
  output$model_ui <- renderUI({
    if (is.null(rv$models)) {
      tags$p(style = "color:#aaa; font-size:12px; font-style:italic;",
             "Model list loads after API key is validated.")
    } else {
      selectInput("model", "GPT model:", choices = rv$models)
    }
  })

  # ── Prompt algorithm UI ───────────────────────────────────────────────
  output$prompt_ui <- renderUI({
    req(rv$models)
    tagList(
      selectInput("prompt_alg", "Prompt algorithm:",
        c("Symmetric single choice" = "v1",
          "Asymmetric 2-choice"     = "v2")),
      conditionalPanel(
        condition = "input.prompt_alg == 'v2'",
        textInput("tissue", "Tissue type:", placeholder = "e.g. spleen")
      )
    )
  })

  # ── Reset ─────────────────────────────────────────────────────────────
  observeEvent(input$reset_button, { js$resetClick() })

  # ── Heatmap: input clusters ───────────────────────────────────────────
  output$plots_ai <- renderPlot({
    spe   <- mydata0()
    group <- input$rb0
    req(spe, length(group) > 0)
    plot_heatmap.mod(x = spe, group_by = group, out_dir = NULL,
                     size.row = 8, size.col = 8)
  })

  # ── Helper: parse v2 GPT response ────────────────────────────────────
  annotations_ <- function(response) {
    nl  <- unlist(strsplit(response, "\n"))
    ci  <- grep("Choice",      nl)
    ei  <- grep("Explanation", nl)
    fmt <- isTRUE(all.equal(ci, ei))
    if (fmt) {
      nl  <- nl[ci]
      pat <- if (length(grep("Choice [0-9]:", nl))) "Choice [0-9]:" else "Choice [A-Z]:"
      cd  <- do.call("rbind", strsplit(nl, pat))[, 2]
      ann <- do.call("rbind", strsplit(cd,  "; Explanation"))[, 1]
      ann <- do.call("rbind", strsplit(ann, "^[ ]"))[, 2]
      exp <- do.call("rbind", strsplit(cd,  "; Explanation"))[, 2]
    } else {
      nl.c <- nl[ci]
      pat  <- if (length(grep("Choice [0-9]:", nl.c))) "Choice [0-9]:" else "Choice [A-Z]:"
      cd   <- do.call("rbind", strsplit(nl.c, pat))[, 2]
      ann  <- do.call("rbind", strsplit(cd,   "^[ ]"))[, 2]
      nl.e <- nl[ei]
      exp  <- do.call("rbind", strsplit(nl.e, "Explanation"))[, 2]
      exp  <- do.call("rbind", strsplit(exp,  "^[ ]"))[, 2]
    }
    list(ann, exp)
  }

  # ── Helper: call GPT with auto temperature fallback ─────────────────
  # use_temp_flag is an environment shared across the session so that
  # once a model rejects temperature=0 we stop trying for all clusters.
  temp_env <- new.env(parent = emptyenv())
  temp_env$use_temp <- TRUE   # TRUE = try 0, FALSE = omit (model default)

  gpt_call <- function(messages_, model_sel, api_key,
                       sys_msg = "You are an expert immunologist.") {
    all_msgs <- c(list(list(role = "system", content = sys_msg)), messages_)
    do_call_ <- function(temp) {
      args <- list(model = model_sel, messages = all_msgs,
                   openai_api_key = api_key)
      if (!is.null(temp)) args$temperature <- temp
      do.call(create_chat_completion, args)
    }

    temp_val <- if (temp_env$use_temp) 0 else NULL
    resp     <- try(do_call_(temp_val))

    if (inherits(resp, "try-error") && temp_env$use_temp &&
        grepl("temperature", conditionMessage(attr(resp, "condition")),
              ignore.case = TRUE)) {
      message(glue("Model {model_sel} does not support temperature=0; switching to default."))
      showNotification(
        glue("Note: {model_sel} does not support temperature=0. Switching to default."),
        type = "warning", duration = 6
      )
      temp_env$use_temp <- FALSE
      resp <- try(do_call_(NULL))
    } else if (inherits(resp, "try-error")) {
      resp <- try(do_call_(temp_val))   # one plain retry
    }
    resp
  }

  # ── Helper: send a single harmonisation/broadening prompt ────────────
  # Returns a named character vector (old_label -> new_label) or NULL on failure.
  gpt_remap <- function(labels, prompt_text, model_sel, api_key) {
    resp <- gpt_call(
      messages_  = list(list(role = "user", content = prompt_text)),
      model_sel  = model_sel,
      api_key    = api_key,
      sys_msg    = "You are an expert immunologist. Return ONLY valid JSON, no markdown, no explanation."
    )
    if (inherits(resp, "try-error")) return(NULL)

    raw <- resp$choices$message.content
    # Strip any accidental ```json ... ``` fencing
    raw <- gsub("^```[[:alpha:]]*\n?|\n?```$", "", trimws(raw))

    mapping <- tryCatch(jsonlite::fromJSON(raw, simplifyVector = TRUE),
                        error = function(e) NULL)
    if (is.null(mapping) || !is.list(mapping) && !is.character(mapping)) return(NULL)

    # fromJSON returns a named list or named character vector
    unlist(mapping)
  }

  # ── MAIN: Phenotyping — lives in observeEvent, fires on button click ──
  observeEvent(input$run, {

    spe   <- mydata0()
    group <- isolate(input$rb0)

    # Guard: validate inputs before doing anything
    if (is.null(spe)) {
      showNotification("Please load an RDS file first.", type = "error"); return()
    }
    if (length(group) == 0) {
      showNotification("Please select a clustering resolution.", type = "error"); return()
    }
    if (is.null(rv$models)) {
      showNotification("Please validate your API key and load models first.", type = "error"); return()
    }
    model_sel   <- isolate(input$model)
    prompt_alg  <- isolate(input$prompt_alg)
    api_key     <- isolate(input$API_KEY)
    tissue_val  <- isolate(input$tissue)

    if (is.null(model_sel) || nchar(trimws(model_sel)) == 0) {
      showNotification("No model selected.", type = "error"); return()
    }

    rv$run_complete <- FALSE
    temp_env$use_temp <- TRUE   # reset per-run

    new_clusters_ <- as.character(spe[[group]])
    new_clusters_ <- gsub("[+]", "pos", new_clusters_)
    new_clusters_ <- gsub("-",   "_",   new_clusters_)

    exprs           <- assays(spe)$exprs
    colnames(exprs) <- new_clusters_

    explanations <- list()
    usage_report <- list()
    clusters_    <- unique(colnames(exprs))

    message("Unique clusters: ", paste(clusters_, collapse = ", "))

    withProgress(message = "Running Phenotyping", value = 0, {

      n_clusters <- length(clusters_)
      # budget: 0.55 for per-cluster calls, 0.10 harmonise, 0.10 broaden,
      #         0.10 plots, 0.15 write/finish
      per_cluster_budget <- 0.55 / n_clusters

      # ── Per-cluster annotation ──────────────────────────────────
      for (idx_i in seq_along(clusters_)) {
        i <- clusters_[[idx_i]]

        incProgress(
          amount  = per_cluster_budget,
          detail  = glue("Annotating cluster {i} ({idx_i}/{n_clusters})…")
        )

        cell_idx         <- colnames(exprs) == i
        cluster_i.avgExp <- if (sum(cell_idx) > 1) rowMeans(exprs[, cell_idx, drop = FALSE]) else exprs[, cell_idx]

        markers  <- do.call("rbind",
                     strsplit(names(cluster_i.avgExp), "_Cytoplasm|_Nucleus_"))[, 1]
        map_vals <- mapply(function(m, v)
          paste0(m, ifelse(v > 0, "+", "-")), markers, cluster_i.avgExp)

        if (prompt_alg == "v1") {
          map_raw  <- sapply(seq_along(cluster_i.avgExp), function(x)
            paste0(c(names(cluster_i.avgExp)[x], cluster_i.avgExp[[x]]), collapse = ":"))
          prompt   <- glue("what celltype is described by {paste0(c(map_raw[cluster_i.avgExp > quantile(cluster_i.avgExp, .95)], map_raw[cluster_i.avgExp < quantile(cluster_i.avgExp, .05)]), collapse = ',')}?")
          content_ <- glue("{prompt} Give me a 3 word response")
        } else {
          prompt   <- glue("what celltype is described by {paste0(c(map_vals[cluster_i.avgExp > quantile(cluster_i.avgExp, .9)], map_vals[cluster_i.avgExp < quantile(cluster_i.avgExp, .05)]), collapse = ',')}?")
          content_ <- glue("{prompt} This is a {tissue_val}. Give me a 3 word response. Give two choices. Explain. Format response like: Choice X: choice; Explanation X: explanation.")
        }

        response <- gpt_call(
          messages_  = list(list(role = "user", content = content_)),
          model_sel  = model_sel,
          api_key    = api_key
        )

        if (!inherits(response, "try-error")) {
          usage_report[[length(usage_report) + 1]] <- data.frame(response$usage)
          msg_content <- response$choices$message.content

          if (prompt_alg != "v1") {
            ann_tmp <- paste0(annotations_(msg_content)[[1]], collapse = ", ")
            ann_tmp <- gsub("-", "_", ann_tmp)
            new_clusters_[new_clusters_ == i] <- glue("{i}-{ann_tmp}")
            explanations[[length(explanations) + 1]] <-
              paste0(annotations_(msg_content)[[2]], collapse = ", ")
          } else {
            ann_tmp <- gsub("-", "_", msg_content)
            new_clusters_[new_clusters_ == i] <- glue("{i}-{ann_tmp}")
          }
        } else {
          message(glue("GPT call failed for cluster {i} after 2 attempts"))
        }
      }

      # ── Extract raw annotation labels (strip cluster-ID prefix) ──
      new_clusters_clean <- do.call("rbind", strsplit(new_clusters_, "[-]"))[, 2]
      unique_raw         <- unique(new_clusters_clean)

      # ── Step 2: Harmonise — collapse near-duplicates ──────────────
      incProgress(0.10, detail = "Harmonising annotation names…")

      harmonise_prompt <- paste0(
        "Below is a list of cell type annotation labels produced by automated immunophenotyping. ",
        "Many entries are near-duplicates that differ only in pluralisation, spacing, capitalisation, ",
        "or minor wording (e.g. 'CD4T cell' and 'CD4T cells', or 'NK cell' and 'Natural Killer cell'). ",
        "Return a JSON object mapping EVERY label in the input list to its canonical standardised form. ",
        "Use Title Case. Do not merge biologically distinct populations. ",
        "Return ONLY the JSON object, no markdown fences, no explanation.\n\n",
        "Labels: ", paste(unique_raw, collapse = ", ")
      )

      harmonise_map <- gpt_remap(unique_raw, harmonise_prompt, model_sel, api_key)

      if (!is.null(harmonise_map)) {
        harmonised_labels <- new_clusters_clean
        for (orig in names(harmonise_map)) {
          harmonised_labels[harmonised_labels == orig] <- harmonise_map[[orig]]
        }
      } else {
        message("Harmonisation step failed; using raw annotations as harmonised.")
        harmonised_labels <- new_clusters_clean
      }

      unique_harmonised <- unique(harmonised_labels)

      # ── Step 3: Broaden — collapse subtypes to lineage ───────────
      incProgress(0.10, detail = "Generating broad lineage annotations…")

      broaden_prompt <- paste0(
        "Below is a list of harmonised cell type labels from an immunofluorescence experiment. ",
        "Map each label to its broad immune lineage category using these rules: ",
        "all T cell subtypes (CD4, CD8, Treg, NKT, gamma-delta, etc.) -> 'T Cell'; ",
        "all B cell subtypes (naive, memory, plasma, plasmablast, etc.) -> 'B Cell'; ",
        "all NK cell variants -> 'NK Cell'; ",
        "all macrophage subtypes (M1, M2, tissue-resident, etc.) -> 'Macrophage'; ",
        "all dendritic cell subtypes (cDC1, cDC2, pDC, moDC, etc.) -> 'Dendritic Cell'; ",
        "all neutrophil variants -> 'Neutrophil'; ",
        "all mast cell variants -> 'Mast Cell'; ",
        "all monocyte subtypes -> 'Monocyte'; ",
        "non-immune stromal, endothelial, or epithelial cells keep their own label. ",
        "If a label is already a broad category, map it to itself. ",
        "Return ONLY a JSON object mapping each input label to its broad category. No markdown, no explanation.\n\n",
        "Labels: ", paste(unique_harmonised, collapse = ", ")
      )

      broaden_map <- gpt_remap(unique_harmonised, broaden_prompt, model_sel, api_key)

      if (!is.null(broaden_map)) {
        broad_labels <- harmonised_labels
        for (orig in names(broaden_map)) {
          broad_labels[broad_labels == orig] <- broaden_map[[orig]]
        }
      } else {
        message("Broadening step failed; skipping broad annotation column.")
        broad_labels <- NULL
      }

      # ── Write usage report ────────────────────────────────────────
      incProgress(0.05, detail = "Writing outputs…")

      usage_report_ <- if (length(usage_report) == 1) {
        unlist(usage_report[[1]])
      } else {
        colSums(do.call("rbind", lapply(usage_report, function(x) unlist(x))))
      }

      annotation.count <- length(grep("annotated_clusters", names(colData(spe))))
      if (annotation.count == 0) {
        spe[["annotated_clusters_marked"]]    <- new_clusters_
        spe[["annotated_clusters"]]           <- new_clusters_clean
        spe[["annotated_clusters_harmonised"]]<- harmonised_labels
        if (!is.null(broad_labels))
          spe[["annotated_clusters_broad"]]   <- broad_labels
        group_out <- "annotated_clusters_harmonised"
      } else {
        v_base                               <- annotation.count + 1
        new.anno                             <- glue("annotated_clusters_v{v_base}")
        spe[[new.anno]]                      <- new_clusters_clean
        spe[[glue("{new.anno}_marked")]]     <- new_clusters_
        spe[[glue("{new.anno}_harmonised")]] <- harmonised_labels
        if (!is.null(broad_labels))
          spe[[glue("{new.anno}_broad")]]    <- broad_labels
        group_out <- glue("{new.anno}_harmonised")
      }

      v <- annotation.count + 1
      write.csv(usage_report_,
                glue("{tempdir0}/usage-report-v{v}.csv"))

      if (length(explanations) > 0)
        write.csv(data.frame(Explanations = unlist(explanations)),
                  glue("{tempdir0}/cluster-explanations-report-v{v}.csv"))

      # Annotation table: raw + harmonised + broad
      ann_tbl <- data.frame(
        Cluster_ID  = unique(spe[[group]]),
        Raw         = do.call("rbind", strsplit(unique(new_clusters_), "[-]"))[, 2],
        Harmonised  = harmonise_map[do.call("rbind", strsplit(unique(new_clusters_), "[-]"))[, 2]],
        stringsAsFactors = FALSE
      )
      if (!is.null(broaden_map))
        ann_tbl$Broad <- broaden_map[ann_tbl$Harmonised]
      write.csv(ann_tbl, glue("{tempdir0}/cluster-annotations-v{v}.csv"),
                row.names = FALSE)

      if (!is.null(harmonise_map))
        write.csv(data.frame(Original   = names(harmonise_map),
                             Harmonised = unname(harmonise_map)),
                  glue("{tempdir0}/harmonisation-map-v{v}.csv"), row.names = FALSE)

      if (!is.null(broaden_map))
        write.csv(data.frame(Harmonised = names(broaden_map),
                             Broad      = unname(broaden_map)),
                  glue("{tempdir0}/broad-map-v{v}.csv"), row.names = FALSE)

      generate_colors <- function(n) {
        if (n <= 102) rainbow(n)
        else hsv(seq(0, 1, length.out = n + 1)[1:n], s = 0.8, v = 0.8)
      }

      incProgress(0.10, detail = "Generating plots…")

      plot_heatmap.mod(x = spe, group_by = group_out, out_dir = tempdir0)
      plot_dr.mod(spe, dr = "UMAP", color_by = group_out,
                  out_dir = tempdir0, h = 20, w = 20)
      plot_spatial.mod(spe, color_by = group_out,
                  out_dir = tempdir0, h = 20, w = 20,
                  colors  = generate_colors(length(unique(spe[[group_out]]))))

      saveRDS(spe, glue("{tempdir0}/spe.rds"))

      spatial_coords         <- data.frame(spatialCoords(spe))
      spatial_coords$cluster <- harmonised_labels
      write.csv(spatial_coords,
                glue("{tempdir0}/annotated_spatial_coords-v{v}.csv"))

      # Store annotated SPE for renderPlots to consume
      rv$annotated_spe  <- spe
      rv$run_complete   <- TRUE

      incProgress(0.10, detail = "Done \u2713")
    })

    showNotification("Phenotyping complete.", type = "message", duration = 5)
  })

  # ── Helper: resolve the latest harmonised annotation column ──────────
  latest_annotation_col <- function(spe) {
    all_cols  <- names(colData(spe))
    # prefer _harmonised suffix; fall back to base annotated_clusters
    harm_cols <- all_cols[grepl("_harmonised$", all_cols)]
    if (length(harm_cols) > 0) {
      # pick highest version number
      versions <- na.omit(as.numeric(sub(".*_v([0-9]+)_harmonised$", "\\1", harm_cols)))
      if (length(versions) == 0) return(harm_cols[1])
      return(glue("annotated_clusters_v{max(versions)}_harmonised"))
    }
    base_cols <- all_cols[grepl("^annotated_clusters", all_cols) & !grepl("_marked|_broad|_harmonised", all_cols)]
    if (length(base_cols) == 0) return(NULL)
    if ("annotated_clusters" %in% base_cols && length(base_cols) == 1)
      return("annotated_clusters")
    versions <- na.omit(as.numeric(sub(".*_v", "", base_cols[grepl("_v[0-9]+$", base_cols)])))
    if (length(versions) == 0) return("annotated_clusters")
    glue("annotated_clusters_v{max(versions)}")
  }

  # ── Heatmap: annotated clusters ───────────────────────────────────────
  output$plots_ai3 <- renderPlot({
    req(rv$run_complete, !is.null(rv$annotated_spe))
    spe       <- rv$annotated_spe
    group_out <- latest_annotation_col(spe)
    req(group_out)
    plot_heatmap.mod(x = spe, group_by = group_out, out_dir = NULL,
                     size.row = 8, size.col = 8)
  })

  # ── UMAP / Spatial: toggled by selectInput ────────────────────────────
  output$plots_ai2 <- renderPlot({
    req(rv$run_complete, !is.null(rv$annotated_spe))
    spe       <- rv$annotated_spe
    group_out <- latest_annotation_col(spe)
    req(group_out)

    generate_colors <- function(n) {
      if (n <= 102) rainbow(n)
      else hsv(seq(0, 1, length.out = n + 1)[1:n], s = 0.8, v = 0.8)
    }

    view <- input$dr_view
    if (is.null(view) || view == "umap") {
      plot_dr.mod(spe, dr = "UMAP", color_by = group_out,
                  out_dir = NULL, h = 20, w = 20)
    } else {
      plot_spatial.mod(spe, color_by = group_out,
                  out_dir = NULL, h = 20, w = 20,
                  colors  = generate_colors(length(unique(spe[[group_out]]))))
    }
  })

  # ── Download ──────────────────────────────────────────────────────────
  output$phenomenalist_download <- downloadHandler(
    filename    = function() glue("phenomenalist-autophenotyping-{Sys.Date()}.zip"),
    content     = function(file) {
      zip::zip(zipfile = file, files = dir(tempdir0), root = tempdir0)
    },
    contentType = "application/zip"
  )
})

shinyApp(ui = ui, server = server)
