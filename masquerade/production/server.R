library(shiny)
library(shinyFiles)
library(stringr)
library(glue)
library(shinyjs)
library(reticulate)
library(plotly)
library(reshape2)
library(dplyr)
library(shinycssloaders)
require(tidyverse)
source('/srv/shiny-server/phenomenalist/utils/provenance.R')

# ── One-time Python environment setup (runs at app start, not per-click) ─────
use_virtualenv("r-reticulate", required = FALSE)
source_python("/srv/shiny-server/phenomenalist/masquerade/production/masquerade.py")
source_python("/srv/shiny-server/phenomenalist/masquerade/production/viewer_utils.py")

server <- function(input, output, session) {

  # ── Reactive values ──────────────────────────────────────────────────
  rv <- reactiveValues(
    data_processed = FALSE,
    outPath        = NULL,
    status_msg     = NULL,
    channel_meta   = NULL,   # list of {name, is_mask, color_r, color_g, color_b}
    img_H          = 0L,
    img_W          = 0L
  )

  # ── File inputs (reactives) ─────────────────────────────────────────
  spatial_metadata <- reactive({
    req(input$spatial_metadata)
    read.csv(input$spatial_metadata$datapath)
  })

  marker_whitelist <- reactive({
    if (is.null(input$marker_whitelist)) return(NULL)
    read.csv(input$marker_whitelist$datapath)
  })

  img_src <- reactive({
    req(input$image_source)
    input$image_source$datapath
  })

  # ── Reset handler ──────────────────────────────────────────────────
  observeEvent(input$reset_button, {
    js$resetClick()
  })

  # ── Expose processed flag to UI conditionalPanel ───────────────────
  output$data_processed <- reactive(rv$data_processed)
  outputOptions(output, "data_processed", suspendWhenHidden = FALSE)

  # ── Session temp directory ─────────────────────────────────────────
  tempdir0 <- file.path(Sys.getenv("PHENOSUITE_TMPDIR", tempdir()), session$token)
  dir.create(tempdir0, showWarnings = FALSE, recursive = TRUE)
  tracker <- ProvenanceTracker$new("masquerade", session, tempdir0)

  # Clean up on session end
  session$onSessionEnded(function() {
    unlink(tempdir0, recursive = TRUE)
  })

  # ── Main processing pipeline ───────────────────────────────────────
  observeEvent(input$run_masquerade, {
    # Validate required inputs
    req(input$spatial_metadata, input$image_source, input$run_label)

    # Provenance: register inputs + capture params
    tracker$register_input(input$image_source, input_id = "image_source")
    tracker$register_input(input$spatial_metadata, input_id = "spatial_metadata")
    if (!is.null(input$marker_whitelist)) {
      tracker$register_input(input$marker_whitelist, input_id = "marker_whitelist")
    }
    tracker$capture_parameters(input)
    tracker$analysis_started()

    rv$data_processed <- FALSE
    rv$outPath <- NULL

    spatial <- spatial_metadata()
    whitelist <- marker_whitelist()
    image_source <- img_src()
    label <- input$run_label

    withProgress(message = "Generating Cluster Masks", value = 0, {

      # Step 1 – Pre-process image
      incProgress(0.05, detail = "Loading & cropping image")
      pre <- PreProcessImage(
        image_source   = image_source,
        spatial_metadata = spatial
      )
      image        <- pre[[1]]
      raw_img_size <- pre[[2]]
      bounds       <- pre[[3]]

      # Step 2 – Build mask channels
      incProgress(0.25, detail = "Generating & compressing mask channels")
      masks <- get_mask_channels(
        image            = image,
        spatial_metadata = spatial,
        raw_img_size     = raw_img_size,
        bounds           = bounds
      )
      channels           <- masks[[1]]
      compression_factor <- masks[[2]]

      # Free the full image immediately
      rm(image, pre, masks)
      gc()

      # Step 3 – Compress biomarker channels
      incProgress(0.50, detail = "Compressing biomarker channels")
      channels <- compress_marker_channels(
        image_source       = image_source,
        channels           = channels,
        compression_factor = compression_factor,
        spatial_metadata   = spatial,
        bounds             = bounds,
        relevant_markers   = whitelist
      )

      # Step 4 – Write output TIFF
      incProgress(0.75, detail = "Writing TIFF")
      out_file <- glue("{tempdir0}/{label}-{round(compression_factor, 2)}x.tiff")
      writeMaskTiff(channels = channels, outPath = out_file)

      rm(channels)
      gc()

      rv$outPath <- out_file

      # Step 5 – Load channels into viewer cache
      incProgress(0.90, detail = "Loading viewer")
      rv$channel_meta <- py$load_tiff_channels(out_file)
      dims            <- py$get_image_dims()
      rv$img_H        <- as.integer(dims[[1]])
      rv$img_W        <- as.integer(dims[[2]])

      rv$data_processed <- TRUE
      incProgress(1.0, detail = "Done")
    })

    # Provenance: write sidecar alongside TIFF in tempdir0
    tracker$analysis_completed()

    showNotification("Masks generated successfully.", type = "message", duration = 4)
  })

  # ── Viewer sidebar controls (rendered after processing) ────────────
  output$viewer_controls <- renderUI({
    req(rv$channel_meta)
    meta <- rv$channel_meta

    channel_rows <- lapply(seq_along(meta), function(i) {
      ch      <- meta[[i]]
      name    <- ch[["name"]]
      safe_id <- gsub("[^A-Za-z0-9]", "_", name)
      hex     <- rgb(ch[["color_r"]], ch[["color_g"]], ch[["color_b"]])

      # Masks hidden by default (toggle on the one(s) you want to isolate);
      # for markers show only the first (typically DAPI)
      default_vis <- if (isTRUE(ch[["is_mask"]])) FALSE else (i == 1L)

      tagList(
        tags$div(
          style = "display:flex; align-items:center; gap:6px; margin-top:6px;",
          tags$div(
            class = "ch-color-picker",
            colourpicker::colourInput(
              inputId    = paste0("ch_color_", safe_id),
              label      = NULL,
              value      = hex,
              showColour = "background"
            )
          ),
          checkboxInput(
            inputId = paste0("ch_vis_", safe_id),
            label   = tags$small(name),
            value   = default_vis
          )
        ),
        tags$div(
          style = "padding-left:20px; margin-bottom:4px;",
          sliderInput(
            inputId = paste0("ch_br_", safe_id),
            label   = NULL,
            min = 0, max = 3, value = 1, step = 0.1,
            ticks = FALSE, width = "100%"
          )
        )
      )
    })

    tagList(
      tags$hr(),
      tags$p(tags$strong("Viewer channels"), style = "margin-bottom:2px; color:#555;"),
      tags$small(style = "color:#999;", "Toggle visibility · adjust brightness"),
      do.call(tagList, channel_rows)
    )
  })

  # ── Composite reactive (recomputes on any channel control change) ───
  composite_b64 <- reactive({
    req(rv$data_processed, rv$channel_meta)
    meta <- rv$channel_meta

    vis_names   <- character(0)
    color_r_vec <- numeric(0)
    color_g_vec <- numeric(0)
    color_b_vec <- numeric(0)
    bright_vec  <- numeric(0)

    for (i in seq_along(meta)) {
      ch      <- meta[[i]]
      name    <- ch[["name"]]
      safe_id <- gsub("[^A-Za-z0-9]", "_", name)

      vis <- input[[paste0("ch_vis_", safe_id)]]
      if (!isTRUE(if (is.null(vis)) TRUE else vis)) next

      br <- input[[paste0("ch_br_", safe_id)]]
      br <- if (is.null(br)) 1.0 else as.numeric(br)

      # Use user-picked colour if available, otherwise fall back to auto-assigned
      color_hex <- input[[paste0("ch_color_", safe_id)]]
      if (!is.null(color_hex) && nchar(trimws(color_hex)) == 7L) {
        cmat <- col2rgb(color_hex) / 255
        cr   <- cmat[1]; cg <- cmat[2]; cb <- cmat[3]
      } else {
        cr <- as.numeric(ch[["color_r"]])
        cg <- as.numeric(ch[["color_g"]])
        cb <- as.numeric(ch[["color_b"]])
      }

      vis_names   <- c(vis_names,   name)
      color_r_vec <- c(color_r_vec, cr)
      color_g_vec <- c(color_g_vec, cg)
      color_b_vec <- c(color_b_vec, cb)
      bright_vec  <- c(bright_vec,  br)
    }

    if (length(vis_names) == 0L) return("")

    py$composite_and_encode(
      visible_names = vis_names,
      color_r       = color_r_vec,
      color_g       = color_g_vec,
      color_b       = color_b_vec,
      brightnesses  = bright_vec
    )
  })

  # ── Composite plotly viewer ─────────────────────────────────────────
  output$composite_viewer <- renderPlotly({
    req(rv$data_processed, rv$img_H > 0)
    b64 <- composite_b64()
    if (nchar(b64) == 0L) return(plotly_empty())

    H <- rv$img_H
    W <- rv$img_W

    plot_ly(
      x    = c(0, W),
      y    = c(0, H),
      type = "scatter",
      mode = "markers",
      marker = list(opacity = 0, size = 1),
      hoverinfo = "none"
    ) %>%
      layout(
        margin = list(l = 0, r = 0, t = 0, b = 0),
        paper_bgcolor = "#111",
        plot_bgcolor  = "#111",
        images = list(list(
          source  = paste0("data:image/png;base64,", b64),
          xref    = "x",   yref    = "y",
          x       = 0,     y       = H,
          sizex   = W,     sizey   = H,
          xanchor = "left", yanchor = "top",
          sizing  = "stretch",
          layer   = "below"
        )),
        xaxis = list(
          range = c(0, W),
          showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE
        ),
        yaxis = list(
          range = c(0, H),
          showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE,
          scaleanchor = "x", scaleratio = 1
        )
      ) %>%
      config(scrollZoom = TRUE, displayModeBar = TRUE,
             modeBarButtonsToRemove = list("lasso2d", "select2d"))
  })

  # ── Status text ─────────────────────────────────────────────────────
  output$status_text <- renderUI({
    if (!rv$data_processed) {
      tags$div(
        style = "padding:40px; text-align:center; color:#888;",
        tags$h4("Upload files and click 'Render Masks' to begin.")
      )
    }
    # Once processed the viewer replaces this; return NULL so the div collapses
  })

  # ── Download handler (zip bundle incl. provenance sidecar) ─────────
  output$mask_download_chunked <- downloadHandler(
    filename = function() {
      req(rv$outPath)
      paste0(tools::file_path_sans_ext(basename(rv$outPath)), ".zip")
    },
    content = function(file) {
      req(rv$outPath)
      source_file <- rv$outPath

      if (!file.exists(source_file)) {
        showNotification("Source file not found.", type = "error", duration = 5)
        stop("Download failed: source file missing.")
      }

      tryCatch({
        showNotification("Preparing download…", type = "message",
                         duration = NULL, id = "dl_note")

        # Ensure provenance sidecar exists alongside the TIFF in tempdir0
        if (is.null(tracker$analysis_end)) tracker$analysis_completed()

        # compression_level 1 (not the package default of 9): mask TIFFs are
        # large (1GB+) but sparse/binary, so low compression captures nearly
        # all the size reduction in a fraction of the time. At level 9 this
        # blocks for minutes with zero bytes reaching the browser, which reads
        # as a hung/failed download long before the zip is ready.
        zip::zip(
          zipfile = file,
          files   = dir(tempdir0),
          root    = tempdir0,
          compression_level = 1
        )

        removeNotification("dl_note")
        showNotification("Download ready!", type = "message", duration = 3)
      }, error = function(e) {
        removeNotification("dl_note")
        showNotification(paste("Download error:", e$message),
                         type = "error", duration = 8)
        # Rethrow so Shiny returns a real failed response instead of quietly
        # serving this error text as if it were the completed .zip.
        stop(e)
      })
    },
    contentType = "application/zip"
  )
}
