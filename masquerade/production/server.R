library(shiny)
library(shinyFiles)
library(stringr)
library(glue)
library(shinyjs)
library(reticulate)
library(reshape2)
library(dplyr)
library(shinycssloaders)
require(tidyverse)
source('/srv/shiny-server/phenomenalist/utils/provenance.R')

# ── One-time Python environment setup (runs at app start, not per-click) ─────
use_virtualenv("r-reticulate", required = FALSE)
source_python("/srv/shiny-server/phenomenalist/masquerade/production/masquerade.py")

server <- function(input, output, session) {

  # ── Reactive values ──────────────────────────────────────────────────
  rv <- reactiveValues(
    data_processed = FALSE,
    outPath        = NULL,
    status_msg     = NULL
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
  tempdir0 <- file.path("/apps/home/rtmp", session$token)
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
      rv$data_processed <- TRUE
      incProgress(1.0, detail = "Done")
    })

    showNotification("Masks generated successfully.", type = "message", duration = 4)
  })

  # ── Status text (optional, replaces the dummy plotOutput) ──────────
  output$status_text <- renderUI({
    if (rv$data_processed) {
      tags$div(
        style = "padding:40px; text-align:center; color:#0dc5c1;",
        tags$h4(icon("check-circle"), " Processing complete"),
        tags$p("Use the download button to retrieve your TIFF.")
      )
    } else {
      tags$div(
        style = "padding:40px; text-align:center; color:#888;",
        tags$h4("Upload files and click 'Render Masks' to begin.")
      )
    }
  })

  # ── Download handler (chunked copy) ────────────────────────────────
  output$mask_download_chunked <- downloadHandler(
    filename = function() {
      req(rv$outPath)
      basename(rv$outPath)
    },
    content = function(file) {
      req(rv$outPath)
      source_file <- rv$outPath

      if (!file.exists(source_file)) {
        showNotification("Source file not found.", type = "error", duration = 5)
        writeLines("Download failed: source file missing.", file)
        return()
      }

      tryCatch({
        showNotification("Preparing download…", type = "message",
                         duration = NULL, id = "dl_note")

        file.copy(source_file, file, overwrite = TRUE)

        if (file.size(file) != file.size(source_file)) {
          stop("File copy incomplete.")
        }

        removeNotification("dl_note")
        showNotification("Download ready!", type = "message", duration = 3)
      }, error = function(e) {
        removeNotification("dl_note")
        showNotification(paste("Download error:", e$message),
                         type = "error", duration = 8)
        writeLines(paste("Download failed:", e$message), file)
      })
    },
    contentType = "image/tiff"
  )
}
