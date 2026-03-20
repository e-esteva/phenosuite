library(shiny)
library(shinyjs)
library(glue)
library(shinycssloaders)
library(shinyWidgets)
require(tidyverse)

options(shiny.maxRequestSize = 4000 * 1024^2)  # 4 GB upload limit

jsResetCode <- "shinyjs.resetClick = function() { history.go(0); }"

ui <- fluidPage(
  useShinyjs(),
  extendShinyjs(text = jsResetCode, functions = c("resetClick")),

  # ── Custom styles ────────────────────────────────────────────────────
  tags$head(tags$style(HTML("
    body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif; }
    .sidebar-panel  { background: #f8f9fa; border-right: 1px solid #dee2e6; }
    .btn-primary    { background-color: #0dc5c1; border-color: #0bb8b4; }
    .btn-primary:hover { background-color: #0aa8a5; }
    .shiny-notification {
      position: fixed; opacity: 1;
      top: calc(50% - 35px); left: calc(50% - 150px);
      height: 70px; width: 300px;
      border-radius: 8px; box-shadow: 0 4px 12px rgba(0,0,0,.15);
    }
    hr { border-top: 1px solid #dee2e6; }
  "))),

  titlePanel(
    tags$span(
      tags$strong("Masquerade"),
      tags$small(style = "color:#888; margin-left:8px;",
                 "Cell-cluster mask builder")
    )
  ),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      # ── Spatial annotation ───────────────────────────────────────────
      fileInput("spatial_metadata",
                "Spatial annotation CSV",
                accept = c(".csv", "text/csv"),
                placeholder = "cluster + x + y columns"),

      # ── Image source ─────────────────────────────────────────────────
      fileInput("image_source",
                "Multiplex TIFF / qpTIFF (≤ 4 GB)",
                accept = c(".tiff", ".tif", ".qptiff", "image/tiff")),

      # ── Optional marker whitelist ────────────────────────────────────
      selectInput("relevant_markers",
                  "Restrict to marker whitelist?",
                  choices = c("No", "Yes"), selected = "No"),

      conditionalPanel(
        condition = "input.relevant_markers == 'Yes'",
        fileInput("marker_whitelist",
                  "Marker whitelist CSV",
                  accept = c(".csv", "text/csv"))
      ),

      tags$hr(),

      # ── Run controls ─────────────────────────────────────────────────
      textInput("run_label", "Output label", placeholder = "e.g. sample-01"),

      actionButton("run_masquerade", "Render Masks",
                   icon = icon("play"),
                   class = "btn-primary btn-block",
                   style = "margin-top:8px;"),

      tags$hr(),

      # ── Download (appears after processing) ──────────────────────────
      conditionalPanel(
        condition = "output.data_processed == true",
        downloadButton("mask_download_chunked", "Download TIFF",
                       icon = icon("file-download"),
                       class = "btn-block",
                       style = "margin-top:4px;")
      )
    ),

    mainPanel(
      width = 9,
      uiOutput("status_text") %>% withSpinner(color = "#0dc5c1")
    )
  )
)
