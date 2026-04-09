library(shiny)
library(DT)
library(readr)
library(ggplot2)
library(scattermore)

options(shiny.maxRequestSize = 2 * 1024^3)

MAX_RENDER <- 200000

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

detect_delim <- function(path) {
  base <- sub("\\.gz$", "", path, ignore.case = TRUE)
  if (grepl("\\.tsv$", base, ignore.case = TRUE)) "\t" else ","
}

normalise_colnames <- function(nms) {
  nms <- tolower(nms)
  nms <- gsub("\u00b5m", "um", nms)
  nms <- gsub("[^a-z0-9]+", "_", nms)
  nms <- gsub("_+", "_", nms)
  nms <- gsub("^_|_$", "", nms)
  nms
}

resolve_centroids <- function(df) {
  nms      <- names(df)
  nms_norm <- normalise_colnames(nms)
  find_col <- function(p) { idx <- grep(p, nms_norm); if (length(idx)) nms[idx[1]] else NULL }

  if (all(c("XMin","XMax","YMin","YMax") %in% nms)) {
    df$x <- (df$XMin + df$XMax) / 2; df$y <- (df$YMin + df$YMax) / 2
    return(list(df=df, method="Bounding box (XMin/XMax/YMin/YMax)", unit=NA))
  }
  if (all(c("Centroid_X","Centroid_Y") %in% nms)) {
    df$x <- df$Centroid_X; df$y <- df$Centroid_Y
    return(list(df=df, method="Centroid_X / Centroid_Y", unit=NA))
  }
  cx_um <- find_col("^centroid_x_um$"); cy_um <- find_col("^centroid_y_um$")
  cx_px <- find_col("^centroid_x_px$"); cy_px <- find_col("^centroid_y_px$")
  cx    <- find_col("^centroid_x$");    cy    <- find_col("^centroid_y$")
  if (!is.null(cx_um) && !is.null(cy_um)) {
    df$x <- df[[cx_um]]; df$y <- df[[cy_um]]
    return(list(df=df, method=paste("QuPath:", cx_um, "/", cy_um), unit="\u00b5m"))
  }
  if (!is.null(cx_px) && !is.null(cy_px)) {
    df$x <- df[[cx_px]]; df$y <- df[[cy_px]]
    return(list(df=df, method=paste("QuPath:", cx_px, "/", cy_px), unit="px"))
  }
  if (!is.null(cx) && !is.null(cy)) {
    df$x <- df[[cx]]; df$y <- df[[cy]]
    return(list(df=df, method=paste("QuPath:", cx, "/", cy), unit=NA))
  }
  if (all(c("x","y") %in% nms)) return(list(df=df, method="x / y (pre-existing)", unit=NA))
  NULL
}

intensity_cols <- function(df) {
  coord_cols <- c("x","y","XMin","XMax","YMin","YMax","Centroid_X","Centroid_Y")
  nms  <- setdiff(names(df), coord_cols)
  hits <- nms[grepl("intensity|mean", nms, ignore.case=TRUE)]
  hits[sapply(hits, function(n) is.numeric(df[[n]]))]
}

#' Vectorised ray-casting point-in-polygon (arbitrary polygon, not just convex hull)
point_in_polygon <- function(poly_x, poly_y, test_x, test_y) {
  n      <- length(poly_x)
  inside <- rep(FALSE, length(test_x))
  j      <- n
  for (k in seq_len(n)) {
    xi <- poly_x[k]; yi <- poly_y[k]
    xj <- poly_x[j]; yj <- poly_y[j]
    cond   <- ((yi > test_y) != (yj > test_y)) &
              (test_x < (xj - xi) * (test_y - yi) / (yj - yi) + xi)
    inside <- xor(inside, cond)
    j      <- k
  }
  inside
}

#' Recompute selection as union of all finalized gates.
#' Stamps a `gate` column recording which gate(s) each cell belongs to.
apply_gates <- function(gates, df) {
  if (length(gates) == 0) return(NULL)

  # Per-gate membership matrix (nrow(df) x n_gates)
  membership <- vapply(seq_along(gates), function(i) {
    point_in_polygon(gates[[i]]$x, gates[[i]]$y, df$x, df$y)
  }, logical(nrow(df)))

  inside <- rowSums(membership) > 0

  # Comma-separated gate label(s) per cell
  gate_labels <- apply(membership[inside, , drop=FALSE], 1, function(row) {
    paste0("G", which(row), collapse=",")
  })

  out      <- df[inside, ]
  out$gate <- gate_labels
  out
}

# Gate fill colours (one per gate, cycling)
GATE_COLORS <- c("#e74c3c","#3498db","#2ecc71","#f39c12","#9b59b6","#1abc9c")

# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------

ui <- fluidPage(
  titlePanel("Interactive Spatial Gating Tool"),

  tags$head(tags$style(HTML("
    .gate-btn { width:100%; margin-bottom:6px; }
    .info-box  { padding:8px 12px; border-radius:4px; font-size:12px; margin-bottom:8px; }
  "))),

  sidebarLayout(
    sidebarPanel(width = 3,

      fileInput("file", "Choose file",
                accept = c(".csv",".tsv",".csv.gz",".tsv.gz",
                           "application/gzip",
                           "text/csv","text/tab-separated-values")),

      helpText(HTML(
        "Accepted formats:<br>
         &bull; <b>QuPath</b> TSV / TSV.gz<br>
         &bull; <b>HALO</b> CSV (XMin/XMax/YMin/YMax)<br>
         &bull; CSV/TSV with Centroid_X / Centroid_Y"
      )),

      hr(),

      conditionalPanel(condition = "output.fileUploaded",

        uiOutput("coordInfoUI"),
        br(),

        # Channel colour overlay
        uiOutput("channelDropdownUI"),
        hr(),

        # Gate drawing controls
        h4("Gates"),
        p(style="font-size:12px; color:#666;",
          "Click on the plot to place polygon vertices.
           Close the polygon to apply the gate.
           Multiple gates are combined (union)."),

        uiOutput("gateStatusUI"),
        br(),

        actionButton("closeGate", "Close & Apply Gate",
                     class="btn-success gate-btn",
                     icon=icon("check")),

        actionButton("undoGate", "\u21A9  Undo Last Gate",
                     class="btn-warning gate-btn",
                     icon=icon("undo")),

        actionButton("clearGates", "Clear All Gates",
                     class="btn-danger gate-btn",
                     icon=icon("trash")),

        br(),
        uiOutput("selectionInfoUI"),
        br(),

        downloadButton("downloadSelected", "Download Gated Cells",
                       class="btn-primary gate-btn"),
        downloadButton("downloadAll", "Download All Cells",
                       class="btn-default gate-btn"),
        hr(),
        h4("Export Plot"),
        fluidRow(
          column(6, selectInput("imgFormat", label=NULL,
                                choices=c("PNG"="png","PDF"="pdf","SVG"="svg","TIFF"="tiff"),
                                selected="png", width="100%")),
          column(6, numericInput("imgDPI", label=NULL,
                                 value=300, min=72, max=600, step=50, width="100%"))
        ),
        helpText(style="font-size:11px; margin-top:-8px;",
                 "Format | DPI (PNG/TIFF only)"),
        downloadButton("downloadPlot", "Download Plot Image",
                       class="btn-info gate-btn")
      )
    ),

    mainPanel(width = 9,
      conditionalPanel(
        condition = "output.fileUploaded",
        tabsetPanel(
          tabPanel("Spatial Plot",
                   br(),
                   plotOutput("scatterPlot",
                              click   = "plot_click",
                              height  = "650px")),
          tabPanel("All Data",
                   DT::dataTableOutput("originalTable")),
          tabPanel("Gated Data",
                   DT::dataTableOutput("selectedTable"))
        )
      ),
      conditionalPanel(
        condition = "!output.fileUploaded",
        div(style="text-align:center; margin-top:120px;",
            h3("Upload a file to begin"),
            p("QuPath TSV.gz, HALO CSV, or any CSV/TSV with centroid columns"))
      )
    )
  )
)

# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------

server <- function(input, output, session) {

  rv <- reactiveValues(
    raw          = NULL,
    processed    = NULL,
    method       = NULL,
    unit         = NA,
    int_cols     = NULL,
    # Gate state
    gates        = list(),   # list of finalized polygons: each list(x, y, color)
    current_poly = NULL,     # in-progress polygon: data.frame(x, y)
    selected     = NULL,     # gated cells (union of all gates)
    last_plot    = NULL      # most recently rendered ggplot object for export
  )

  # ---- File upload --------------------------------------------------------
  observeEvent(input$file, {
    req(input$file)
    tryCatch({
      delim  <- detect_delim(input$file$name)
      df     <- read_delim(input$file$datapath, delim=delim,
                           show_col_types=FALSE, name_repair="minimal")
      df     <- as.data.frame(df)
      result <- resolve_centroids(df)

      if (is.null(result)) {
        showNotification(
          HTML("Could not find coordinate columns.<br>
               Expected: XMin/XMax/YMin/YMax, Centroid_X/Y, or Centroid X µm/px"),
          type="error", duration=10)
        return()
      }

      rv$raw          <- df
      rv$processed    <- result$df
      rv$method       <- result$method
      rv$unit         <- result$unit
      rv$int_cols     <- intensity_cols(result$df)
      rv$gates        <- list()
      rv$current_poly <- NULL
      rv$selected     <- NULL

      showNotification(
        paste0("Loaded ", format(nrow(df), big.mark=","), " cells. ",
               length(rv$int_cols), " intensity/mean channel(s) found."),
        type="message")
    }, error = function(e) {
      showNotification(paste("Error reading file:", e$message),
                       type="error", duration=10)
    })
  })

  output$fileUploaded <- reactive(!is.null(rv$processed))
  outputOptions(output, "fileUploaded", suspendWhenHidden=FALSE)

  # ---- Coord info banner --------------------------------------------------
  output$coordInfoUI <- renderUI({
    req(rv$method)
    unit_str <- if (!is.na(rv$unit)) paste0(" (", rv$unit, ")") else ""
    div(class="info-box",
        style="background:#1a1a2e; border-left:3px solid #7b61ff;",
        HTML(paste0("<b>Coord source:</b> ", rv$method, unit_str, "<br>",
                    "<b>Total cells:</b> ",
                    format(nrow(rv$processed), big.mark=","))))
  })

  # ---- Channel dropdown ---------------------------------------------------
  output$channelDropdownUI <- renderUI({
    req(rv$int_cols)
    if (length(rv$int_cols) == 0) return(NULL)
    tagList(
      h4("Colour by channel"),
      selectInput("colorChannel", label=NULL,
                  choices  = c("— gate colouring —"="", rv$int_cols),
                  selected = "", width="100%")
    )
  })

  # ---- Gate status banner -------------------------------------------------
  output$gateStatusUI <- renderUI({
    n_vert  <- if (is.null(rv$current_poly)) 0 else nrow(rv$current_poly)
    n_gates <- length(rv$gates)

    if (n_vert == 0 && n_gates == 0) {
      div(class="info-box",
          style="background:#2c2c2c; border-left:3px solid #aaa;",
          HTML("<b>No gate in progress.</b><br>
               Click on the plot to start drawing."))
    } else if (n_vert > 0) {
      gate_col <- GATE_COLORS[[(length(rv$gates) %% length(GATE_COLORS)) + 1]]
      div(class="info-box",
          style=paste0("background:#1a2a1a; border-left:3px solid ", gate_col, ";"),
          HTML(paste0("<b>Drawing gate ", length(rv$gates)+1, "</b><br>",
                      n_vert, " vertices placed<br>",
                      "<span style='color:#aaa; font-size:11px;'>",
                      "Click <b>Close & Apply Gate</b> when done</span>")))
    } else {
      div(class="info-box",
          style="background:#1a2a1a; border-left:3px solid #2ecc71;",
          HTML(paste0("<b>", n_gates, " gate(s) active</b><br>",
                      "<span style='color:#aaa; font-size:11px;'>",
                      "Click plot to draw another gate</span>")))
    }
  })

  # ---- Selection info banner ----------------------------------------------
  output$selectionInfoUI <- renderUI({
    if (is.null(rv$selected)) return(NULL)
    n_all <- nrow(rv$processed); n_sel <- nrow(rv$selected)
    pct   <- round(100 * n_sel / n_all, 1)
    div(class="info-box",
        style="background:#0d2b1a; border-left:3px solid #28a745;",
        HTML(paste0("<b>Gated total:</b> ", format(n_sel, big.mark=","),
                    " / ", format(n_all, big.mark=","),
                    " cells (", pct, "%)")))
  })

  # ---- Plot click: add vertex to current polygon --------------------------
  observeEvent(input$plot_click, {
    req(rv$processed)
    cx <- input$plot_click$x
    cy <- input$plot_click$y
    if (is.null(cx) || is.null(cy)) return()

    new_pt <- data.frame(x=cx, y=cy)
    if (is.null(rv$current_poly)) {
      rv$current_poly <- new_pt
    } else {
      rv$current_poly <- rbind(rv$current_poly, new_pt)
    }
  })

  # ---- Close gate button --------------------------------------------------
  observeEvent(input$closeGate, {
    if (is.null(rv$current_poly) || nrow(rv$current_poly) < 3) {
      showNotification("Need at least 3 vertices to close a gate.", type="warning")
      return()
    }

    gate_col <- GATE_COLORS[[(length(rv$gates) %% length(GATE_COLORS)) + 1]]
    new_gate <- list(x=rv$current_poly$x, y=rv$current_poly$y, color=gate_col)

    rv$gates        <- c(rv$gates, list(new_gate))
    rv$current_poly <- NULL

    withProgress(message="Applying gates\u2026", value=0.5, {
      rv$selected <- apply_gates(rv$gates, rv$processed)
    })

    showNotification(
      paste0("Gate ", length(rv$gates), " applied. ",
             format(nrow(rv$selected), big.mark=","), " cells selected."),
      type="message")
  })

  # ---- Undo last gate -----------------------------------------------------
  observeEvent(input$undoGate, {
    if (!is.null(rv$current_poly)) {
      # If drawing in progress, just clear that first
      rv$current_poly <- NULL
      showNotification("In-progress gate cleared.", type="message")
    } else if (length(rv$gates) > 0) {
      rv$gates <- rv$gates[-length(rv$gates)]
      rv$selected <- if (length(rv$gates) > 0) apply_gates(rv$gates, rv$processed) else NULL
      showNotification(
        paste0("Last gate removed. ", length(rv$gates), " gate(s) remaining."),
        type="message")
    }
  })

  # ---- Clear all gates ----------------------------------------------------
  observeEvent(input$clearGates, {
    rv$gates        <- list()
    rv$current_poly <- NULL
    rv$selected     <- NULL
    showNotification("All gates cleared.", type="message")
  })

  # ---- Main plot ----------------------------------------------------------
  output$scatterPlot <- renderPlot({
    req(rv$processed)
    df    <- rv$processed
    x_lab <- if (!is.na(rv$unit)) paste0("X (", rv$unit, ")") else "X"
    y_lab <- if (!is.na(rv$unit)) paste0("Y (", rv$unit, ")") else "Y"

    # Downsample deterministically for rendering
    if (nrow(df) > MAX_RENDER) {
      set.seed(42)
      df_plot <- df[sort(sample(nrow(df), MAX_RENDER)), ]
    } else {
      df_plot <- df
    }

    # ---- Point colour ----
    ch <- if (!is.null(input$colorChannel)) input$colorChannel else ""

    if (nchar(ch) > 0 && ch %in% names(df_plot)) {
      # Continuous intensity overlay
      vals      <- df_plot[[ch]]
      lo        <- quantile(vals, 0.01, na.rm=TRUE)
      hi        <- quantile(vals, 0.99, na.rm=TRUE)
      df_plot$.color_val <- pmax(pmin(vals, hi), lo)
      pt_aes    <- aes(x=x, y=y, color=.color_val)
      color_scale <- scale_color_gradientn(
        colours = c("#0000FF","#FFFFFF","#FFFF00","#FF0000"),
        name    = ch,
        guide   = guide_colorbar(barwidth=0.8, barheight=8)
      )
      use_color_scale <- TRUE

    } else if (length(rv$gates) > 0) {
      # Gate membership colouring
      inside <- rep(FALSE, nrow(df_plot))
      for (g in rv$gates)
        inside <- inside | point_in_polygon(g$x, g$y, df_plot$x, df_plot$y)
      df_plot$.color_val <- ifelse(inside, "Gated", "Other")
      pt_aes  <- aes(x=x, y=y, color=.color_val)
      color_scale <- scale_color_manual(
        values = c("Gated"="#28a745", "Other"="grey70"),
        name   = NULL,
        guide  = guide_legend(override.aes=list(size=3)))
      use_color_scale <- TRUE

    } else {
      pt_aes          <- aes(x=x, y=y)
      use_color_scale <- FALSE
    }

    # ---- Build plot ----
    p <- ggplot(df_plot, pt_aes)

    # Fast rasterized rendering via scattermore
    p <- p + geom_scattermore(
      pointsize = 1.5,
      alpha     = 0.6,
      pixels    = c(1200, 800)
    )

    if (use_color_scale) p <- p + color_scale

    # ---- Overlay finalized gate polygons ----
    for (i in seq_along(rv$gates)) {
      g    <- rv$gates[[i]]
      gdf  <- data.frame(x=c(g$x, g$x[1]), y=c(g$y, g$y[1]))
      p <- p +
        geom_polygon(data=gdf, aes(x=x, y=y),
                     fill=g$color, alpha=0.08, color=g$color,
                     linewidth=0.8, linetype="solid", inherit.aes=FALSE) +
        annotate("text", x=mean(g$x), y=mean(g$y),
                 label=paste0("G", i), color=g$color,
                 size=4, fontface="bold")
    }

    # ---- Overlay in-progress polygon ----
    if (!is.null(rv$current_poly) && nrow(rv$current_poly) >= 1) {
      gate_col <- GATE_COLORS[[(length(rv$gates) %% length(GATE_COLORS)) + 1]]
      cpoly    <- rv$current_poly

      if (nrow(cpoly) >= 2) {
        p <- p + geom_path(data=cpoly, aes(x=x, y=y),
                           color=gate_col, linewidth=0.9,
                           linetype="dashed", inherit.aes=FALSE)
      }
      p <- p + geom_point(data=cpoly, aes(x=x, y=y),
                           color=gate_col, size=2.5,
                           shape=21, fill="white",
                           stroke=1.2, inherit.aes=FALSE)
    }

    # ---- Theme ----
    p <- p +
      scale_y_reverse() +
      labs(x=x_lab, y=y_lab) +
      theme_minimal(base_size=13) +
      theme(
        panel.background = element_rect(fill="#111111", color=NA),
        plot.background  = element_rect(fill="#111111", color=NA),
        panel.grid.major = element_line(color="#2a2a2a"),
        panel.grid.minor = element_blank(),
        axis.text        = element_text(color="#aaaaaa"),
        axis.title       = element_text(color="#cccccc"),
        legend.background = element_rect(fill="#111111", color=NA),
        legend.text      = element_text(color="#cccccc"),
        legend.title     = element_text(color="#cccccc")
      )

    rv$last_plot <- p
    p
  }, bg="transparent")

  # ---- Tables -------------------------------------------------------------
  output$originalTable <- DT::renderDataTable({
    req(rv$raw)
    DT::datatable(rv$raw,
                  options=list(scrollX=TRUE, pageLength=10),
                  caption=paste("All data —", format(nrow(rv$raw), big.mark=","), "rows"))
  })

  output$selectedTable <- DT::renderDataTable({
    if (is.null(rv$selected)) {
      return(DT::datatable(
        data.frame(Message="No gates applied yet."),
        options=list(dom="t"), rownames=FALSE))
    }
    DT::datatable(rv$selected,
                  options=list(scrollX=TRUE, pageLength=10),
                  caption=paste("Gated —", format(nrow(rv$selected), big.mark=","), "rows"))
  })

  # ---- Downloads ----------------------------------------------------------
  strip_xy <- function(df) {
    orig_nms <- names(rv$raw)
    internal <- c("x", "y", ".color_val")
    to_drop  <- setdiff(internal, orig_nms)  # keep gate column always
    df[, setdiff(names(df), to_drop), drop=FALSE]
  }

  out_ext <- reactive({
    req(input$file)
    if (grepl("\\.tsv", input$file$name, ignore.case=TRUE)) ".tsv" else ".csv"
  })
  write_out <- function(df, file) {
    if (out_ext() == ".tsv") write_tsv(df, file) else write_csv(df, file)
  }

  output$downloadSelected <- downloadHandler(
    filename = function() paste0("gated_cells_", Sys.Date(), out_ext()),
    content  = function(file) { req(rv$selected); write_out(strip_xy(rv$selected), file) }
  )
  output$downloadAll <- downloadHandler(
    filename = function() paste0("all_cells_", Sys.Date(), out_ext()),
    content  = function(file) { req(rv$raw); write_out(rv$raw, file) }
  )

  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste0("spatial_gates_", Sys.Date(), ".", input$imgFormat)
    },
    content = function(file) {
      req(rv$last_plot)
      fmt <- input$imgFormat
      dpi <- as.integer(input$imgDPI)
      if (fmt == "pdf") {
        ggsave(file, plot=rv$last_plot, device="pdf",
               width=12, height=8, units="in")
      } else if (fmt == "svg") {
        ggsave(file, plot=rv$last_plot, device="svg",
               width=12, height=8, units="in")
      } else if (fmt == "tiff") {
        ggsave(file, plot=rv$last_plot, device="tiff",
               width=12, height=8, units="in", dpi=dpi,
               compression="lzw")
      } else {
        # PNG default
        ggsave(file, plot=rv$last_plot, device="png",
               width=12, height=8, units="in", dpi=dpi,
               bg="#111111")
      }
    }
  )
}

shinyApp(ui=ui, server=server)
