renderCircos <- function(logOdds, label, p1, p2, out_dir,
                         continuous_color_scheme = TRUE,
                         scale = FALSE,
                         discontinuity = FALSE,
                         col.fun = NULL,
                         label_size.cex = 0.9,
                         transformed = FALSE,
                         self_interactions = FALSE) {
  require(colorspace)
  require(ComplexHeatmap)
  require(glue)
  require(circlize)

  logOdds <- as.matrix(logOdds)
  if (!self_interactions) {
    diag(logOdds) <- 0
  }

  # --- detect discontinuity (infinities present) ---
  if (sum(is.infinite(logOdds)) > 0) {
    discontinuity <- TRUE
    # Defense in depth: a caller that forgets to stabilise will otherwise
    # hand -Inf into chordDiagramFromMatrix, where `rowSums(abs(mat))`
    # produces Inf, `xlim/sum(xlim)` becomes NaN, `%in% NA` drops every
    # sector, and the whole thing explodes with a cryptic
    # "undefined columns selected" from deep inside [.data.frame.
    # Apply log(exp(x)+1) — same transform circos-artist/builder do — so
    # the matrix is always finite by the time circlize sees it.
    logOdds <- log(exp(logOdds) + 1)
    transformed <- TRUE
  }

  # --- helper: columns that actually carry signal ---
  active_cols <- abs(colSums(logOdds)) > 0
  active_vals <- logOdds[, active_cols]

  # --- build colour function ---
  #
  # Midpoint logic:
  #   - transformed or discontinuity (no natural zero): use data midpoint
  #   - normal:                                          use 0
  #
  # Palette logic:
  #   - continuous_color_scheme:  white -> yellow -> red
  #   - divergent:               blue  -> white  -> red
  #
  if (is.null(col.fun)) {
    # Defense in depth: strip non-finite values before picking ramp endpoints.
    # If a caller forgets to log(exp(x)+1)-stabilise a log-odds matrix, raw
    # -Inf / Inf values would otherwise propagate into colorRamp2 and emit
    # invalid hex strings like "NAFF", which then crash chordDiagramFromMatrix
    # deep inside col2rgb with a misleading "undefined columns selected".
    finite_vals <- active_vals[is.finite(active_vals)]
    if (length(finite_vals) == 0) {
      # Degenerate: nothing finite to interpolate across. Pick a trivial
      # ramp and let the downstream chord diagram handle the empty plot.
      finite_vals <- c(0, 1)
    }
    lo <- min(finite_vals)
    hi <- max(finite_vals)
    if (lo == hi) {
      # Flat matrix — widen slightly so colorRamp2 has a non-zero interval.
      lo <- lo - 0.5
      hi <- hi + 0.5
    }

    # Pick the centre of the colour ramp
    centre <- if (transformed || discontinuity) (lo + hi) / 2 else 0

    if (continuous_color_scheme) {
      col_fun <- colorRamp2(c(lo, centre, hi),
                            c("white", "yellow", "red"))
    } else {
      col_fun <- colorRamp2(c(lo, centre, hi),
                            c("blue", "white", "red"))
    }
  } else {
    col_fun <- col.fun
  }

  # ------------------------------------------------------------------
  #  Internal helper that draws the chord diagram with directional arcs
  # ------------------------------------------------------------------
  .draw_circos <- function(add_title = FALSE) {
    par(cex = label_size.cex, mar = c(1.25, 1.25, 1.25, 1.25))

    chordDiagram(
      logOdds,
      annotationTrack  = "grid",
      preAllocateTracks = list(track.height = 0.1),
      scale            = TRUE,
      col              = col_fun,
      grid.col         = seq(ncol(logOdds)),

      # --- KEY CHANGE: show directionality via diffHeight ---
      # Source end is taller; target (incoming) end is shorter,
      # giving the classic "incoming arc" visual cue.
      directional       = 1,
      direction.type    = c("diffHeight", "arrows"),
      link.arr.type     = "big.arrow",
      diffHeight        = mm_h(3),        # 3 mm height difference
      link.sort         = TRUE,
      link.largest.ontop = TRUE
    )

    if (add_title && !is.null(label) && label != "") {
      title(glue("{label} | {p1};{p2}"))
    }

    # Sector labels (clockwise for narrow sectors, inside for wide ones)
    circos.trackPlotRegion(
      track.index = 1,
      panel.fun = function(x, y) {
        xlim        <- get.cell.meta.data("xlim")
        xplot       <- get.cell.meta.data("xplot")
        ylim        <- get.cell.meta.data("ylim")
        sector.name <- get.cell.meta.data("sector.index")

        if (abs(xplot[2] - xplot[1]) < 20) {
          circos.text(mean(xlim), ylim[1], sector.name,
                      facing = "clockwise", niceFacing = TRUE,
                      adj = c(0, 0.5))
        } else {
          circos.text(mean(xlim), ylim[1], sector.name,
                      facing = "inside", niceFacing = TRUE,
                      adj = c(0.5, 0))
        }
      },
      bg.border = NA
    )

    # Legend
    lgd_links <- Legend(
      at             = round(as.vector(quantile(active_vals)), 4),
      col_fun        = col_fun,
      title_position = "topleft",
      title          = "log-odds"
    )
    lgd_list_vertical <- packLegend(lgd_links)
    draw(lgd_list_vertical,
         x = unit(4, "mm"), y = unit(4, "mm"),
         just = c("left", "bottom"))
  }

  # --- draw to current device ---
  circos.clear()
  .draw_circos(add_title = FALSE)

  # --- save PDF if out_dir supplied ---
  if (!is.null(out_dir)) {
    pdf(glue("{out_dir}/{label}-circos.pdf"))
    circos.clear()
    .draw_circos(add_title = TRUE)
    dev.off()
  }

  circos.clear()
}
