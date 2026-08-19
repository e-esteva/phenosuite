renderCircos <- function(logOdds, label, p1, p2, out_dir,
                         continuous_color_scheme = TRUE,
                         scale = FALSE,
                         discontinuity = FALSE,
                         col.fun = NULL,
                         label_size.cex = 0.9,
                         transformed = FALSE,
                         self_interactions = FALSE,
                         grid.col = NULL,
                         legend_title = "log-odds") {
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
  # Continuous palette:  white (min) -> yellow (median) -> red (max)
  #   Always data-driven so the gradient "heats up" from the lowest
  #   value to the highest regardless of whether the data is 0+ or +/-.
  #
  # Divergent palette:   blue (min) -> white (0) -> red (max)
  #   Centre is 0 for untransformed data (natural divergence point)
  #   or the data midpoint for transformed / 0+ data.
  #
  if (is.null(col.fun)) {
    finite_vals <- active_vals[is.finite(active_vals)]
    if (length(finite_vals) == 0) {
      finite_vals <- c(0, 1)
    }
    lo <- min(finite_vals)
    hi <- max(finite_vals)
    if (lo == hi) {
      lo <- lo - 0.5
      hi <- hi + 0.5
    }

    if (continuous_color_scheme) {
      # Continuous "heat-up": white -> yellow -> red
      # For 0+ (transformed) data: use median of non-zero values as midpoint
      #   so yellow is visible in the gradient (most matrix entries are 0).
      # For +/- data: anchor yellow at 0.
      if (transformed || discontinuity || lo >= 0) {
        nonzero_vals <- finite_vals[finite_vals > 0]
        centre <- if (length(nonzero_vals) > 0) median(nonzero_vals) else (lo + hi) / 2
      } else {
        centre <- 0
      }
      col_fun <- colorRamp2(c(lo, centre, hi),
                            c("white", "yellow", "red"))
    } else {
      # Divergent: blue -> white -> red
      # Use 0 as natural centre for untransformed +/- data;
      # use data midpoint for transformed 0+ data.
      centre <- if (transformed || discontinuity) (lo + hi) / 2 else 0
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

    # Use custom grid colors if provided, otherwise default integer palette
    grid_colors <- if (!is.null(grid.col)) grid.col else seq(ncol(logOdds))

    chordDiagram(
      logOdds,
      annotationTrack  = "grid",
      preAllocateTracks = list(track.height = 0.1),
      scale            = TRUE,
      col              = col_fun,
      grid.col         = grid_colors,

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
      title          = legend_title
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
