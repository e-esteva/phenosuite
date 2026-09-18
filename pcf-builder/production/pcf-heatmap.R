# pcf-heatmap.R
# ─────────────────────────────────────────────────────────────────────────────
# Mean-interaction heatmap for the PCF builder: one row per available cell
# type, one column per labelled sample group, each cell the mean normalised
# PCF of that cell type's interactions with the reference within that group.
#
# Kept free of any Shiny dependency so it can be sourced and tested on its
# own; server.R supplies the already-renamed data frame.
# ─────────────────────────────────────────────────────────────────────────────

require(ggplot2)

# global_pcf: rows = observations, one numeric column per cell type, plus a
#             `Sample` column holding the (renamed) group label.
# Returns a cell type x group matrix of means, NA-safe.
pcf_heatmap_matrix <- function(global_pcf, celltypes, group_order = NULL) {
  groups <- unique(as.character(global_pcf$Sample))
  if (!is.null(group_order) && setequal(group_order, groups)) groups <- group_order

  m <- vapply(groups, function(g) {
    rows <- as.character(global_pcf$Sample) == g
    vapply(celltypes, function(ct) {
      v <- suppressWarnings(as.numeric(global_pcf[rows, ct]))
      if (all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE)
    }, numeric(1))
  }, numeric(length(celltypes)))

  m <- matrix(m, nrow = length(celltypes),
              dimnames = list(celltypes, groups))
  m
}

# Rows ordered by similarity so related cell types sit together; 'All' is the
# pooled baseline rather than a cell type, so it is pinned to the top instead
# of being clustered in among them.
pcf_heatmap_row_order <- function(m) {
  has_all <- "All" %in% rownames(m)
  body <- m[setdiff(rownames(m), "All"), , drop = FALSE]
  ord <- rownames(body)
  if (nrow(body) >= 3) {
    complete <- body[stats::complete.cases(body), , drop = FALSE]
    if (nrow(complete) >= 3) {
      d <- stats::dist(complete)
      if (all(is.finite(d))) {
        hc <- stats::hclust(d, method = "average")
        ord <- c(rownames(complete)[hc$order],
                 setdiff(rownames(body), rownames(complete)))
      }
    }
  }
  c(if (has_all) "All", ord)
}

# `cap_quantile` squishes the colour scale at an upper quantile: normalised PCF
# is heavily right-skewed, and one extreme cell would otherwise flatten every
# other tile to the same colour. Values above the cap keep their label and take
# the end colour.
pcf_heatmap_plot <- function(m,
                             title        = "Mean PCF interactions",
                             subtitle     = NULL,
                             show_values  = TRUE,
                             digits       = 2,
                             cap_quantile = 0.98,
                             midpoint     = 1) {
  df <- expand.grid(celltype = rownames(m), sample = colnames(m),
                    stringsAsFactors = FALSE)
  df$value <- as.vector(m)
  df$celltype <- factor(df$celltype, levels = rev(rownames(m)))
  df$sample   <- factor(df$sample,   levels = colnames(m))

  finite <- df$value[is.finite(df$value)]
  hi  <- if (length(finite)) stats::quantile(finite, cap_quantile, names = FALSE) else 1
  lo  <- if (length(finite)) min(finite) else 0
  if (!is.finite(hi) || hi <= lo) hi <- lo + 1
  mid <- min(max(midpoint, lo), hi)

  p <- ggplot(df, aes(x = sample, y = celltype, fill = value)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    scale_fill_gradient2(
      low = "#2c7bb6", mid = "#f7f7f7", high = "#d7191c",
      midpoint = mid, limits = c(lo, hi), oob = scales::squish,
      na.value = "grey92", name = "mean\nnorm PCF"
    ) +
    scale_x_discrete(position = "top", expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL,
         caption = paste0("1 = random expectation; colour capped at the ",
                          round(cap_quantile * 100), "th percentile",
                          "  |  grey = no interaction measured")) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid       = element_blank(),
      axis.text.x.top  = element_text(angle = 45, hjust = 0, face = "bold",
                                    margin = margin(b = 4)),
      axis.text.y      = element_text(face = "bold"),
      axis.ticks       = element_blank(),
      plot.title       = element_text(face = "bold", size = 15),
      # Generous gap: the column labels are angled and drawn above the panel,
    # so a tight subtitle margin lets them run into the subtitle.
    plot.subtitle    = element_text(colour = "grey35", margin = margin(b = 26)),
      plot.caption     = element_text(colour = "grey45", hjust = 0, size = 9),
      legend.key.height = unit(1.4, "lines"),
      plot.margin      = margin(12, 18, 12, 12)
    )

  if (show_values) {
    lab <- formatC(df$value, format = "f", digits = digits)
    lab[!is.finite(df$value)] <- ""
    # Dark text on pale tiles, white on saturated ones. Intensity on a
    # diverging scale is distance from the *midpoint*, not from the scale
    # minimum — measuring from `lo` paints near-midpoint (i.e. near-white)
    # tiles as if they were saturated and makes their labels invisible.
    span   <- max(mid - lo, hi - mid)
    strong <- is.finite(df$value) & span > 0 &
      (abs(pmin(pmax(df$value, lo), hi) - mid) / span >= 0.72)
    p <- p + geom_text(aes(label = lab),
                       colour = ifelse(strong, "white", "grey15"),
                       size = 3.1, fontface = "bold")
  }
  p
}

# Sizing that keeps cells legible regardless of matrix shape. Width has to
# clear three things, not just the columns: the longest row label, the legend,
# and the subtitle/caption strip — a narrow matrix (2 groups) otherwise
# renders a figure too narrow for its own caption, which then clips.
pcf_heatmap_size <- function(m) {
  label_in <- 0.075 * max(nchar(rownames(m)), 0)
  list(width  = max(9.5, 1.15 * ncol(m) + label_in + 4.2),
       height = max(4.5, 0.42 * nrow(m) + 2.4))
}
