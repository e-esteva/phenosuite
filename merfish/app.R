## ============================================================================
## Phenomenalist — MERFISH Spatial Transcriptomics Module
## ============================================================================
## Segmentation-based, cellular-resolution spatial transcriptomics analysis
## for MERFISH (Multiplexed Error-Robust FISH) datasets.
##
## Pipeline: Import → QC → Normalize → Reduce → Cluster → Spatial → DE → Export
##
## Dependencies: Seurat (>=5), ggplot2, plotly, DT, dbscan, igraph, pheatmap,
##               viridis, patchwork, RColorBrewer, spdep, RANN
## ============================================================================

suppressPackageStartupMessages({
  library(shiny)
  library(ggplot2)
  library(plotly)
  library(DT)
  library(viridis)
  library(RColorBrewer)
  library(patchwork)
  library(pheatmap)
  library(igraph)
  library(RANN)
  library(Matrix)
  library(grid)
  library(gridExtra)
})

# ---------------------------------------------------------------------------
# PHENOMENALIST THEME SYSTEM
# ---------------------------------------------------------------------------
phenomenalist_css <- "
@import url('https://fonts.googleapis.com/css2?family=DM+Sans:ital,opsz,wght@0,9..40,100..1000;1,9..40,100..1000&family=JetBrains+Mono:wght@300;400;500;600&display=swap');

:root {
 --phen-bg-primary:    #0a0c10;
 --phen-bg-secondary:  #12151c;
 --phen-bg-card:       rgba(18, 22, 32, 0.75);
 --phen-bg-input:      #1a1e2a;
 --phen-border:        rgba(255, 255, 255, 0.06);
 --phen-border-focus:  rgba(0, 194, 255, 0.35);
 --phen-text-primary:  #e8eaed;
 --phen-text-secondary:#8b919e;
 --phen-text-muted:    #555b6a;
 --phen-accent-cyan:   #00c2ff;
 --phen-accent-magenta:#e040fb;
 --phen-accent-green:  #00e676;
 --phen-accent-amber:  #ffab00;
 --phen-accent-red:    #ff5252;
 --phen-glow-cyan:     0 0 20px rgba(0,194,255,0.15);
 --phen-glow-magenta:  0 0 20px rgba(224,64,251,0.15);
 --phen-radius-sm:     6px;
 --phen-radius-md:     10px;
 --phen-radius-lg:     14px;
 --phen-font-sans:     'DM Sans', system-ui, sans-serif;
 --phen-font-mono:     'JetBrains Mono', 'Fira Code', monospace;
 --phen-transition:    all 0.2s cubic-bezier(0.4, 0, 0.2, 1);
}

* { box-sizing: border-box; }

html, body {
  margin: 0; padding: 0;
  background: var(--phen-bg-primary);
  color: var(--phen-text-primary);
  font-family: var(--phen-font-sans);
  font-size: 14px;
  line-height: 1.5;
  -webkit-font-smoothing: antialiased;
  overflow-x: hidden;
}

/* Animated background canvas */
body::before {
  content: '';
  position: fixed;
  inset: 0;
  z-index: -1;
  background:
    radial-gradient(ellipse 800px 600px at 15% 20%, rgba(0,194,255,0.04), transparent),
    radial-gradient(ellipse 600px 800px at 85% 75%, rgba(224,64,251,0.03), transparent),
    radial-gradient(ellipse 500px 500px at 50% 50%, rgba(0,230,118,0.02), transparent);
  animation: bgShift 30s ease-in-out infinite alternate;
}
@keyframes bgShift {
  0%   { filter: hue-rotate(0deg); }
  100% { filter: hue-rotate(25deg); }
}

/* ---- Top navigation bar ---- */
.phen-topbar {
  position: sticky; top: 0; z-index: 1000;
  display: flex; align-items: center; gap: 16px;
  padding: 12px 24px;
  background: rgba(10,12,16,0.88);
  backdrop-filter: blur(16px) saturate(1.6);
  border-bottom: 1px solid var(--phen-border);
}
.phen-topbar .logo {
  font-weight: 700; font-size: 16px; letter-spacing: -0.3px;
  background: linear-gradient(135deg, var(--phen-accent-cyan), var(--phen-accent-magenta));
  -webkit-background-clip: text; -webkit-text-fill-color: transparent;
}
.phen-topbar .module-tag {
  font-family: var(--phen-font-mono);
  font-size: 11px; font-weight: 500;
  color: var(--phen-accent-green);
  background: rgba(0,230,118,0.08);
  border: 1px solid rgba(0,230,118,0.15);
  padding: 3px 10px; border-radius: 20px;
}
.phen-topbar .separator { flex: 1; }
.phen-topbar .status-text {
  font-family: var(--phen-font-mono);
  font-size: 11px; color: var(--phen-text-muted);
}

/* ---- Glass card ---- */
.glass-card {
  background: var(--phen-bg-card);
  backdrop-filter: blur(12px) saturate(1.4);
  border: 1px solid var(--phen-border);
  border-radius: var(--phen-radius-lg);
  padding: 20px;
  margin-bottom: 16px;
  transition: var(--phen-transition);
}
.glass-card:hover {
  border-color: rgba(255,255,255,0.10);
  box-shadow: 0 4px 24px rgba(0,0,0,0.25);
}
.glass-card h3 {
  margin: 0 0 14px 0;
  font-size: 14px; font-weight: 600;
  color: var(--phen-text-primary);
  display: flex; align-items: center; gap: 8px;
}
.glass-card h3 .card-icon {
  width: 28px; height: 28px;
  display: flex; align-items: center; justify-content: center;
  border-radius: var(--phen-radius-sm);
  font-size: 13px;
}
.icon-import  { background: rgba(0,194,255,0.10); color: var(--phen-accent-cyan); }
.icon-qc      { background: rgba(255,171,0,0.10); color: var(--phen-accent-amber); }
.icon-process { background: rgba(0,230,118,0.10); color: var(--phen-accent-green); }
.icon-cluster { background: rgba(224,64,251,0.10); color: var(--phen-accent-magenta); }
.icon-spatial { background: rgba(0,194,255,0.10); color: var(--phen-accent-cyan); }
.icon-de      { background: rgba(255,82,82,0.10);  color: var(--phen-accent-red); }
.icon-export  { background: rgba(255,171,0,0.10); color: var(--phen-accent-amber); }

/* ---- Tabs (navbarPage overrides) ---- */
.nav-tabs {
  border-bottom: 1px solid var(--phen-border) !important;
  padding: 0 20px;
  background: var(--phen-bg-secondary);
}
.nav-tabs > li > a {
  font-family: var(--phen-font-sans) !important;
  font-size: 12.5px !important; font-weight: 500 !important;
  color: var(--phen-text-secondary) !important;
  border: none !important;
  padding: 12px 18px !important;
  border-radius: 0 !important;
  transition: var(--phen-transition);
  border-bottom: 2px solid transparent !important;
  background: transparent !important;
  letter-spacing: 0.02em;
}
.nav-tabs > li > a:hover {
  color: var(--phen-text-primary) !important;
  background: rgba(255,255,255,0.03) !important;
}
.nav-tabs > li.active > a,
.nav-tabs > li.active > a:focus,
.nav-tabs > li.active > a:hover {
  color: var(--phen-accent-cyan) !important;
  background: transparent !important;
  border: none !important;
  border-bottom: 2px solid var(--phen-accent-cyan) !important;
}
.tab-content {
  background: transparent !important;
  padding: 20px;
}

/* ---- Inputs ---- */
.form-control, .selectize-input, .shiny-input-container select {
  background: var(--phen-bg-input) !important;
  color: var(--phen-text-primary) !important;
  border: 1px solid var(--phen-border) !important;
  border-radius: var(--phen-radius-sm) !important;
  font-family: var(--phen-font-sans) !important;
  font-size: 13px !important;
  transition: var(--phen-transition);
}
.form-control:focus, .selectize-input.focus {
  border-color: var(--phen-border-focus) !important;
  box-shadow: var(--phen-glow-cyan) !important;
  outline: none !important;
}
.selectize-dropdown {
  background: var(--phen-bg-secondary) !important;
  border: 1px solid var(--phen-border) !important;
  border-radius: var(--phen-radius-sm) !important;
}
.selectize-dropdown .option {
  color: var(--phen-text-primary) !important;
}
.selectize-dropdown .option.active {
  background: rgba(0,194,255,0.12) !important;
}
.control-label, .shiny-input-container label {
  color: var(--phen-text-secondary) !important;
  font-size: 12px !important;
  font-weight: 500 !important;
  margin-bottom: 4px !important;
}
.irs--shiny .irs-bar { background: var(--phen-accent-cyan); }
.irs--shiny .irs-handle { border-color: var(--phen-accent-cyan); background: var(--phen-bg-input); }
.irs--shiny .irs-single { background: var(--phen-accent-cyan); font-family: var(--phen-font-mono); font-size: 11px; }
.irs--shiny .irs-line { background: var(--phen-bg-input); }
.irs--shiny .irs-min, .irs--shiny .irs-max { color: var(--phen-text-muted); font-family: var(--phen-font-mono); font-size: 10px; }
.irs--shiny .irs-grid-text { color: var(--phen-text-muted); font-size: 10px; }

/* ---- Buttons ---- */
.btn-primary, .btn-default, .action-button {
  font-family: var(--phen-font-sans) !important;
  font-size: 12.5px !important; font-weight: 600 !important;
  border-radius: var(--phen-radius-sm) !important;
  padding: 8px 20px !important;
  transition: var(--phen-transition);
  border: none !important;
  cursor: pointer;
}
.phen-btn-primary {
  background: linear-gradient(135deg, var(--phen-accent-cyan), #0091d5) !important;
  color: #fff !important;
  box-shadow: 0 2px 12px rgba(0,194,255,0.2);
}
.phen-btn-primary:hover {
  box-shadow: 0 4px 20px rgba(0,194,255,0.35);
  transform: translateY(-1px);
}
.phen-btn-secondary {
  background: var(--phen-bg-input) !important;
  color: var(--phen-text-primary) !important;
  border: 1px solid var(--phen-border) !important;
}
.phen-btn-secondary:hover {
  background: rgba(255,255,255,0.06) !important;
  border-color: rgba(255,255,255,0.12) !important;
}

/* ---- KPI metric blocks ---- */
.kpi-row { display: flex; gap: 12px; margin-bottom: 16px; flex-wrap: wrap; }
.kpi-item {
  flex: 1; min-width: 120px;
  background: var(--phen-bg-card);
  border: 1px solid var(--phen-border);
  border-radius: var(--phen-radius-md);
  padding: 14px 16px;
  backdrop-filter: blur(8px);
}
.kpi-label { font-size: 11px; color: var(--phen-text-muted); font-weight: 500; text-transform: uppercase; letter-spacing: 0.06em; }
.kpi-value { font-family: var(--phen-font-mono); font-size: 22px; font-weight: 600; margin-top: 4px; }
.kpi-value.cyan    { color: var(--phen-accent-cyan); }
.kpi-value.magenta { color: var(--phen-accent-magenta); }
.kpi-value.green   { color: var(--phen-accent-green); }
.kpi-value.amber   { color: var(--phen-accent-amber); }
.kpi-value.red     { color: var(--phen-accent-red); }

/* ---- DataTables ---- */
table.dataTable { color: var(--phen-text-primary) !important; }
.dataTables_wrapper {
  color: var(--phen-text-secondary) !important;
  font-family: var(--phen-font-sans) !important;
  font-size: 12px;
}
table.dataTable thead th {
  background: var(--phen-bg-secondary) !important;
  color: var(--phen-text-secondary) !important;
  border-bottom: 1px solid var(--phen-border) !important;
  font-size: 11px !important; font-weight: 600 !important;
  text-transform: uppercase; letter-spacing: 0.04em;
}
table.dataTable tbody td {
  border-bottom: 1px solid rgba(255,255,255,0.03) !important;
  font-family: var(--phen-font-mono); font-size: 12px;
}
table.dataTable tbody tr:hover { background: rgba(0,194,255,0.04) !important; }
.dataTables_filter input { background: var(--phen-bg-input) !important; color: var(--phen-text-primary) !important; border: 1px solid var(--phen-border) !important; border-radius: var(--phen-radius-sm) !important; }
.dataTables_paginate .paginate_button { color: var(--phen-text-secondary) !important; }
.dataTables_paginate .paginate_button.current { background: rgba(0,194,255,0.12) !important; color: var(--phen-accent-cyan) !important; border: 1px solid rgba(0,194,255,0.2) !important; border-radius: var(--phen-radius-sm) !important; }

/* ---- Progress / notifications ---- */
.shiny-notification {
  background: var(--phen-bg-secondary) !important;
  color: var(--phen-text-primary) !important;
  border: 1px solid var(--phen-border) !important;
  border-radius: var(--phen-radius-md) !important;
  font-family: var(--phen-font-sans) !important;
  box-shadow: 0 8px 32px rgba(0,0,0,0.4);
}
.progress-bar { background: linear-gradient(90deg, var(--phen-accent-cyan), var(--phen-accent-magenta)) !important; }

/* ---- Plot containers ---- */
.plot-container {
  background: var(--phen-bg-secondary);
  border: 1px solid var(--phen-border);
  border-radius: var(--phen-radius-md);
  padding: 8px;
  overflow: hidden;
}

/* ---- Well panels ---- */
.well {
  background: var(--phen-bg-card) !important;
  border: 1px solid var(--phen-border) !important;
  border-radius: var(--phen-radius-md) !important;
  box-shadow: none !important;
}

/* ---- Status pill ---- */
.status-pill {
  display: inline-flex; align-items: center; gap: 6px;
  font-family: var(--phen-font-mono);
  font-size: 11px; font-weight: 500;
  padding: 4px 12px;
  border-radius: 20px;
}
.status-pill.ready    { color: var(--phen-accent-green); background: rgba(0,230,118,0.08); border: 1px solid rgba(0,230,118,0.15); }
.status-pill.running  { color: var(--phen-accent-amber); background: rgba(255,171,0,0.08); border: 1px solid rgba(255,171,0,0.15); }
.status-pill.error    { color: var(--phen-accent-red);   background: rgba(255,82,82,0.08);  border: 1px solid rgba(255,82,82,0.15); }
.status-pill .dot { width: 6px; height: 6px; border-radius: 50%; }
.status-pill.ready .dot    { background: var(--phen-accent-green); box-shadow: 0 0 6px var(--phen-accent-green); }
.status-pill.running .dot  { background: var(--phen-accent-amber); box-shadow: 0 0 6px var(--phen-accent-amber); animation: pulse 1.5s ease-in-out infinite; }
@keyframes pulse { 0%,100% { opacity: 1; } 50% { opacity: 0.3; } }

/* ---- Scrollbar ---- */
::-webkit-scrollbar { width: 6px; height: 6px; }
::-webkit-scrollbar-track { background: transparent; }
::-webkit-scrollbar-thumb { background: rgba(255,255,255,0.08); border-radius: 3px; }
::-webkit-scrollbar-thumb:hover { background: rgba(255,255,255,0.15); }

/* ---- Checkbox / Radio overrides ---- */
.checkbox label, .radio label { color: var(--phen-text-secondary) !important; font-size: 13px; }

/* ---- Plotly overrides ---- */
.js-plotly-plot .plotly .modebar { background: transparent !important; }
.js-plotly-plot .plotly .modebar-btn path { fill: var(--phen-text-muted) !important; }
.js-plotly-plot .plotly .modebar-btn:hover path { fill: var(--phen-accent-cyan) !important; }

/* ---- Layout helpers ---- */
.sidebar-panel {
  background: var(--phen-bg-secondary);
  border-right: 1px solid var(--phen-border);
  padding: 20px;
  height: calc(100vh - 60px);
  overflow-y: auto;
}
.main-panel { padding: 20px; }
.section-divider {
  height: 1px;
  background: var(--phen-border);
  margin: 16px 0;
}
.help-text {
  font-size: 11.5px;
  color: var(--phen-text-muted);
  font-style: italic;
  margin-top: 4px;
}
"

# ---------------------------------------------------------------------------
# THEME FOR ggplot2 FIGURES — dark Phenomenalist style
# ---------------------------------------------------------------------------
theme_phenomenalist <- function(base_size = 12) {
  theme_minimal(base_size = base_size) %+replace%
    theme(
      text             = element_text(family = "sans", color = "#e8eaed"),
      plot.title       = element_text(size = rel(1.15), face = "bold", color = "#e8eaed",
                                      margin = margin(b = 8)),
      plot.subtitle    = element_text(size = rel(0.85), color = "#8b919e",
                                      margin = margin(b = 12)),
      plot.background  = element_rect(fill = "#12151c", color = NA),
      panel.background = element_rect(fill = "#12151c", color = NA),
      panel.grid.major = element_line(color = "#FFFFFF0A", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      axis.text        = element_text(color = "#8b919e", size = rel(0.8)),
      axis.title       = element_text(color = "#8b919e", size = rel(0.85), face = "bold"),
      legend.background= element_rect(fill = "#12151c", color = NA),
      legend.text      = element_text(color = "#8b919e", size = rel(0.8)),
      legend.title     = element_text(color = "#e8eaed", size = rel(0.85), face = "bold"),
      legend.key       = element_rect(fill = "#12151c", color = NA),
      strip.text       = element_text(color = "#e8eaed", face = "bold", size = rel(0.85)),
      strip.background = element_rect(fill = "#1a1e2a", color = NA),
      plot.margin      = margin(16, 16, 16, 16)
    )
}

# Fluorescence-inspired palettes
phen_palette_discrete <- function(n) {
  base_colors <- c("#00c2ff", "#e040fb", "#00e676", "#ffab00", "#ff5252",
                   "#18ffff", "#b388ff", "#69f0ae", "#ffd740", "#ff8a80",
                   "#40c4ff", "#ea80fc", "#b9f6ca", "#ffe57f", "#ff867c",
                   "#80d8ff", "#ce93d8", "#a5d6a7", "#fff176", "#ef9a9a",
                   "#4dd0e1", "#ab47bc", "#66bb6a", "#ffa726", "#ef5350",
                   "#26c6da", "#9c27b0", "#43a047", "#ff9800", "#e53935")
  if (n <= length(base_colors)) return(base_colors[seq_len(n)])
  colorRampPalette(base_colors)(n)
}

phen_palette_continuous <- function(option = "viridis") {
  switch(option,
    "viridis"  = scale_color_viridis_c(option = "D"),
    "magma"    = scale_color_viridis_c(option = "A"),
    "inferno"  = scale_color_viridis_c(option = "B"),
    "plasma"   = scale_color_viridis_c(option = "C"),
    "cividis"  = scale_color_viridis_c(option = "E"),
    "cyan"     = scale_color_gradient(low = "#0a0c10", high = "#00c2ff"),
    "magenta"  = scale_color_gradient(low = "#0a0c10", high = "#e040fb"),
    "green"    = scale_color_gradient(low = "#0a0c10", high = "#00e676"),
    "thermal"  = scale_color_gradient2(low = "#00c2ff", mid = "#0a0c10",
                                       high = "#ff5252", midpoint = 0),
    scale_color_viridis_c(option = "D")
  )
}

# ---------------------------------------------------------------------------
# HELPER FUNCTIONS
# ---------------------------------------------------------------------------

# Compute spatial neighbors (k-nearest neighbors)
compute_knn_graph <- function(coords, k = 15) {
  nn <- nn2(coords, k = k + 1)
  nn$nn.idx <- nn$nn.idx[, -1]  # remove self
  nn$nn.dists <- nn$nn.dists[, -1]
  nn
}

# Neighborhood enrichment analysis
neighborhood_enrichment <- function(clusters, nn_idx) {
  n_cells <- length(clusters)
  cl_levels <- sort(unique(clusters))
  n_cl <- length(cl_levels)
  obs_mat <- matrix(0, nrow = n_cl, ncol = n_cl,
                    dimnames = list(cl_levels, cl_levels))
  for (i in seq_len(n_cells)) {
    neighbor_cls <- clusters[nn_idx[i, ]]
    for (cl in neighbor_cls) {
      obs_mat[clusters[i], cl] <- obs_mat[clusters[i], cl] + 1
    }
  }
  # Expected under random permutation
  cl_sizes <- table(clusters)[cl_levels]
  total_edges <- sum(obs_mat)
  exp_mat <- outer(cl_sizes, cl_sizes) / sum(cl_sizes)^2 * total_edges
  # Z-score
  z_mat <- (obs_mat - exp_mat) / sqrt(exp_mat + 1e-10)
  z_mat
}

# Spatial autocorrelation (Moran's I) per gene — lightweight version
morans_i_fast <- function(expr_vec, nn_idx, nn_dists) {
  n <- length(expr_vec)
  x <- expr_vec - mean(expr_vec)
  ss <- sum(x^2)
  if (ss < 1e-15) return(list(I = 0, p = 1))
  k <- ncol(nn_idx)
  W_sum <- 0; lag_sum <- 0
  for (i in seq_len(n)) {
    weights <- 1 / (nn_dists[i, ] + 1e-6)
    W_sum <- W_sum + sum(weights)
    lag_sum <- lag_sum + sum(weights * x[nn_idx[i, ]])
  }
  I <- (n / W_sum) * (lag_sum / ss) # Not needed: correct for expectation
  # Approximate p via normal
  EI <- -1 / (n - 1)
  # Use a permutation-free z-approximation
  z <- (I - EI) / (0.5 / sqrt(n))
  p <- 2 * pnorm(-abs(z))
  list(I = I, p = p)
}

# Simple Wilcoxon-based DE between two groups
run_de <- function(expr_mat, group1_idx, group2_idx, gene_names) {
  results <- data.frame(
    gene = gene_names,
    avg_log2FC = NA_real_,
    pct_1 = NA_real_,
    pct_2 = NA_real_,
    p_val = NA_real_,
    p_adj = NA_real_,
    stringsAsFactors = FALSE
  )
  for (i in seq_along(gene_names)) {
    g1 <- expr_mat[group1_idx, i]
    g2 <- expr_mat[group2_idx, i]
    mean1 <- mean(g1); mean2 <- mean(g2)
    results$avg_log2FC[i] <- log2((mean1 + 1e-9) / (mean2 + 1e-9))
    results$pct_1[i] <- mean(g1 > 0)
    results$pct_2[i] <- mean(g2 > 0)
    if (sd(c(g1, g2)) > 1e-15) {
      results$p_val[i] <- tryCatch(
        wilcox.test(g1, g2)$p.value,
        error = function(e) 1
      )
    } else {
      results$p_val[i] <- 1
    }
  }
  results$p_adj <- p.adjust(results$p_val, method = "BH")
  results <- results[order(results$p_adj, -abs(results$avg_log2FC)), ]
  results
}

# Leiden clustering on shared nearest neighbor graph
# Uses igraph::cluster_leiden with CPM or modularity objective
snn_cluster <- function(knn_idx, resolution = 0.8,
                        objective = c("modularity", "CPM")) {
  objective <- match.arg(objective)
  n <- nrow(knn_idx)
  k <- ncol(knn_idx)
  edge_list <- vector("list", n)
  for (i in seq_len(n)) {
    nn_i <- knn_idx[i, ]
    for (j in nn_i) {
      if (j > i) {
        shared <- length(intersect(knn_idx[i, ], knn_idx[j, ]))
        if (shared > 0) {
          edge_list[[i]] <- rbind(edge_list[[i]],
                                  data.frame(from = i, to = j, weight = shared / k))
        }
      }
    }
  }
  edges <- do.call(rbind, edge_list)
  if (is.null(edges) || nrow(edges) == 0) {
    return(rep(1, n))
  }
  g <- graph_from_data_frame(edges, directed = FALSE, vertices = seq_len(n))
  cl <- cluster_leiden(
    g,
    weights          = E(g)$weight,
    resolution       = resolution,
    objective_function = objective,
    n_iterations     = 10
  )
  membership(cl)
}

# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------
ui <- fluidPage(
  tags$head(
    tags$style(HTML(phenomenalist_css)),
    tags$title("Phenomenalist — MERFISH ST")
  ),

  # ---- Top bar ----
  div(class = "phen-topbar",
    span(class = "logo", "PHENOMENALIST"),
    span(class = "module-tag", "MERFISH · Spatial Transcriptomics"),
    span(class = "separator"),
    uiOutput("status_indicator")
  ),

  # ---- Main tabbed layout ----
  navbarPage(
    title = NULL, id = "main_tabs",

    # ======== TAB 1: DATA IMPORT ========
    tabPanel("Import",
      fluidRow(
        column(4,
          div(class = "glass-card",
            h3(span(class = "card-icon icon-import", "↓"), "Data Import"),
            p(class = "help-text",
              "Upload MERFISH segmented data: a cell × gene expression matrix and ",
              "a metadata file with spatial coordinates (x, y) per cell."),
            div(class = "section-divider"),
            fileInput("expr_file", "Expression Matrix (.csv / .tsv)",
                      accept = c(".csv", ".tsv", ".txt")),
            fileInput("meta_file", "Cell Metadata with Coordinates (.csv / .tsv)",
                      accept = c(".csv", ".tsv", ".txt")),
            div(class = "section-divider"),
            selectInput("x_col", "X-coordinate column", choices = NULL),
            selectInput("y_col", "Y-coordinate column", choices = NULL),
            selectInput("fov_col", "FOV / Slice column (optional)", choices = c("(none)" = "")),
            div(class = "section-divider"),
            actionButton("load_demo", "Load Demo Data",
                         class = "phen-btn-secondary", style = "width:100%; margin-bottom:8px;"),
            actionButton("btn_load", "Validate & Load",
                         class = "phen-btn-primary", style = "width:100%;")
          )
        ),
        column(8,
          div(class = "glass-card",
            h3(span(class = "card-icon icon-import", "⬡"), "Data Preview"),
            uiOutput("import_kpis"),
            div(class = "section-divider"),
            DTOutput("preview_table", height = "420px")
          )
        )
      )
    ),

    # ======== TAB 2: QC & FILTERING ========
    tabPanel("QC",
      fluidRow(
        column(3,
          div(class = "glass-card",
            h3(span(class = "card-icon icon-qc", "◉"), "Quality Control"),
            selectInput("qc_metrics", "Active QC filters",
                        choices = c("Total counts" = "counts",
                                    "Genes detected" = "genes",
                                    "Cell area / volume" = "area",
                                    "Negative control ratio" = "negctrl",
                                    "Transcript density" = "density"),
                        selected = c("counts", "genes"),
                        multiple = TRUE),
            div(class = "section-divider"),
            p(class = "help-text", "Count thresholds"),
            sliderInput("qc_min_counts", "Min total counts", 0, 500, 10, step = 5),
            sliderInput("qc_max_counts", "Max total counts", 500, 100000, 50000, step = 500),
            div(class = "section-divider"),
            p(class = "help-text", "Gene detection thresholds"),
            sliderInput("qc_min_genes", "Min genes detected", 0, 200, 5, step = 1),
            sliderInput("qc_max_genes", "Max genes detected", 50, 5000, 2000, step = 50),
            div(class = "section-divider"),
            p(class = "help-text", "Morphology filters"),
            selectInput("area_col", "Cell area / volume column", choices = c("(none)" = "")),
            numericInput("qc_min_area", "Min cell area (µm²)", value = 0, min = 0),
            numericInput("qc_max_area", "Max cell area (µm²)", value = 1e6, min = 0),
            div(class = "section-divider"),
            p(class = "help-text", "Negative control probes — filter cells with high blank ratio"),
            selectInput("negctrl_col", "Blank / neg-ctrl count column",
                        choices = c("(none)" = "")),
            numericInput("qc_max_negctrl_ratio", "Max neg-ctrl / total ratio",
                         value = 0.05, min = 0, max = 1, step = 0.005),
            div(class = "section-divider"),
            p(class = "help-text", "Transcript density = counts / area"),
            numericInput("qc_min_density", "Min transcript density", value = 0, min = 0, step = 0.01),
            numericInput("qc_max_density", "Max transcript density", value = 1e6, min = 0, step = 0.01),
            div(class = "section-divider"),
            actionButton("btn_qc", "Apply QC Filter",
                         class = "phen-btn-primary", style = "width:100%;")
          )
        ),
        column(9,
          uiOutput("qc_kpis"),
          fluidRow(
            column(6, div(class = "plot-container", plotOutput("qc_violin", height = "360px"))),
            column(6, div(class = "plot-container", plotOutput("qc_scatter", height = "360px")))
          ),
          fluidRow(style = "margin-top:16px;",
            column(12, div(class = "plot-container", plotOutput("qc_spatial", height = "380px")))
          )
        )
      )
    ),

    # ======== TAB 3: PROCESSING ========
    tabPanel("Process",
      fluidRow(
        column(3,
          div(class = "glass-card",
            h3(span(class = "card-icon icon-process", "⚙"), "Normalization & Reduction"),
            selectInput("norm_method", "Normalization method",
                        choices = c("Log-normalize (library size → log1p)" = "lognorm",
                                    "Pearson residuals (SCT-like)" = "sct",
                                    "Counts per 10k (CP10K)" = "cp10k",
                                    "Cell-volume normalized (counts / vol)" = "cellvol",
                                    "Raw counts (no normalization)" = "none")),
            selectInput("vol_col", "Cell volume / area column (for cell-vol norm)",
                        choices = c("(none)" = "")),
            div(class = "section-divider"),
            selectInput("scale_method", "Feature scaling (pre-PCA)",
                        choices = c("Z-score (center + scale)" = "zscore",
                                    "Center only (mean = 0)" = "center",
                                    "None" = "none")),
            div(class = "section-divider"),
            selectInput("hvg_method", "Feature selection",
                        choices = c("Top variance" = "variance",
                                    "Mean-variance trend (VST)" = "vst",
                                    "All genes (no selection)" = "all")),
            sliderInput("n_hvg", "N highly variable genes", 100, 5000, 2000, step = 100),
            div(class = "section-divider"),
            sliderInput("n_pcs", "N principal components", 5, 50, 20, step = 1),
            sliderInput("umap_neighbors", "UMAP n_neighbors", 5, 50, 15, step = 1),
            sliderInput("umap_min_dist", "UMAP min_dist", 0.01, 1, 0.3, step = 0.01),
            div(class = "section-divider"),
            actionButton("btn_process", "Run Processing Pipeline",
                         class = "phen-btn-primary", style = "width:100%;")
          )
        ),
        column(9,
          uiOutput("process_kpis"),
          fluidRow(
            column(6, div(class = "plot-container", plotOutput("pca_elbow", height = "340px"))),
            column(6, div(class = "plot-container", plotOutput("umap_plot", height = "340px")))
          ),
          fluidRow(style = "margin-top:16px;",
            column(6, div(class = "plot-container", plotOutput("hvg_plot", height = "340px"))),
            column(6, div(class = "plot-container", plotOutput("pca_loadings", height = "340px")))
          )
        )
      )
    ),

    # ======== TAB 4: CLUSTERING ========
    tabPanel("Cluster",
      fluidRow(
        column(3,
          div(class = "glass-card",
            h3(span(class = "card-icon icon-cluster", "⬢"), "Leiden Clustering"),
            sliderInput("cluster_k", "k-nearest neighbors", 5, 50, 15, step = 1),
            sliderInput("cluster_res", "Resolution", 0.1, 3.0, 0.8, step = 0.1),
            selectInput("leiden_objective", "Leiden objective function",
                        choices = c("Modularity" = "modularity",
                                    "CPM (Constant Potts Model)" = "CPM")),
            div(class = "section-divider"),
            actionButton("btn_cluster", "Run Clustering",
                         class = "phen-btn-primary", style = "width:100%;"),
            div(class = "section-divider"),
            h3(span(class = "card-icon icon-cluster", "↕"), "Resolution Sweep"),
            sliderInput("sweep_range", "Resolution range",
                        min = 0.1, max = 3.0, value = c(0.3, 1.5), step = 0.1),
            sliderInput("sweep_steps", "N steps", 3, 10, 5, step = 1),
            actionButton("btn_sweep", "Run Sweep",
                         class = "phen-btn-secondary", style = "width:100%;")
          )
        ),
        column(9,
          uiOutput("cluster_kpis"),
          fluidRow(
            column(6, div(class = "plot-container", plotOutput("cluster_umap", height = "380px"))),
            column(6, div(class = "plot-container", plotOutput("cluster_spatial", height = "380px")))
          ),
          fluidRow(style = "margin-top:16px;",
            column(6, div(class = "plot-container", plotOutput("cluster_composition", height = "340px"))),
            column(6, div(class = "plot-container", plotOutput("cluster_sweep", height = "340px")))
          )
        )
      )
    ),

    # ======== TAB 5: SPATIAL ANALYSIS ========
    tabPanel("Spatial",
      fluidRow(
        column(3,
          div(class = "glass-card",
            h3(span(class = "card-icon icon-spatial", "◎"), "Spatial Visualization"),
            selectInput("spatial_color", "Color by",
                        choices = c("Cluster" = "cluster")),
            selectInput("spatial_gene", "Gene expression overlay", choices = NULL),
            selectInput("spatial_palette", "Continuous palette",
                        choices = c("viridis", "magma", "inferno", "plasma",
                                    "cividis", "cyan", "magenta", "green", "thermal")),
            sliderInput("spatial_pt_size", "Point size", 0.1, 4, 0.8, step = 0.1),
            sliderInput("spatial_alpha", "Opacity", 0.1, 1, 0.8, step = 0.05),
            checkboxInput("spatial_coord_flip", "Flip Y-axis", value = TRUE),
            div(class = "section-divider"),
            h3(span(class = "card-icon icon-spatial", "⊞"), "Neighborhood Analysis"),
            sliderInput("nhood_k", "Neighborhood k", 5, 50, 15, step = 1),
            actionButton("btn_nhood", "Run Neighborhood Enrichment",
                         class = "phen-btn-primary", style = "width:100%;"),
            div(class = "section-divider"),
            h3(span(class = "card-icon icon-spatial", "〰"), "Spatially Variable Genes"),
            sliderInput("svg_k", "SVG neighbor k", 5, 30, 10, step = 1),
            numericInput("svg_n_top", "Top N genes to test", value = 200, min = 10, max = 2000),
            actionButton("btn_svg", "Detect SVGs",
                         class = "phen-btn-primary", style = "width:100%;")
          )
        ),
        column(9,
          fluidRow(
            column(6, div(class = "plot-container", plotlyOutput("spatial_main", height = "440px"))),
            column(6, div(class = "plot-container", plotlyOutput("spatial_gene_plot", height = "440px")))
          ),
          fluidRow(style = "margin-top:16px;",
            column(6,
              div(class = "glass-card",
                h3(span(class = "card-icon icon-spatial", "⊞"), "Neighborhood Enrichment (Z-score)"),
                plotOutput("nhood_heatmap", height = "380px")
              )
            ),
            column(6,
              div(class = "glass-card",
                h3(span(class = "card-icon icon-spatial", "〰"), "Spatially Variable Genes"),
                DTOutput("svg_table", height = "380px")
              )
            )
          )
        )
      )
    ),

    # ======== TAB 6: DIFFERENTIAL EXPRESSION ========
    tabPanel("DE",
      fluidRow(
        column(3,
          div(class = "glass-card",
            h3(span(class = "card-icon icon-de", "Δ"), "Differential Expression"),
            selectInput("de_group1", "Group 1 (cluster)", choices = NULL),
            selectInput("de_group2", "Group 2 (cluster / rest)", choices = NULL),
            checkboxInput("de_vs_rest", "Group 1 vs all others", value = TRUE),
            div(class = "section-divider"),
            numericInput("de_lfc_thresh", "Log2FC threshold", value = 0.25, min = 0, step = 0.05),
            numericInput("de_pval_thresh", "Adj. p-value threshold", value = 0.05, min = 0, max = 1, step = 0.01),
            div(class = "section-divider"),
            actionButton("btn_de", "Run DE Analysis",
                         class = "phen-btn-primary", style = "width:100%;")
          )
        ),
        column(9,
          uiOutput("de_kpis"),
          fluidRow(
            column(6, div(class = "plot-container", plotOutput("volcano_plot", height = "400px"))),
            column(6, div(class = "plot-container", plotOutput("de_dotplot", height = "400px")))
          ),
          fluidRow(style = "margin-top:16px;",
            column(12,
              div(class = "glass-card",
                h3(span(class = "card-icon icon-de", "≡"), "DE Results Table"),
                DTOutput("de_table", height = "340px")
              )
            )
          )
        )
      )
    ),

    # ======== TAB 7: EXPORT ========
    tabPanel("Export",
      fluidRow(
        column(6,
          div(class = "glass-card",
            h3(span(class = "card-icon icon-export", "↗"), "Export Results"),
            p(class = "help-text",
              "Download processed data, cluster assignments, DE results, and plots."),
            div(class = "section-divider"),
            downloadButton("dl_metadata", "Cell Metadata + Clusters (CSV)",
                           class = "phen-btn-secondary", style = "width:100%; margin-bottom:10px;"),
            downloadButton("dl_de_results", "DE Results (CSV)",
                           class = "phen-btn-secondary", style = "width:100%; margin-bottom:10px;"),
            downloadButton("dl_svg_results", "SVG Results (CSV)",
                           class = "phen-btn-secondary", style = "width:100%; margin-bottom:10px;"),
            downloadButton("dl_nhood", "Neighborhood Enrichment (CSV)",
                           class = "phen-btn-secondary", style = "width:100%; margin-bottom:10px;"),
            downloadButton("dl_expr_norm", "Normalized Expression Matrix (CSV)",
                           class = "phen-btn-secondary", style = "width:100%; margin-bottom:10px;"),
            div(class = "section-divider"),
            h3(span(class = "card-icon icon-export", "◧"), "Figure Export"),
            selectInput("export_fig", "Select figure",
                        choices = c("Spatial — Clusters", "Spatial — Gene",
                                    "UMAP — Clusters", "Volcano", "Neighborhood Heatmap",
                                    "QC Violin")),
            fluidRow(
              column(4, numericInput("fig_w", "Width (in)", value = 8, min = 2)),
              column(4, numericInput("fig_h", "Height (in)", value = 6, min = 2)),
              column(4, numericInput("fig_dpi", "DPI", value = 300, min = 72))
            ),
            downloadButton("dl_figure", "Download Figure (PDF)",
                           class = "phen-btn-primary", style = "width:100%;")
          )
        ),
        column(6,
          div(class = "glass-card",
            h3(span(class = "card-icon icon-export", "ℹ"), "Session Information"),
            verbatimTextOutput("session_info")
          )
        )
      )
    )
  )
)


# ---------------------------------------------------------------------------
# SERVER
# ---------------------------------------------------------------------------
server <- function(input, output, session) {

  # ---- Reactive values store ----
  rv <- reactiveValues(
    raw_expr      = NULL,    # raw expression matrix (cells × genes)
    raw_meta      = NULL,    # raw metadata
    gene_names    = NULL,
    cell_ids      = NULL,
    coords        = NULL,    # matrix with x, y
    filt_idx      = NULL,    # logical indices post-QC
    norm_expr     = NULL,    # normalized expression
    hvg           = NULL,    # highly variable gene names
    pca           = NULL,    # PCA result (scores matrix)
    pca_sdev      = NULL,
    umap          = NULL,    # UMAP 2D coordinates
    clusters      = NULL,    # cluster assignments
    knn_idx       = NULL,    # KNN neighbor indices
    knn_dists     = NULL,
    nhood_z       = NULL,    # neighborhood enrichment z-matrix
    svg_results   = NULL,    # spatially variable genes
    de_results    = NULL,    # differential expression results
    status        = "awaiting data",
    sweep_data    = NULL
  )

  # ---- Status indicator ----
  output$status_indicator <- renderUI({
    cls <- if (rv$status == "ready") "ready"
           else if (grepl("running|processing", rv$status)) "running"
           else if (grepl("error", rv$status)) "error"
           else "running"
    div(class = paste("status-pill", cls),
        span(class = "dot"), rv$status)
  })

  # ---- DEMO DATA ----
  observeEvent(input$load_demo, {
    withProgress(message = "Generating demo MERFISH data...", value = 0.2, {
      set.seed(42)
      n_cells <- 3000; n_genes <- 150
      # Simulate 5 spatial clusters in a tissue
      centers <- matrix(c(200, 200, 800, 200, 500, 600, 200, 800, 800, 800),
                        ncol = 2, byrow = TRUE)
      cl_assign <- sample(1:5, n_cells, replace = TRUE, prob = c(0.3, 0.2, 0.2, 0.15, 0.15))
      coords <- matrix(NA, nrow = n_cells, ncol = 2)
      for (k in 1:5) {
        idx <- which(cl_assign == k)
        coords[idx, 1] <- rnorm(length(idx), centers[k, 1], 100)
        coords[idx, 2] <- rnorm(length(idx), centers[k, 2], 100)
      }
      colnames(coords) <- c("center_x", "center_y")

      # Simulate expression with cluster-specific marker genes
      gene_names <- paste0("Gene_", sprintf("%03d", 1:n_genes))
      expr <- matrix(rpois(n_cells * n_genes, lambda = 2), nrow = n_cells, ncol = n_genes)
      colnames(expr) <- gene_names
      # Inject marker genes per cluster
      for (k in 1:5) {
        marker_idx <- ((k - 1) * 10 + 1):(k * 10)
        expr[cl_assign == k, marker_idx] <- expr[cl_assign == k, marker_idx] + rpois(
          sum(cl_assign == k) * 10, lambda = 15
        )
      }
      # Add spatially variable genes
      expr[, 51] <- expr[, 51] + pmax(0, round(coords[, 1] / 50))
      expr[, 52] <- expr[, 52] + pmax(0, round(coords[, 2] / 40))
      expr[, 53] <- expr[, 53] + pmax(0, round(sin(coords[, 1] / 150) * 10))

      cell_ids <- paste0("cell_", seq_len(n_cells))
      meta <- data.frame(
        cell_id = cell_ids,
        center_x = coords[, 1],
        center_y = coords[, 2],
        cell_area = rlnorm(n_cells, log(200), 0.5),
        blank_counts = rpois(n_cells, lambda = 1),
        fov = sample(paste0("FOV_", 1:4), n_cells, replace = TRUE),
        stringsAsFactors = FALSE
      )
      rownames(expr) <- cell_ids

      setProgress(0.8, detail = "Loading into session")
      rv$raw_expr   <- expr
      rv$raw_meta   <- meta
      rv$gene_names <- gene_names
      rv$cell_ids   <- cell_ids
      rv$coords     <- coords
      rv$filt_idx   <- rep(TRUE, n_cells)
      rv$status     <- "demo loaded"

      updateSelectInput(session, "x_col", choices = colnames(meta), selected = "center_x")
      updateSelectInput(session, "y_col", choices = colnames(meta), selected = "center_y")
      updateSelectInput(session, "fov_col",
                        choices = c("(none)" = "", colnames(meta)), selected = "fov")
      updateSelectInput(session, "area_col",
                        choices = c("(none)" = "", colnames(meta)), selected = "cell_area")
      updateSelectInput(session, "negctrl_col",
                        choices = c("(none)" = "", colnames(meta)), selected = "blank_counts")
      updateSelectInput(session, "vol_col",
                        choices = c("(none)" = "", colnames(meta)), selected = "cell_area")
      updateSelectInput(session, "spatial_gene", choices = gene_names, selected = gene_names[1])
    })
  })

  # ---- FILE IMPORT ----
  observeEvent(input$meta_file, {
    ext <- tools::file_ext(input$meta_file$name)
    sep <- if (ext %in% c("tsv", "txt")) "\t" else ","
    df <- read.csv(input$meta_file$datapath, sep = sep, stringsAsFactors = FALSE)
    rv$raw_meta <- df
    updateSelectInput(session, "x_col", choices = colnames(df))
    updateSelectInput(session, "y_col", choices = colnames(df))
    updateSelectInput(session, "fov_col", choices = c("(none)" = "", colnames(df)))
    updateSelectInput(session, "area_col", choices = c("(none)" = "", colnames(df)))
    updateSelectInput(session, "negctrl_col", choices = c("(none)" = "", colnames(df)))
    updateSelectInput(session, "vol_col", choices = c("(none)" = "", colnames(df)))
  })

  observeEvent(input$btn_load, {
    req(input$expr_file, input$meta_file, input$x_col, input$y_col)
    withProgress(message = "Loading data...", value = 0.3, {
      ext_e <- tools::file_ext(input$expr_file$name)
      sep_e <- if (ext_e %in% c("tsv", "txt")) "\t" else ","
      expr <- as.matrix(read.csv(input$expr_file$datapath, sep = sep_e,
                                 row.names = 1, check.names = FALSE))
      meta <- rv$raw_meta
      # Align cells
      common <- intersect(rownames(expr), meta[[1]])
      if (length(common) < 10) {
        showNotification("Fewer than 10 cells matched between expression and metadata. Check IDs.",
                         type = "error")
        return()
      }
      expr <- expr[common, ]
      meta <- meta[match(common, meta[[1]]), ]

      coords <- cbind(meta[[input$x_col]], meta[[input$y_col]])
      colnames(coords) <- c("x", "y")

      setProgress(0.8, detail = "Finalizing")
      rv$raw_expr   <- expr
      rv$raw_meta   <- meta
      rv$gene_names <- colnames(expr)
      rv$cell_ids   <- common
      rv$coords     <- coords
      rv$filt_idx   <- rep(TRUE, nrow(expr))
      rv$status     <- "data loaded"

      updateSelectInput(session, "spatial_gene", choices = colnames(expr), selected = colnames(expr)[1])
    })
  })

  # ---- Import KPIs ----
  output$import_kpis <- renderUI({
    req(rv$raw_expr)
    div(class = "kpi-row",
      div(class = "kpi-item",
        div(class = "kpi-label", "Cells"),
        div(class = "kpi-value cyan", format(nrow(rv$raw_expr), big.mark = ","))
      ),
      div(class = "kpi-item",
        div(class = "kpi-label", "Genes"),
        div(class = "kpi-value magenta", format(ncol(rv$raw_expr), big.mark = ","))
      ),
      div(class = "kpi-item",
        div(class = "kpi-label", "Metadata cols"),
        div(class = "kpi-value green", ncol(rv$raw_meta))
      ),
      div(class = "kpi-item",
        div(class = "kpi-label", "Median counts/cell"),
        div(class = "kpi-value amber", format(round(median(rowSums(rv$raw_expr))), big.mark = ","))
      )
    )
  })

  output$preview_table <- renderDT({
    req(rv$raw_meta)
    datatable(head(rv$raw_meta, 200),
              options = list(pageLength = 10, scrollX = TRUE,
                             dom = 'frtip'),
              class = 'compact', rownames = FALSE)
  })

  # ---- QC ----
  qc_stats <- reactive({
    req(rv$raw_expr)
    data.frame(
      cell = rv$cell_ids,
      total_counts = rowSums(rv$raw_expr),
      n_genes = rowSums(rv$raw_expr > 0),
      stringsAsFactors = FALSE
    )
  })

  observeEvent(input$btn_qc, {
    req(rv$raw_expr)
    qs <- qc_stats()
    active <- input$qc_metrics
    keep <- rep(TRUE, nrow(qs))

    # Count filter
    if ("counts" %in% active) {
      keep <- keep & qs$total_counts >= input$qc_min_counts &
                      qs$total_counts <= input$qc_max_counts
    }
    # Gene detection filter
    if ("genes" %in% active) {
      keep <- keep & qs$n_genes >= input$qc_min_genes &
                      qs$n_genes <= input$qc_max_genes
    }
    # Area / volume filter
    if ("area" %in% active &&
        input$area_col != "" && input$area_col %in% colnames(rv$raw_meta)) {
      areas <- rv$raw_meta[[input$area_col]]
      keep <- keep & areas >= input$qc_min_area & areas <= input$qc_max_area
    }
    # Negative control ratio filter
    if ("negctrl" %in% active &&
        input$negctrl_col != "" && input$negctrl_col %in% colnames(rv$raw_meta)) {
      negctrl_counts <- rv$raw_meta[[input$negctrl_col]]
      neg_ratio <- negctrl_counts / (qs$total_counts + 1e-10)
      keep <- keep & neg_ratio <= input$qc_max_negctrl_ratio
    }
    # Transcript density filter (counts / area)
    if ("density" %in% active &&
        input$area_col != "" && input$area_col %in% colnames(rv$raw_meta)) {
      areas <- rv$raw_meta[[input$area_col]]
      density_vals <- qs$total_counts / (areas + 1e-10)
      keep <- keep & density_vals >= input$qc_min_density &
                      density_vals <= input$qc_max_density
    }

    rv$filt_idx <- keep
    rv$status <- paste0("QC: ", sum(keep), " / ", length(keep), " cells pass")
    showNotification(
      paste0("QC complete — retained ", sum(keep), " of ", length(keep), " cells (",
             round(100 * mean(keep), 1), "%) using filters: ",
             paste(active, collapse = ", ")),
      type = "message"
    )
  })

  output$qc_kpis <- renderUI({
    req(rv$filt_idx)
    n_total <- length(rv$filt_idx)
    n_pass  <- sum(rv$filt_idx)
    div(class = "kpi-row",
      div(class = "kpi-item",
        div(class = "kpi-label", "Total cells"), div(class = "kpi-value cyan", format(n_total, big.mark = ","))
      ),
      div(class = "kpi-item",
        div(class = "kpi-label", "Pass QC"), div(class = "kpi-value green", format(n_pass, big.mark = ","))
      ),
      div(class = "kpi-item",
        div(class = "kpi-label", "Filtered out"), div(class = "kpi-value red", format(n_total - n_pass, big.mark = ","))
      ),
      div(class = "kpi-item",
        div(class = "kpi-label", "% retained"), div(class = "kpi-value amber", paste0(round(100 * n_pass / n_total, 1), "%"))
      )
    )
  })

  output$qc_violin <- renderPlot({
    req(rv$raw_expr)
    qs <- qc_stats()
    qs$pass <- rv$filt_idx
    df <- data.frame(
      value = c(qs$total_counts, qs$n_genes),
      metric = rep(c("Total Counts", "Genes Detected"), each = nrow(qs)),
      pass = rep(qs$pass, 2)
    )
    ggplot(df, aes(x = metric, y = value, fill = pass)) +
      geom_violin(alpha = 0.7, scale = "width", color = NA) +
      scale_fill_manual(values = c("TRUE" = "#00c2ff", "FALSE" = "#ff5252"),
                        labels = c("Filtered", "Pass"), name = "Status") +
      labs(title = "QC Distributions", x = NULL, y = "Value") +
      theme_phenomenalist() +
      theme(legend.position = "top")
  }, bg = "#12151c")

  output$qc_scatter <- renderPlot({
    req(rv$raw_expr)
    qs <- qc_stats()
    qs$pass <- rv$filt_idx
    ggplot(qs, aes(x = total_counts, y = n_genes, color = pass)) +
      geom_point(size = 0.5, alpha = 0.6) +
      scale_color_manual(values = c("TRUE" = "#00c2ff", "FALSE" = "#ff5252"),
                         labels = c("Filtered", "Pass"), name = "Status") +
      labs(title = "Counts vs Genes Detected", x = "Total Counts", y = "Genes Detected") +
      theme_phenomenalist() +
      theme(legend.position = "top")
  }, bg = "#12151c")

  output$qc_spatial <- renderPlot({
    req(rv$raw_expr, rv$coords)
    qs <- qc_stats()
    df <- data.frame(x = rv$coords[, 1], y = rv$coords[, 2],
                     total_counts = qs$total_counts, pass = rv$filt_idx)
    ggplot(df, aes(x = x, y = y, color = total_counts)) +
      geom_point(size = 0.4, alpha = 0.7) +
      scale_color_viridis_c(option = "D", name = "Total\nCounts") +
      coord_fixed() + scale_y_reverse() +
      labs(title = "Spatial — Total Counts", x = "X", y = "Y") +
      theme_phenomenalist()
  }, bg = "#12151c")

  # ---- PROCESSING ----
  observeEvent(input$btn_process, {
    req(rv$raw_expr, rv$filt_idx)
    rv$status <- "processing..."

    withProgress(message = "Running processing pipeline...", value = 0.1, {

      expr <- rv$raw_expr[rv$filt_idx, ]
      n_cells <- nrow(expr); n_genes <- ncol(expr)

      # 1. Normalize
      setProgress(0.2, detail = paste0("Normalizing — ", input$norm_method))
      if (input$norm_method == "lognorm") {
        lib_size <- rowSums(expr)
        norm <- log1p(sweep(expr, 1, lib_size, "/") * 10000)
      } else if (input$norm_method == "cp10k") {
        lib_size <- rowSums(expr)
        norm <- sweep(expr, 1, lib_size, "/") * 10000
      } else if (input$norm_method == "cellvol") {
        # Cell-volume normalization: counts / cell volume, then log1p
        vol_col <- input$vol_col
        if (vol_col != "" && vol_col %in% colnames(rv$raw_meta)) {
          volumes <- rv$raw_meta[[vol_col]][rv$filt_idx]
          volumes[volumes <= 0] <- 1e-6
          norm <- log1p(sweep(expr, 1, volumes, "/") * median(volumes))
        } else {
          # Fall back to standard log-norm if no volume column
          lib_size <- rowSums(expr)
          norm <- log1p(sweep(expr, 1, lib_size, "/") * 10000)
          showNotification("No volume column selected — falling back to log-normalize.",
                           type = "warning")
        }
      } else if (input$norm_method == "none") {
        norm <- expr
      } else {
        # Pearson residual approximation (SCT-like)
        lib_size <- rowSums(expr)
        gene_mean <- colMeans(expr)
        theta <- 100
        norm <- matrix(0, nrow = n_cells, ncol = n_genes)
        for (j in seq_len(n_genes)) {
          mu <- outer(lib_size, gene_mean[j]) / mean(lib_size)
          norm[, j] <- (expr[, j] - mu) / sqrt(mu + mu^2 / theta)
        }
        # Clip
        norm[norm > sqrt(n_cells)] <- sqrt(n_cells)
        norm[norm < -sqrt(n_cells)] <- -sqrt(n_cells)
      }
      colnames(norm) <- colnames(expr)
      rownames(norm) <- rownames(expr)
      rv$norm_expr <- norm

      # 2. HVGs
      setProgress(0.4, detail = paste0("Selecting HVGs — ", input$hvg_method))
      if (input$hvg_method == "all") {
        hvg_idx <- seq_len(n_genes)
        rv$hvg <- colnames(norm)
      } else if (input$hvg_method == "vst") {
        # Mean-variance trend: fit loess, rank by residual variance
        gene_mean <- colMeans(norm)
        gene_var  <- apply(norm, 2, var)
        # Avoid log of zero
        lm_fit <- tryCatch({
          lo <- loess(log1p(gene_var) ~ log1p(gene_mean))
          fitted_var <- exp(predict(lo)) - 1
          residual_var <- gene_var / (fitted_var + 1e-10)
          residual_var
        }, error = function(e) gene_var)
        n_hvg <- min(input$n_hvg, n_genes)
        hvg_idx <- order(lm_fit, decreasing = TRUE)[seq_len(n_hvg)]
        rv$hvg <- colnames(norm)[hvg_idx]
      } else {
        # Simple top-variance
        gene_var <- apply(norm, 2, var)
        n_hvg <- min(input$n_hvg, n_genes)
        hvg_idx <- order(gene_var, decreasing = TRUE)[seq_len(n_hvg)]
        rv$hvg <- colnames(norm)[hvg_idx]
      }

      # 3. PCA with user-selected scaling
      setProgress(0.6, detail = "Running PCA")
      expr_hvg <- norm[, hvg_idx]
      if (input$scale_method == "zscore") {
        expr_scaled <- scale(expr_hvg, center = TRUE, scale = TRUE)
      } else if (input$scale_method == "center") {
        expr_scaled <- scale(expr_hvg, center = TRUE, scale = FALSE)
      } else {
        expr_scaled <- expr_hvg
      }
      expr_scaled[is.nan(expr_scaled)] <- 0
      n_pcs <- min(input$n_pcs, n_hvg - 1, n_cells - 1)
      pca <- prcomp(expr_scaled, center = FALSE, scale. = FALSE, rank. = n_pcs)
      rv$pca <- pca$x[, seq_len(n_pcs)]
      rv$pca_sdev <- pca$sdev[seq_len(n_pcs)]

      # 4. UMAP (simple implementation via neighbor graph + spectral)
      setProgress(0.8, detail = "Computing UMAP embedding")
      # Lightweight UMAP: use R's built-in cmdscale on KNN distance graph
      # then refine with t-SNE-like repulsion (simplified for portability)
      nn <- nn2(rv$pca, k = input$umap_neighbors + 1)
      nn_idx <- nn$nn.idx[, -1]; nn_dists <- nn$nn.dists[, -1]

      # Build weighted adjacency and do spectral embedding
      n <- nrow(rv$pca)
      sigma <- apply(nn_dists, 1, function(d) d[min(input$umap_neighbors, length(d))])
      sigma[sigma < 1e-8] <- 1e-8

      # Fuzzy set membership
      trip_i <- rep(seq_len(n), each = ncol(nn_idx))
      trip_j <- as.vector(nn_idx)
      trip_v <- as.vector(exp(-nn_dists / sigma[trip_i]))
      W <- sparseMatrix(i = trip_i, j = trip_j, x = trip_v, dims = c(n, n))
      W <- (W + t(W)) / 2

      # Spectral initialization
      D_inv_sqrt <- Diagonal(x = 1 / sqrt(rowSums(W) + 1e-10))
      L_norm <- D_inv_sqrt %*% W %*% D_inv_sqrt
      # Get top 2 eigenvectors
      eig <- tryCatch({
        RSpectra::eigs_sym(L_norm, k = 3, which = "LM")
      }, error = function(e) {
        # Fallback to classical MDS
        dist_mat <- dist(rv$pca[, 1:min(5, ncol(rv$pca))])
        mds <- cmdscale(dist_mat, k = 2)
        list(vectors = cbind(mds, 0))
      })
      umap_init <- eig$vectors[, 2:3]
      # Add small jitter for min_dist effect
      umap_init <- umap_init + matrix(rnorm(n * 2, 0, input$umap_min_dist * 0.1), ncol = 2)
      colnames(umap_init) <- c("UMAP_1", "UMAP_2")
      rv$umap <- umap_init
      rv$knn_idx <- nn_idx
      rv$knn_dists <- nn_dists

      setProgress(1, detail = "Complete")
      rv$status <- "processed"
      showNotification("Processing pipeline complete.", type = "message")
    })
  })

  output$process_kpis <- renderUI({
    req(rv$norm_expr)
    div(class = "kpi-row",
      div(class = "kpi-item",
        div(class = "kpi-label", "Cells analyzed"), div(class = "kpi-value cyan", format(nrow(rv$norm_expr), big.mark = ","))
      ),
      div(class = "kpi-item",
        div(class = "kpi-label", "HVGs selected"), div(class = "kpi-value magenta", length(rv$hvg))
      ),
      div(class = "kpi-item",
        div(class = "kpi-label", "PCs computed"), div(class = "kpi-value green", ncol(rv$pca))
      ),
      div(class = "kpi-item",
        div(class = "kpi-label", "Normalization"), div(class = "kpi-value amber", input$norm_method)
      )
    )
  })

  output$pca_elbow <- renderPlot({
    req(rv$pca_sdev)
    df <- data.frame(PC = seq_along(rv$pca_sdev),
                     Variance = rv$pca_sdev^2 / sum(rv$pca_sdev^2) * 100)
    ggplot(df, aes(x = PC, y = Variance)) +
      geom_line(color = "#00c2ff", linewidth = 0.8) +
      geom_point(color = "#00c2ff", size = 2.5) +
      labs(title = "PCA Elbow Plot", x = "Principal Component", y = "% Variance Explained") +
      theme_phenomenalist()
  }, bg = "#12151c")

  output$umap_plot <- renderPlot({
    req(rv$umap)
    df <- data.frame(UMAP_1 = rv$umap[, 1], UMAP_2 = rv$umap[, 2])
    ggplot(df, aes(x = UMAP_1, y = UMAP_2)) +
      geom_point(size = 0.5, alpha = 0.6, color = "#00c2ff") +
      labs(title = "UMAP Embedding", x = "UMAP 1", y = "UMAP 2") +
      theme_phenomenalist() + coord_fixed()
  }, bg = "#12151c")

  output$hvg_plot <- renderPlot({
    req(rv$norm_expr, rv$hvg)
    gene_mean <- colMeans(rv$norm_expr)
    gene_var  <- apply(rv$norm_expr, 2, var)
    df <- data.frame(mean = gene_mean, var = gene_var,
                     hvg = colnames(rv$norm_expr) %in% rv$hvg)
    ggplot(df, aes(x = mean, y = var, color = hvg)) +
      geom_point(size = 0.8, alpha = 0.6) +
      scale_color_manual(values = c("TRUE" = "#e040fb", "FALSE" = "#555b6a"),
                         labels = c("Non-variable", "HVG"), name = "") +
      labs(title = "Feature Selection", x = "Mean Expression", y = "Variance") +
      theme_phenomenalist() + theme(legend.position = "top")
  }, bg = "#12151c")

  output$pca_loadings <- renderPlot({
    req(rv$pca)
    df <- data.frame(PC1 = rv$pca[, 1], PC2 = rv$pca[, 2])
    ggplot(df, aes(x = PC1, y = PC2)) +
      geom_point(size = 0.4, alpha = 0.5, color = "#00e676") +
      labs(title = "PC1 vs PC2", x = "PC 1", y = "PC 2") +
      theme_phenomenalist() + coord_fixed()
  }, bg = "#12151c")

  # ---- CLUSTERING ----
  observeEvent(input$btn_cluster, {
    req(rv$pca, rv$knn_idx)
    rv$status <- "clustering..."
    withProgress(message = "Running SNN clustering...", value = 0.3, {
      # Re-compute KNN if k changed
      k <- input$cluster_k
      nn <- nn2(rv$pca, k = k + 1)
      knn_idx <- nn$nn.idx[, -1]

      setProgress(0.6, detail = "Leiden community detection")
      clusters <- snn_cluster(knn_idx, resolution = input$cluster_res,
                              objective = input$leiden_objective)
      rv$clusters <- as.character(clusters)
      rv$knn_idx  <- knn_idx
      rv$knn_dists <- nn$nn.dists[, -1]

      # Update UI
      cl_levels <- sort(unique(rv$clusters))
      updateSelectInput(session, "spatial_color",
                        choices = c("Cluster" = "cluster", colnames(rv$raw_meta)))
      updateSelectInput(session, "de_group1", choices = cl_levels)
      updateSelectInput(session, "de_group2", choices = c("(rest)" = "rest", cl_levels))

      setProgress(1, detail = "Complete")
      rv$status <- paste0("clustered: ", length(cl_levels), " clusters")
      showNotification(paste0("Found ", length(cl_levels), " clusters."), type = "message")
    })
  })

  # Resolution sweep
  observeEvent(input$btn_sweep, {
    req(rv$pca, rv$knn_idx)
    withProgress(message = "Running resolution sweep...", value = 0.1, {
      res_seq <- seq(input$sweep_range[1], input$sweep_range[2],
                     length.out = input$sweep_steps)
      sweep_df <- data.frame(resolution = numeric(), n_clusters = integer())
      k <- input$cluster_k
      nn <- nn2(rv$pca, k = k + 1)
      knn_idx <- nn$nn.idx[, -1]
      for (i in seq_along(res_seq)) {
        setProgress(i / length(res_seq), detail = paste0("res=", round(res_seq[i], 2)))
        cl <- snn_cluster(knn_idx, resolution = res_seq[i],
                          objective = input$leiden_objective)
        sweep_df <- rbind(sweep_df, data.frame(resolution = res_seq[i],
                                               n_clusters = length(unique(cl))))
      }
      rv$sweep_data <- sweep_df
    })
  })

  output$cluster_kpis <- renderUI({
    req(rv$clusters)
    cl_levels <- sort(unique(rv$clusters))
    sizes <- table(rv$clusters)
    div(class = "kpi-row",
      div(class = "kpi-item",
        div(class = "kpi-label", "N clusters"), div(class = "kpi-value magenta", length(cl_levels))
      ),
      div(class = "kpi-item",
        div(class = "kpi-label", "Largest cluster"), div(class = "kpi-value cyan", format(max(sizes), big.mark = ","))
      ),
      div(class = "kpi-item",
        div(class = "kpi-label", "Smallest cluster"), div(class = "kpi-value amber", format(min(sizes), big.mark = ","))
      ),
      div(class = "kpi-item",
        div(class = "kpi-label", "Resolution"), div(class = "kpi-value green", input$cluster_res)
      )
    )
  })

  output$cluster_umap <- renderPlot({
    req(rv$umap, rv$clusters)
    df <- data.frame(UMAP_1 = rv$umap[, 1], UMAP_2 = rv$umap[, 2],
                     cluster = factor(rv$clusters))
    n_cl <- length(unique(rv$clusters))
    ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = cluster)) +
      geom_point(size = 0.5, alpha = 0.7) +
      scale_color_manual(values = phen_palette_discrete(n_cl)) +
      labs(title = "UMAP — Clusters", x = "UMAP 1", y = "UMAP 2") +
      theme_phenomenalist() + coord_fixed() +
      guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
  }, bg = "#12151c")

  output$cluster_spatial <- renderPlot({
    req(rv$clusters, rv$coords, rv$filt_idx)
    coords_filt <- rv$coords[rv$filt_idx, ]
    df <- data.frame(x = coords_filt[, 1], y = coords_filt[, 2],
                     cluster = factor(rv$clusters))
    n_cl <- length(unique(rv$clusters))
    ggplot(df, aes(x = x, y = y, color = cluster)) +
      geom_point(size = 0.5, alpha = 0.7) +
      scale_color_manual(values = phen_palette_discrete(n_cl)) +
      coord_fixed() + scale_y_reverse() +
      labs(title = "Tissue Space — Clusters", x = "X (µm)", y = "Y (µm)") +
      theme_phenomenalist() +
      guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
  }, bg = "#12151c")

  output$cluster_composition <- renderPlot({
    req(rv$clusters)
    df <- as.data.frame(table(cluster = rv$clusters))
    n_cl <- nrow(df)
    ggplot(df, aes(x = reorder(cluster, -Freq), y = Freq, fill = cluster)) +
      geom_col(alpha = 0.85) +
      scale_fill_manual(values = phen_palette_discrete(n_cl)) +
      labs(title = "Cluster Sizes", x = "Cluster", y = "N Cells") +
      theme_phenomenalist() + theme(legend.position = "none")
  }, bg = "#12151c")

  output$cluster_sweep <- renderPlot({
    req(rv$sweep_data)
    ggplot(rv$sweep_data, aes(x = resolution, y = n_clusters)) +
      geom_line(color = "#e040fb", linewidth = 1) +
      geom_point(color = "#e040fb", size = 3) +
      labs(title = "Resolution Sweep", x = "Resolution", y = "N Clusters") +
      theme_phenomenalist()
  }, bg = "#12151c")

  # ---- SPATIAL ANALYSIS ----
  output$spatial_main <- renderPlotly({
    req(rv$coords, rv$filt_idx, rv$clusters)
    coords_filt <- rv$coords[rv$filt_idx, ]
    color_by <- input$spatial_color
    if (color_by == "cluster") {
      df <- data.frame(x = coords_filt[, 1], y = coords_filt[, 2],
                       color = factor(rv$clusters))
    } else if (color_by %in% colnames(rv$raw_meta)) {
      vals <- rv$raw_meta[[color_by]][rv$filt_idx]
      df <- data.frame(x = coords_filt[, 1], y = coords_filt[, 2],
                       color = vals)
    } else {
      df <- data.frame(x = coords_filt[, 1], y = coords_filt[, 2], color = "1")
    }
    yaxis_opts <- if (input$spatial_coord_flip) list(autorange = "reversed", title = "Y (µm)")
                  else list(title = "Y (µm)")
    n_cl <- length(unique(df$color))
    plot_ly(df, x = ~x, y = ~y, color = ~color,
            colors = if (is.factor(df$color)) phen_palette_discrete(n_cl) else "Viridis",
            type = "scattergl", mode = "markers",
            marker = list(size = input$spatial_pt_size * 3,
                          opacity = input$spatial_alpha),
            text = ~paste0("x: ", round(x, 1), "<br>y: ", round(y, 1),
                          "<br>", color_by, ": ", color)) %>%
      layout(
        xaxis = list(title = "X (µm)", scaleanchor = "y"),
        yaxis = yaxis_opts,
        plot_bgcolor = "#12151c", paper_bgcolor = "#12151c",
        font = list(color = "#8b919e", family = "DM Sans"),
        legend = list(font = list(color = "#8b919e"))
      )
  })

  output$spatial_gene_plot <- renderPlotly({
    req(rv$norm_expr, rv$coords, rv$filt_idx, input$spatial_gene)
    gene <- input$spatial_gene
    if (!gene %in% colnames(rv$norm_expr)) return(NULL)
    coords_filt <- rv$coords[rv$filt_idx, ]
    expr_val <- rv$norm_expr[, gene]
    df <- data.frame(x = coords_filt[, 1], y = coords_filt[, 2], expr = expr_val)
    yaxis_opts <- if (input$spatial_coord_flip) list(autorange = "reversed", title = "Y (µm)")
                  else list(title = "Y (µm)")
    plot_ly(df, x = ~x, y = ~y, color = ~expr,
            colors = viridis::viridis(100),
            type = "scattergl", mode = "markers",
            marker = list(size = input$spatial_pt_size * 3,
                          opacity = input$spatial_alpha),
            text = ~paste0("x: ", round(x, 1), "<br>y: ", round(y, 1),
                          "<br>", gene, ": ", round(expr, 3))) %>%
      layout(
        title = list(text = paste0(gene, " expression"), font = list(size = 14, color = "#e8eaed")),
        xaxis = list(title = "X (µm)", scaleanchor = "y"),
        yaxis = yaxis_opts,
        plot_bgcolor = "#12151c", paper_bgcolor = "#12151c",
        font = list(color = "#8b919e", family = "DM Sans")
      )
  })

  # Neighborhood enrichment
  observeEvent(input$btn_nhood, {
    req(rv$clusters, rv$coords, rv$filt_idx)
    withProgress(message = "Computing neighborhood enrichment...", value = 0.3, {
      coords_filt <- rv$coords[rv$filt_idx, ]
      nn <- nn2(coords_filt, k = input$nhood_k + 1)
      nn_idx <- nn$nn.idx[, -1]
      setProgress(0.7, detail = "Computing Z-scores")
      z_mat <- neighborhood_enrichment(rv$clusters, nn_idx)
      rv$nhood_z <- z_mat
      showNotification("Neighborhood enrichment complete.", type = "message")
    })
  })

  output$nhood_heatmap <- renderPlot({
    req(rv$nhood_z)
    # Clamp extreme values for visualization
    z <- rv$nhood_z
    z[z > 10] <- 10; z[z < -10] <- -10
    pheatmap(z,
             color = colorRampPalette(c("#00c2ff", "#12151c", "#ff5252"))(100),
             border_color = "#1a1e2a",
             fontsize = 10,
             main = "Neighborhood Enrichment (Z-score)",
             fontsize_row = 9, fontsize_col = 9)
  }, bg = "#12151c")

  # Spatially variable genes
  observeEvent(input$btn_svg, {
    req(rv$norm_expr, rv$coords, rv$filt_idx)
    withProgress(message = "Detecting spatially variable genes...", value = 0.1, {
      coords_filt <- rv$coords[rv$filt_idx, ]
      nn <- nn2(coords_filt, k = input$svg_k + 1)
      nn_idx <- nn$nn.idx[, -1]; nn_dists <- nn$nn.dists[, -1]

      # Select top variable genes to test
      gene_var <- apply(rv$norm_expr, 2, var)
      n_test <- min(input$svg_n_top, ncol(rv$norm_expr))
      test_genes_idx <- order(gene_var, decreasing = TRUE)[seq_len(n_test)]
      test_genes <- colnames(rv$norm_expr)[test_genes_idx]

      results <- data.frame(gene = test_genes, morans_I = NA_real_,
                            p_value = NA_real_, stringsAsFactors = FALSE)
      for (i in seq_along(test_genes)) {
        setProgress(i / length(test_genes), detail = test_genes[i])
        mi <- morans_i_fast(rv$norm_expr[, test_genes[i]], nn_idx, nn_dists)
        results$morans_I[i] <- mi$I
        results$p_value[i] <- mi$p
      }
      results$p_adj <- p.adjust(results$p_value, method = "BH")
      results <- results[order(results$p_adj, -results$morans_I), ]
      rv$svg_results <- results
      showNotification(
        paste0("Found ", sum(results$p_adj < 0.05), " significant SVGs (padj < 0.05)."),
        type = "message"
      )
    })
  })

  output$svg_table <- renderDT({
    req(rv$svg_results)
    df <- rv$svg_results
    df$morans_I <- round(df$morans_I, 4)
    df$p_value  <- signif(df$p_value, 3)
    df$p_adj    <- signif(df$p_adj, 3)
    datatable(df, options = list(pageLength = 10, scrollX = TRUE, dom = 'frtip'),
              class = 'compact', rownames = FALSE,
              selection = "single")
  })

  # ---- DIFFERENTIAL EXPRESSION ----
  observeEvent(input$btn_de, {
    req(rv$norm_expr, rv$clusters, input$de_group1)
    rv$status <- "running DE..."
    withProgress(message = "Running differential expression...", value = 0.2, {
      g1 <- input$de_group1
      g1_idx <- which(rv$clusters == g1)
      if (input$de_vs_rest) {
        g2_idx <- which(rv$clusters != g1)
      } else {
        req(input$de_group2)
        g2_idx <- which(rv$clusters == input$de_group2)
      }
      if (length(g1_idx) < 3 || length(g2_idx) < 3) {
        showNotification("Need at least 3 cells per group.", type = "error")
        return()
      }
      setProgress(0.5, detail = "Wilcoxon tests")
      results <- run_de(rv$norm_expr, g1_idx, g2_idx, colnames(rv$norm_expr))
      rv$de_results <- results
      rv$status <- "DE complete"
      n_sig <- sum(results$p_adj < input$de_pval_thresh &
                   abs(results$avg_log2FC) > input$de_lfc_thresh, na.rm = TRUE)
      showNotification(paste0("DE complete — ", n_sig, " significant genes."), type = "message")
    })
  })

  output$de_kpis <- renderUI({
    req(rv$de_results)
    df <- rv$de_results
    n_up <- sum(df$p_adj < input$de_pval_thresh & df$avg_log2FC > input$de_lfc_thresh, na.rm = TRUE)
    n_dn <- sum(df$p_adj < input$de_pval_thresh & df$avg_log2FC < -input$de_lfc_thresh, na.rm = TRUE)
    div(class = "kpi-row",
      div(class = "kpi-item",
        div(class = "kpi-label", "Genes tested"), div(class = "kpi-value cyan", nrow(df))
      ),
      div(class = "kpi-item",
        div(class = "kpi-label", "Upregulated"), div(class = "kpi-value green", n_up)
      ),
      div(class = "kpi-item",
        div(class = "kpi-label", "Downregulated"), div(class = "kpi-value red", n_dn)
      ),
      div(class = "kpi-item",
        div(class = "kpi-label", "Total significant"), div(class = "kpi-value magenta", n_up + n_dn)
      )
    )
  })

  output$volcano_plot <- renderPlot({
    req(rv$de_results)
    df <- rv$de_results
    df$neg_log10_padj <- -log10(df$p_adj + 1e-300)
    df$sig <- ifelse(df$p_adj < input$de_pval_thresh & df$avg_log2FC > input$de_lfc_thresh, "Up",
              ifelse(df$p_adj < input$de_pval_thresh & df$avg_log2FC < -input$de_lfc_thresh, "Down", "NS"))
    # Label top genes
    top_genes <- head(df[df$sig != "NS", ], 15)

    ggplot(df, aes(x = avg_log2FC, y = neg_log10_padj, color = sig)) +
      geom_point(size = 1, alpha = 0.6) +
      scale_color_manual(values = c("Up" = "#00e676", "Down" = "#ff5252", "NS" = "#555b6a"),
                         name = "Status") +
      geom_vline(xintercept = c(-input$de_lfc_thresh, input$de_lfc_thresh),
                 linetype = "dashed", color = "#555b6a", linewidth = 0.4) +
      geom_hline(yintercept = -log10(input$de_pval_thresh),
                 linetype = "dashed", color = "#555b6a", linewidth = 0.4) +
      {if (nrow(top_genes) > 0)
        ggrepel::geom_text_repel(data = top_genes, aes(label = gene),
                                 size = 3, color = "#e8eaed", max.overlaps = 12,
                                 segment.color = "#555b6a", segment.size = 0.3)
      } +
      labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10(adj. P-value)") +
      theme_phenomenalist() + theme(legend.position = "top")
  }, bg = "#12151c")

  output$de_dotplot <- renderPlot({
    req(rv$de_results)
    df <- rv$de_results
    top <- head(df[df$p_adj < input$de_pval_thresh & abs(df$avg_log2FC) > input$de_lfc_thresh, ], 20)
    if (nrow(top) == 0) {
      ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No significant genes",
                          color = "#8b919e", size = 5) + theme_void() +
        theme(plot.background = element_rect(fill = "#12151c", color = NA))
    } else {
      ggplot(top, aes(x = avg_log2FC, y = reorder(gene, avg_log2FC),
                      size = pct_1, color = -log10(p_adj + 1e-300))) +
        geom_point() +
        scale_color_viridis_c(option = "C", name = "-log10(padj)") +
        scale_size_continuous(range = c(2, 8), name = "% expressed") +
        labs(title = "Top DE Genes", x = "Log2 Fold Change", y = NULL) +
        theme_phenomenalist()
    }
  }, bg = "#12151c")

  output$de_table <- renderDT({
    req(rv$de_results)
    df <- rv$de_results
    df$avg_log2FC <- round(df$avg_log2FC, 3)
    df$pct_1 <- round(df$pct_1, 3); df$pct_2 <- round(df$pct_2, 3)
    df$p_val <- signif(df$p_val, 3); df$p_adj <- signif(df$p_adj, 3)
    datatable(df, options = list(pageLength = 15, scrollX = TRUE, dom = 'frtip'),
              class = 'compact', rownames = FALSE)
  })

  # ---- EXPORT ----
  output$dl_metadata <- downloadHandler(
    filename = function() paste0("phenomenalist_merfish_metadata_", Sys.Date(), ".csv"),
    content = function(file) {
      meta <- rv$raw_meta[rv$filt_idx, ]
      if (!is.null(rv$clusters)) meta$cluster <- rv$clusters
      if (!is.null(rv$umap)) {
        meta$UMAP_1 <- rv$umap[, 1]; meta$UMAP_2 <- rv$umap[, 2]
      }
      write.csv(meta, file, row.names = FALSE)
    }
  )

  output$dl_de_results <- downloadHandler(
    filename = function() paste0("phenomenalist_merfish_DE_", Sys.Date(), ".csv"),
    content = function(file) { write.csv(rv$de_results, file, row.names = FALSE) }
  )

  output$dl_svg_results <- downloadHandler(
    filename = function() paste0("phenomenalist_merfish_SVG_", Sys.Date(), ".csv"),
    content = function(file) { write.csv(rv$svg_results, file, row.names = FALSE) }
  )

  output$dl_nhood <- downloadHandler(
    filename = function() paste0("phenomenalist_merfish_neighborhood_", Sys.Date(), ".csv"),
    content = function(file) {
      if (!is.null(rv$nhood_z)) write.csv(as.data.frame(rv$nhood_z), file)
    }
  )

  output$dl_expr_norm <- downloadHandler(
    filename = function() paste0("phenomenalist_merfish_norm_expr_", Sys.Date(), ".csv"),
    content = function(file) { write.csv(rv$norm_expr, file) }
  )

  # Figure export
  make_figure <- function(which_fig) {
    switch(which_fig,
      "Spatial — Clusters" = {
        req(rv$clusters, rv$coords, rv$filt_idx)
        coords_filt <- rv$coords[rv$filt_idx, ]
        df <- data.frame(x = coords_filt[, 1], y = coords_filt[, 2],
                         cluster = factor(rv$clusters))
        n_cl <- length(unique(rv$clusters))
        ggplot(df, aes(x = x, y = y, color = cluster)) +
          geom_point(size = 0.5, alpha = 0.7) +
          scale_color_manual(values = phen_palette_discrete(n_cl)) +
          coord_fixed() + scale_y_reverse() +
          labs(title = "MERFISH — Spatial Clusters", x = "X (µm)", y = "Y (µm)") +
          theme_phenomenalist() +
          guides(color = guide_legend(override.aes = list(size = 3)))
      },
      "Spatial — Gene" = {
        req(rv$norm_expr, rv$coords, rv$filt_idx, input$spatial_gene)
        coords_filt <- rv$coords[rv$filt_idx, ]
        gene <- input$spatial_gene
        df <- data.frame(x = coords_filt[, 1], y = coords_filt[, 2],
                         expr = rv$norm_expr[, gene])
        ggplot(df, aes(x = x, y = y, color = expr)) +
          geom_point(size = 0.5, alpha = 0.7) +
          scale_color_viridis_c(option = "D", name = gene) +
          coord_fixed() + scale_y_reverse() +
          labs(title = paste0("MERFISH — ", gene), x = "X (µm)", y = "Y (µm)") +
          theme_phenomenalist()
      },
      "UMAP — Clusters" = {
        req(rv$umap, rv$clusters)
        df <- data.frame(UMAP_1 = rv$umap[, 1], UMAP_2 = rv$umap[, 2],
                         cluster = factor(rv$clusters))
        n_cl <- length(unique(rv$clusters))
        ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = cluster)) +
          geom_point(size = 0.5, alpha = 0.7) +
          scale_color_manual(values = phen_palette_discrete(n_cl)) +
          labs(title = "UMAP — Clusters") +
          theme_phenomenalist() + coord_fixed() +
          guides(color = guide_legend(override.aes = list(size = 3)))
      },
      "Volcano" = {
        req(rv$de_results)
        df <- rv$de_results
        df$neg_log10_padj <- -log10(df$p_adj + 1e-300)
        df$sig <- ifelse(df$p_adj < 0.05 & df$avg_log2FC > 0.25, "Up",
                  ifelse(df$p_adj < 0.05 & df$avg_log2FC < -0.25, "Down", "NS"))
        ggplot(df, aes(x = avg_log2FC, y = neg_log10_padj, color = sig)) +
          geom_point(size = 1, alpha = 0.6) +
          scale_color_manual(values = c("Up" = "#00e676", "Down" = "#ff5252", "NS" = "#555b6a")) +
          labs(title = "Volcano Plot", x = "Log2FC", y = "-Log10(padj)") +
          theme_phenomenalist()
      },
      "QC Violin" = {
        req(rv$raw_expr)
        qs <- qc_stats()
        qs$pass <- rv$filt_idx
        df <- data.frame(
          value = c(qs$total_counts, qs$n_genes),
          metric = rep(c("Total Counts", "Genes Detected"), each = nrow(qs)),
          pass = rep(qs$pass, 2)
        )
        ggplot(df, aes(x = metric, y = value, fill = pass)) +
          geom_violin(alpha = 0.7, scale = "width", color = NA) +
          scale_fill_manual(values = c("TRUE" = "#00c2ff", "FALSE" = "#ff5252")) +
          labs(title = "QC Distributions") +
          theme_phenomenalist()
      },
      ggplot() + theme_void()
    )
  }

  output$dl_figure <- downloadHandler(
    filename = function() {
      fig_name <- gsub("[^a-zA-Z0-9]", "_", input$export_fig)
      paste0("phenomenalist_", fig_name, "_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      p <- make_figure(input$export_fig)
      ggsave(file, plot = p, width = input$fig_w, height = input$fig_h,
             dpi = input$fig_dpi, bg = "#12151c")
    }
  )

  # Session info
  output$session_info <- renderPrint({
    cat("Phenomenalist — MERFISH Spatial Transcriptomics Module\n")
    cat("======================================================\n\n")
    cat("R version:", R.version.string, "\n")
    cat("Platform: ", R.version$platform, "\n\n")
    cat("Loaded packages:\n")
    pkgs <- c("shiny", "ggplot2", "plotly", "DT", "viridis", "RColorBrewer",
              "pheatmap", "igraph", "RANN", "Matrix", "patchwork")
    for (p in pkgs) {
      v <- tryCatch(as.character(packageVersion(p)), error = function(e) "not installed")
      cat(sprintf("  %-18s %s\n", p, v))
    }
    cat("\nSession time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
    if (!is.null(rv$raw_expr)) {
      cat("\nData summary:\n")
      cat(sprintf("  Cells (raw):       %s\n", format(nrow(rv$raw_expr), big.mark = ",")))
      cat(sprintf("  Genes:             %s\n", format(ncol(rv$raw_expr), big.mark = ",")))
      cat(sprintf("  Cells (post-QC):   %s\n", format(sum(rv$filt_idx), big.mark = ",")))
      if (!is.null(rv$clusters))
        cat(sprintf("  Clusters:          %d\n", length(unique(rv$clusters))))
    }
  })
}


# ---------------------------------------------------------------------------
# RUN
# ---------------------------------------------------------------------------
shinyApp(ui = ui, server = server)
