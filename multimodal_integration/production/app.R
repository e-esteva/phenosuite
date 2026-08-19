## ============================================================================
## Phenomenalist — Multi-Modal Spatial Integration Module
## ============================================================================
## Integrates CODEX (protein) and MERFISH (transcript) SpatialExperiment
## objects from the same tissue by spatially matching cells and producing
## a unified multi-modal SPE for downstream analysis.
##
## Pipeline: Import → Align → Match → Integrate → Visualize → Export
##
## Dependencies: SpatialExperiment, SingleCellExperiment, RANN, uwot,
##               ggplot2, plotly, DT, shiny, shinyjs
## ============================================================================

suppressPackageStartupMessages({
  library(shiny)
  library(shinyjs)
  library(ggplot2)
  library(plotly)
  library(DT)
  library(viridis)
  library(RColorBrewer)
  library(patchwork)
  library(RANN)
  library(Matrix)
  library(SpatialExperiment)
  library(SingleCellExperiment)
  library(uwot)
  library(igraph)
  library(scales)
})
source('/srv/shiny-server/phenomenalist/utils/provenance.R')

options(shiny.maxRequestSize = 1000 * 1024^2)

# ── PHENOMENALIST THEME SYSTEM ──────────────────────────────────────────────
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
  font-size: 14px; line-height: 1.5;
  -webkit-font-smoothing: antialiased;
  overflow-x: hidden;
}
body::before {
  content: '';
  position: fixed; inset: 0; z-index: -1;
  background:
    radial-gradient(ellipse 800px 600px at 15% 20%, rgba(0,194,255,0.04), transparent),
    radial-gradient(ellipse 600px 800px at 85% 75%, rgba(224,64,251,0.03), transparent),
    radial-gradient(ellipse 500px 500px at 50% 50%, rgba(0,230,118,0.02), transparent);
  animation: bgShift 30s ease-in-out infinite alternate;
}
@keyframes bgShift { 0% { filter: hue-rotate(0deg); } 100% { filter: hue-rotate(25deg); } }

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

.glass-card {
  background: var(--phen-bg-card);
  backdrop-filter: blur(12px) saturate(1.4);
  border: 1px solid var(--phen-border);
  border-radius: var(--phen-radius-lg);
  padding: 20px; margin-bottom: 16px;
  transition: var(--phen-transition);
}
.glass-card:hover {
  border-color: rgba(255,255,255,0.10);
  box-shadow: 0 4px 24px rgba(0,0,0,0.25);
}
.glass-card h3 {
  margin: 0 0 14px 0; font-size: 14px; font-weight: 600;
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
.icon-align   { background: rgba(224,64,251,0.10); color: var(--phen-accent-magenta); }
.icon-match   { background: rgba(0,230,118,0.10); color: var(--phen-accent-green); }
.icon-integrate { background: rgba(255,171,0,0.10); color: var(--phen-accent-amber); }
.icon-export  { background: rgba(255,171,0,0.10); color: var(--phen-accent-amber); }

.nav-tabs {
  border-bottom: 1px solid var(--phen-border) !important;
  padding: 0 20px; background: var(--phen-bg-secondary);
}
.nav-tabs > li > a {
  font-family: var(--phen-font-sans) !important;
  font-size: 12.5px !important; font-weight: 500 !important;
  color: var(--phen-text-secondary) !important;
  border: none !important; padding: 12px 18px !important;
  border-radius: 0 !important;
  transition: var(--phen-transition);
  border-bottom: 2px solid transparent !important;
  background: transparent !important;
}
.nav-tabs > li > a:hover { color: var(--phen-text-primary) !important; background: rgba(255,255,255,0.03) !important; }
.nav-tabs > li.active > a,
.nav-tabs > li.active > a:focus,
.nav-tabs > li.active > a:hover {
  color: var(--phen-accent-cyan) !important;
  background: transparent !important; border: none !important;
  border-bottom: 2px solid var(--phen-accent-cyan) !important;
}
.tab-content { background: transparent !important; padding: 20px; }

.form-control, .selectize-input, .shiny-input-container select {
  background: var(--phen-bg-input) !important;
  color: var(--phen-text-primary) !important;
  border: 1px solid var(--phen-border) !important;
  border-radius: var(--phen-radius-sm) !important;
  font-family: var(--phen-font-sans) !important; font-size: 13px !important;
  transition: var(--phen-transition);
}
.form-control:focus, .selectize-input.focus {
  border-color: var(--phen-border-focus) !important;
  box-shadow: var(--phen-glow-cyan) !important; outline: none !important;
}
.selectize-dropdown { background: var(--phen-bg-secondary) !important; border: 1px solid var(--phen-border) !important; border-radius: var(--phen-radius-sm) !important; }
.selectize-dropdown .option { color: var(--phen-text-primary) !important; }
.selectize-dropdown .option.active { background: rgba(0,194,255,0.12) !important; }
.control-label, .shiny-input-container label {
  color: var(--phen-text-secondary) !important;
  font-size: 12px !important; font-weight: 500 !important; margin-bottom: 4px !important;
}

.btn-primary, .btn-default, .action-button {
  font-family: var(--phen-font-sans) !important;
  font-size: 12.5px !important; font-weight: 600 !important;
  border-radius: var(--phen-radius-sm) !important;
  padding: 8px 20px !important;
  transition: var(--phen-transition); border: none !important; cursor: pointer;
}
.phen-btn-primary {
  background: linear-gradient(135deg, var(--phen-accent-cyan), #0091d5) !important;
  color: #fff !important; box-shadow: 0 2px 12px rgba(0,194,255,0.2);
}
.phen-btn-primary:hover { box-shadow: 0 4px 20px rgba(0,194,255,0.35); transform: translateY(-1px); }
.phen-btn-secondary {
  background: var(--phen-bg-input) !important;
  color: var(--phen-text-primary) !important;
  border: 1px solid var(--phen-border) !important;
}
.phen-btn-secondary:hover { background: rgba(255,255,255,0.06) !important; border-color: rgba(255,255,255,0.12) !important; }

.kpi-row { display: flex; gap: 12px; margin-bottom: 16px; flex-wrap: wrap; }
.kpi-item {
  flex: 1; min-width: 120px;
  background: var(--phen-bg-card); border: 1px solid var(--phen-border);
  border-radius: var(--phen-radius-md); padding: 14px 16px;
  backdrop-filter: blur(8px);
}
.kpi-label { font-size: 11px; color: var(--phen-text-muted); font-weight: 500; text-transform: uppercase; letter-spacing: 0.06em; }
.kpi-value { font-family: var(--phen-font-mono); font-size: 22px; font-weight: 600; margin-top: 4px; }
.kpi-value.cyan    { color: var(--phen-accent-cyan); }
.kpi-value.magenta { color: var(--phen-accent-magenta); }
.kpi-value.green   { color: var(--phen-accent-green); }
.kpi-value.amber   { color: var(--phen-accent-amber); }

table.dataTable { color: var(--phen-text-primary) !important; }
.dataTables_wrapper { color: var(--phen-text-secondary) !important; font-family: var(--phen-font-sans) !important; font-size: 12px; }
table.dataTable thead th {
  background: var(--phen-bg-secondary) !important; color: var(--phen-text-secondary) !important;
  border-bottom: 1px solid var(--phen-border) !important;
  font-size: 11px !important; font-weight: 600 !important; text-transform: uppercase; letter-spacing: 0.04em;
}
table.dataTable tbody td { border-bottom: 1px solid rgba(255,255,255,0.03) !important; font-family: var(--phen-font-mono); font-size: 12px; }
table.dataTable tbody tr:hover { background: rgba(0,194,255,0.04) !important; }

.shiny-notification {
  background: var(--phen-bg-secondary) !important; color: var(--phen-text-primary) !important;
  border: 1px solid var(--phen-border) !important; border-radius: var(--phen-radius-md) !important;
  font-family: var(--phen-font-sans) !important; box-shadow: 0 8px 32px rgba(0,0,0,0.4);
}
.progress-bar { background: linear-gradient(90deg, var(--phen-accent-cyan), var(--phen-accent-magenta)) !important; }

.plot-container {
  background: var(--phen-bg-secondary); border: 1px solid var(--phen-border);
  border-radius: var(--phen-radius-md); padding: 8px; overflow: hidden;
}
.well { background: var(--phen-bg-card) !important; border: 1px solid var(--phen-border) !important; border-radius: var(--phen-radius-md) !important; box-shadow: none !important; }
.status-pill {
  display: inline-flex; align-items: center; gap: 6px;
  font-family: var(--phen-font-mono); font-size: 11px; font-weight: 500;
  padding: 4px 12px; border-radius: 20px;
}
.status-pill.ready   { color: var(--phen-accent-green); background: rgba(0,230,118,0.08); border: 1px solid rgba(0,230,118,0.15); }
.status-pill.running { color: var(--phen-accent-amber); background: rgba(255,171,0,0.08); border: 1px solid rgba(255,171,0,0.15); }
.status-pill.error   { color: var(--phen-accent-red);   background: rgba(255,82,82,0.08);  border: 1px solid rgba(255,82,82,0.15); }
.status-pill .dot { width: 6px; height: 6px; border-radius: 50%; }
.status-pill.ready .dot   { background: var(--phen-accent-green); box-shadow: 0 0 6px var(--phen-accent-green); }
.status-pill.running .dot { background: var(--phen-accent-amber); box-shadow: 0 0 6px var(--phen-accent-amber); animation: pulse 1.5s ease-in-out infinite; }
@keyframes pulse { 0%,100% { opacity: 1; } 50% { opacity: 0.3; } }

::-webkit-scrollbar { width: 6px; height: 6px; }
::-webkit-scrollbar-track { background: transparent; }
::-webkit-scrollbar-thumb { background: rgba(255,255,255,0.08); border-radius: 3px; }
::-webkit-scrollbar-thumb:hover { background: rgba(255,255,255,0.15); }

.js-plotly-plot .plotly .modebar { background: transparent !important; }
.js-plotly-plot .plotly .modebar-btn path { fill: var(--phen-text-muted) !important; }
.js-plotly-plot .plotly .modebar-btn:hover path { fill: var(--phen-accent-cyan) !important; }

.section-divider { height: 1px; background: var(--phen-border); margin: 16px 0; }
.help-text { font-size: 11.5px; color: var(--phen-text-muted); font-style: italic; margin-top: 4px; }
.checkbox label, .radio label { color: var(--phen-text-secondary) !important; font-size: 13px; }
"

# ── THEME + PALETTE HELPERS ─────────────────────────────────────────────────

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

# ═══════════════════════════════════════════════════════════════════════════
# UI
# ═══════════════════════════════════════════════════════════════════════════

ui <- fluidPage(
  useShinyjs(),
  tags$head(tags$style(HTML(phenomenalist_css))),

  # ── Top bar ──
  div(class = "phen-topbar",
    span(class = "logo", "PHENOMENALIST"),
    span(class = "module-tag", "Multi-Modal Integration · CODEX + MERFISH"),
    span(class = "separator"),
    uiOutput("status_indicator")
  ),

  # ── Main tabs ──
  navbarPage(
    title = NULL, id = "main_tabs",

    # ━━━━━━━━━━━ Tab 1: Import & Align ━━━━━━━━━━━
    tabPanel("Import & Align",
      fluidRow(
        column(4,
          div(class = "glass-card",
            h3(span(class = "card-icon icon-import", "\u2193"), "CODEX SPE (Protein)"),
            p(class = "help-text",
              "Upload the SpatialExperiment (.rds) from your CODEX/PhenoCycler analysis."),
            fileInput("codex_file", "CODEX SPE (.rds)", accept = ".rds",
                      width = "100%")
          ),
          div(class = "glass-card",
            h3(span(class = "card-icon icon-import", "\u2193"), "MERFISH SPE (Transcript)"),
            p(class = "help-text",
              "Upload the SpatialExperiment (.rds) from your MERFISH analysis."),
            fileInput("merfish_file", "MERFISH SPE (.rds)", accept = ".rds",
                      width = "100%")
          ),
          div(class = "glass-card",
            h3(span(class = "card-icon icon-align", "\u21c4"), "Spatial Alignment"),
            p(class = "help-text",
              "If coordinates are in different reference frames, apply translation/scaling ",
              "to align the MERFISH coordinates to the CODEX coordinate space."),
            div(class = "section-divider"),
            numericInput("offset_x", "X offset (MERFISH)", value = 0, step = 1),
            numericInput("offset_y", "Y offset (MERFISH)", value = 0, step = 1),
            numericInput("scale_factor", "Scale factor (MERFISH)", value = 1.0,
                         min = 0.01, step = 0.01),
            checkboxInput("flip_y", "Flip MERFISH Y-axis", value = FALSE),
            div(class = "section-divider"),
            actionButton("apply_alignment", "Apply Alignment",
                         class = "phen-btn-primary", style = "width:100%;")
          )
        ),
        column(8,
          uiOutput("import_kpis"),
          fluidRow(
            column(6,
              div(class = "plot-container",
                plotOutput("spatial_codex", height = "400px")
              )
            ),
            column(6,
              div(class = "plot-container",
                plotOutput("spatial_merfish", height = "400px")
              )
            )
          ),
          div(class = "plot-container", style = "margin-top: 16px;",
            plotOutput("spatial_overlay", height = "450px")
          )
        )
      )
    ),

    # ━━━━━━━━━━━ Tab 2: Cell Matching ━━━━━━━━━━━
    tabPanel("Cell Matching",
      fluidRow(
        column(4,
          div(class = "glass-card",
            h3(span(class = "card-icon icon-match", "\u2194"), "Nearest-Neighbor Matching"),
            p(class = "help-text",
              "Match each CODEX cell to the nearest MERFISH cell by spatial distance. ",
              "Cells beyond the distance threshold remain unmatched."),
            div(class = "section-divider"),
            numericInput("match_threshold", "Max distance threshold",
                         value = 15, min = 1, step = 1),
            selectInput("match_direction", "Match direction",
                        choices = c("CODEX → MERFISH" = "codex_to_merfish",
                                    "MERFISH → CODEX" = "merfish_to_codex",
                                    "Mutual nearest neighbors" = "mutual")),
            p(class = "help-text",
              "Mutual NN: only keep pairs where each cell is the other's nearest neighbor. ",
              "More conservative but reduces false matches."),
            div(class = "section-divider"),
            actionButton("run_matching", "Run Cell Matching",
                         class = "phen-btn-primary", style = "width:100%;")
          )
        ),
        column(8,
          uiOutput("matching_kpis"),
          fluidRow(
            column(6,
              div(class = "plot-container",
                plotOutput("match_distances_hist", height = "350px")
              )
            ),
            column(6,
              div(class = "plot-container",
                plotOutput("match_spatial_map", height = "350px")
              )
            )
          ),
          div(class = "plot-container", style = "margin-top: 16px;",
            plotOutput("match_links", height = "400px")
          )
        )
      )
    ),

    # ━━━━━━━━━━━ Tab 3: Integration & Visualization ━━━━━━━━━━━
    tabPanel("Integration",
      fluidRow(
        column(4,
          div(class = "glass-card",
            h3(span(class = "card-icon icon-integrate", "\u2726"), "Joint Embedding"),
            p(class = "help-text",
              "Compute a joint UMAP from concatenated protein + transcript features. ",
              "Optionally weight each modality."),
            div(class = "section-divider"),
            sliderInput("weight_codex", "CODEX weight", min = 0, max = 1,
                        value = 0.5, step = 0.05),
            numericInput("n_neighbors", "UMAP neighbors", value = 30, min = 5),
            numericInput("min_dist", "UMAP min_dist", value = 0.3, min = 0.0,
                         max = 1.0, step = 0.05),
            div(class = "section-divider"),
            actionButton("run_integration", "Compute Joint UMAP",
                         class = "phen-btn-primary", style = "width:100%;"),
            div(class = "section-divider"),
            h3(span(class = "card-icon icon-integrate", "\u25cb"), "Joint Clustering"),
            numericInput("resolution", "Leiden resolution", value = 0.8,
                         min = 0.1, max = 5.0, step = 0.1),
            numericInput("k_snn", "SNN neighbors (k)", value = 20, min = 5),
            actionButton("run_clustering", "Run Joint Clustering",
                         class = "phen-btn-primary", style = "width:100%;")
          ),
          div(class = "glass-card",
            h3(span(class = "card-icon icon-integrate", "\u25ce"), "Visualization"),
            selectInput("color_by", "Color by",
                        choices = c("joint_cluster", "modality")),
            selectInput("plot_feature", "Feature overlay (gene/protein)",
                        choices = NULL)
          )
        ),
        column(8,
          uiOutput("integration_kpis"),
          fluidRow(
            column(6,
              div(class = "plot-container",
                plotOutput("joint_umap", height = "400px")
              )
            ),
            column(6,
              div(class = "plot-container",
                plotOutput("joint_spatial", height = "400px")
              )
            )
          ),
          fluidRow(style = "margin-top: 16px;",
            column(6,
              div(class = "plot-container",
                plotOutput("feature_umap", height = "350px")
              )
            ),
            column(6,
              div(class = "plot-container",
                plotOutput("feature_spatial", height = "350px")
              )
            )
          )
        )
      )
    ),

    # ━━━━━━━━━━━ Tab 4: Export ━━━━━━━━━━━
    tabPanel("Export",
      fluidRow(
        column(6,
          div(class = "glass-card",
            h3(span(class = "card-icon icon-export", "\u2197"), "Export Integrated SPE"),
            p(class = "help-text",
              "Download the integrated multi-modal SpatialExperiment object. ",
              "Contains both protein (CODEX) and transcript (MERFISH) assays, ",
              "joint UMAP, and cluster assignments. Compatible with all downstream ",
              "Phenomenalist modules."),
            div(class = "section-divider"),
            downloadButton("dl_integrated_spe", "Download Integrated SPE (.rds)",
                           class = "phen-btn-primary", style = "width:100%; margin-bottom:10px;"),
            div(class = "section-divider"),
            downloadButton("dl_metadata_csv", "Download Metadata (CSV)",
                           class = "phen-btn-secondary", style = "width:100%; margin-bottom:10px;"),
            downloadButton("dl_matching_csv", "Download Cell Matching Table (CSV)",
                           class = "phen-btn-secondary", style = "width:100%; margin-bottom:10px;")
          )
        ),
        column(6,
          div(class = "glass-card",
            h3(span(class = "card-icon icon-export", "\u25a7"), "Figure Export"),
            selectInput("export_fig", "Select figure",
                        choices = c("Spatial Overlay", "Joint UMAP",
                                    "Joint Spatial Clusters",
                                    "Matching Links", "Distance Histogram")),
            fluidRow(
              column(4, numericInput("fig_w", "Width (in)", value = 8, min = 2)),
              column(4, numericInput("fig_h", "Height (in)", value = 6, min = 2)),
              column(4, numericInput("fig_dpi", "DPI", value = 300, min = 72))
            ),
            downloadButton("dl_figure", "Download Figure (PDF)",
                           class = "phen-btn-primary", style = "width:100%;")
          )
        )
      )
    )
  )
)

# ═══════════════════════════════════════════════════════════════════════════
# SERVER
# ═══════════════════════════════════════════════════════════════════════════

server <- function(input, output, session) {
  prov_dir <- file.path(tempdir(), paste0("multimodal_", session$token))
  dir.create(prov_dir, showWarnings = FALSE)
  tracker <- ProvenanceTracker$new("multimodal_integration", session, prov_dir)

  rv <- reactiveValues(
    codex_spe       = NULL,
    merfish_spe     = NULL,
    codex_coords    = NULL,   # aligned coordinates
    merfish_coords  = NULL,   # aligned coordinates
    match_idx       = NULL,   # data.frame: codex_idx, merfish_idx, distance
    integrated_spe  = NULL,
    joint_umap      = NULL,
    joint_clusters  = NULL,
    all_features    = NULL,   # combined feature names
    status          = "awaiting data"
  )

  # ── Status indicator ──────────────────────────────────────────────────
  output$status_indicator <- renderUI({
    cls <- if (rv$status == "awaiting data") "running"
           else if (grepl("error", rv$status, ignore.case = TRUE)) "error"
           else "ready"
    div(class = paste("status-pill", cls),
        span(class = "dot"), rv$status)
  })

  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  # TAB 1: Import & Align
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  # ── Load CODEX SPE ──
  observeEvent(input$codex_file, {
    tracker$register_input(input$codex_file, input_id = "codex_file")
    tryCatch({
      obj <- readRDS(input$codex_file$datapath)
      if (!is(obj, "SpatialExperiment") && !is(obj, "SingleCellExperiment")) {
        showNotification("CODEX file must be a SpatialExperiment or SingleCellExperiment.",
                         type = "error")
        return()
      }
      rv$codex_spe <- obj
      rv$codex_coords <- as.matrix(spatialCoords(obj))
      rv$status <- "CODEX loaded"
      showNotification(paste("CODEX SPE loaded:", ncol(obj), "cells,",
                             nrow(obj), "features"), type = "message")
    }, error = function(e) {
      showNotification(paste("Error loading CODEX:", e$message), type = "error")
    })
  })

  # ── Load MERFISH SPE ──
  observeEvent(input$merfish_file, {
    tracker$register_input(input$merfish_file, input_id = "merfish_file")
    tryCatch({
      obj <- readRDS(input$merfish_file$datapath)
      if (!is(obj, "SpatialExperiment") && !is(obj, "SingleCellExperiment")) {
        showNotification("MERFISH file must be a SpatialExperiment or SingleCellExperiment.",
                         type = "error")
        return()
      }
      rv$merfish_spe <- obj
      rv$merfish_coords <- as.matrix(spatialCoords(obj))
      rv$status <- "MERFISH loaded"
      showNotification(paste("MERFISH SPE loaded:", ncol(obj), "cells,",
                             nrow(obj), "features"), type = "message")

      # Populate feature selector with combined features
      codex_feats <- if (!is.null(rv$codex_spe)) rownames(rv$codex_spe) else character(0)
      merfish_feats <- rownames(obj)
      all_feats <- union(codex_feats, merfish_feats)
      rv$all_features <- all_feats
      updateSelectInput(session, "plot_feature", choices = all_feats,
                        selected = all_feats[1])
    }, error = function(e) {
      showNotification(paste("Error loading MERFISH:", e$message), type = "error")
    })
  })

  # ── Import KPIs ──
  output$import_kpis <- renderUI({
    codex_n   <- if (!is.null(rv$codex_spe))   ncol(rv$codex_spe)   else "—"
    codex_f   <- if (!is.null(rv$codex_spe))   nrow(rv$codex_spe)   else "—"
    merfish_n <- if (!is.null(rv$merfish_spe))  ncol(rv$merfish_spe) else "—"
    merfish_f <- if (!is.null(rv$merfish_spe))  nrow(rv$merfish_spe) else "—"
    div(class = "kpi-row",
      div(class = "kpi-item",
        div(class = "kpi-label", "CODEX Cells"),
        div(class = "kpi-value cyan", format(codex_n, big.mark = ","))
      ),
      div(class = "kpi-item",
        div(class = "kpi-label", "CODEX Proteins"),
        div(class = "kpi-value cyan", codex_f)
      ),
      div(class = "kpi-item",
        div(class = "kpi-label", "MERFISH Cells"),
        div(class = "kpi-value magenta", format(merfish_n, big.mark = ","))
      ),
      div(class = "kpi-item",
        div(class = "kpi-label", "MERFISH Genes"),
        div(class = "kpi-value magenta", merfish_f)
      )
    )
  })

  # ── Apply spatial alignment to MERFISH coords ──
  observeEvent(input$apply_alignment, {
    req(rv$merfish_spe)
    coords <- as.matrix(spatialCoords(rv$merfish_spe))
    # Scale
    coords <- coords * input$scale_factor
    # Flip Y
    if (input$flip_y) coords[, 2] <- -coords[, 2]
    # Translate
    coords[, 1] <- coords[, 1] + input$offset_x
    coords[, 2] <- coords[, 2] + input$offset_y
    rv$merfish_coords <- coords
    rv$status <- "aligned"
    showNotification("Alignment applied to MERFISH coordinates.", type = "message")
  })

  # ── Spatial plots ──
  output$spatial_codex <- renderPlot({
    req(rv$codex_coords)
    df <- data.frame(x = rv$codex_coords[, 1], y = rv$codex_coords[, 2])
    ggplot(df, aes(x, y)) +
      geom_point(size = 0.3, alpha = 0.5, color = "#00c2ff") +
      coord_fixed() + scale_y_reverse() +
      labs(title = "CODEX — Spatial Positions", x = "X", y = "Y") +
      theme_phenomenalist()
  })

  output$spatial_merfish <- renderPlot({
    req(rv$merfish_coords)
    df <- data.frame(x = rv$merfish_coords[, 1], y = rv$merfish_coords[, 2])
    ggplot(df, aes(x, y)) +
      geom_point(size = 0.3, alpha = 0.5, color = "#e040fb") +
      coord_fixed() + scale_y_reverse() +
      labs(title = "MERFISH — Spatial Positions", x = "X", y = "Y") +
      theme_phenomenalist()
  })

  output$spatial_overlay <- renderPlot({
    req(rv$codex_coords, rv$merfish_coords)
    df <- rbind(
      data.frame(x = rv$codex_coords[, 1], y = rv$codex_coords[, 2],
                 modality = "CODEX"),
      data.frame(x = rv$merfish_coords[, 1], y = rv$merfish_coords[, 2],
                 modality = "MERFISH")
    )
    ggplot(df, aes(x, y, color = modality)) +
      geom_point(size = 0.3, alpha = 0.4) +
      scale_color_manual(values = c("CODEX" = "#00c2ff", "MERFISH" = "#e040fb")) +
      coord_fixed() + scale_y_reverse() +
      labs(title = "Spatial Overlay — CODEX + MERFISH",
           x = "X", y = "Y") +
      theme_phenomenalist() +
      guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
  })

  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  # TAB 2: Cell Matching
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  observeEvent(input$run_matching, {
    req(rv$codex_coords, rv$merfish_coords)
    withProgress(message = "Matching cells...", value = 0.3, {

      threshold <- input$match_threshold

      if (input$match_direction == "codex_to_merfish") {
        # For each CODEX cell, find nearest MERFISH cell
        nn <- nn2(rv$merfish_coords, rv$codex_coords, k = 1)
        matches <- data.frame(
          codex_idx   = seq_len(nrow(rv$codex_coords)),
          merfish_idx = as.integer(nn$nn.idx[, 1]),
          distance    = nn$nn.dists[, 1]
        )
        matches <- matches[matches$distance <= threshold, ]

      } else if (input$match_direction == "merfish_to_codex") {
        nn <- nn2(rv$codex_coords, rv$merfish_coords, k = 1)
        matches <- data.frame(
          codex_idx   = as.integer(nn$nn.idx[, 1]),
          merfish_idx = seq_len(nrow(rv$merfish_coords)),
          distance    = nn$nn.dists[, 1]
        )
        matches <- matches[matches$distance <= threshold, ]

      } else {
        # Mutual nearest neighbors (vectorized)
        nn_c2m <- nn2(rv$merfish_coords, rv$codex_coords, k = 1)
        nn_m2c <- nn2(rv$codex_coords, rv$merfish_coords, k = 1)
        codex_idx <- seq_len(nrow(rv$codex_coords))
        merfish_target <- as.integer(nn_c2m$nn.idx[, 1])
        # Keep only pairs where MERFISH's NN points back to the same CODEX cell
        mutual_mask <- nn_m2c$nn.idx[merfish_target, 1] == codex_idx &
                       nn_c2m$nn.dists[, 1] <= threshold
        matches <- data.frame(
          codex_idx   = codex_idx[mutual_mask],
          merfish_idx = merfish_target[mutual_mask],
          distance    = nn_c2m$nn.dists[mutual_mask, 1]
        )
      }

      # Remove duplicate MERFISH assignments (keep closest)
      if (nrow(matches) > 0) {
        matches <- matches[order(matches$distance), ]
        matches <- matches[!duplicated(matches$merfish_idx), ]
        matches <- matches[!duplicated(matches$codex_idx), ]
      }

      rv$match_idx <- matches
      setProgress(1, detail = paste(nrow(matches), "pairs matched"))
      rv$status <- paste(nrow(matches), "cells matched")
      showNotification(paste("Matched", nrow(matches), "cell pairs"),
                       type = "message")
    })
  })

  # ── Matching KPIs ──
  output$matching_kpis <- renderUI({
    m <- rv$match_idx
    if (is.null(m)) return(NULL)
    codex_n   <- if (!is.null(rv$codex_spe)) ncol(rv$codex_spe) else 0
    merfish_n <- if (!is.null(rv$merfish_spe)) ncol(rv$merfish_spe) else 0
    div(class = "kpi-row",
      div(class = "kpi-item",
        div(class = "kpi-label", "Matched Pairs"),
        div(class = "kpi-value green", format(nrow(m), big.mark = ","))
      ),
      div(class = "kpi-item",
        div(class = "kpi-label", "CODEX Match Rate"),
        div(class = "kpi-value cyan",
            paste0(round(100 * nrow(m) / max(codex_n, 1), 1), "%"))
      ),
      div(class = "kpi-item",
        div(class = "kpi-label", "MERFISH Match Rate"),
        div(class = "kpi-value magenta",
            paste0(round(100 * nrow(m) / max(merfish_n, 1), 1), "%"))
      ),
      div(class = "kpi-item",
        div(class = "kpi-label", "Median Distance"),
        div(class = "kpi-value amber",
            round(median(m$distance), 2))
      )
    )
  })

  # ── Matching plots ──
  output$match_distances_hist <- renderPlot({
    req(rv$match_idx)
    df <- rv$match_idx
    ggplot(df, aes(x = distance)) +
      geom_histogram(bins = 50, fill = "#00c2ff", color = NA, alpha = 0.8) +
      geom_vline(xintercept = input$match_threshold, linetype = "dashed",
                 color = "#ff5252", linewidth = 0.8) +
      labs(title = "Matching Distance Distribution",
           x = "Distance (coordinate units)", y = "Count") +
      theme_phenomenalist()
  })

  output$match_spatial_map <- renderPlot({
    req(rv$match_idx, rv$codex_coords, rv$merfish_coords)
    m <- rv$match_idx
    codex_matched <- rv$codex_coords[m$codex_idx, ]
    df <- data.frame(x = codex_matched[, 1], y = codex_matched[, 2],
                     distance = m$distance)
    ggplot(df, aes(x, y, color = distance)) +
      geom_point(size = 0.5, alpha = 0.7) +
      scale_color_viridis_c(option = "C", name = "Dist") +
      coord_fixed() + scale_y_reverse() +
      labs(title = "Matched Cells — Distance Map", x = "X", y = "Y") +
      theme_phenomenalist()
  })

  output$match_links <- renderPlot({
    req(rv$match_idx, rv$codex_coords, rv$merfish_coords)
    m <- rv$match_idx
    # Sample links for visualization (max 2000)
    idx <- if (nrow(m) > 2000) sample(nrow(m), 2000) else seq_len(nrow(m))
    m_sub <- m[idx, ]
    df_links <- data.frame(
      x    = rv$codex_coords[m_sub$codex_idx, 1],
      y    = rv$codex_coords[m_sub$codex_idx, 2],
      xend = rv$merfish_coords[m_sub$merfish_idx, 1],
      yend = rv$merfish_coords[m_sub$merfish_idx, 2]
    )
    ggplot() +
      geom_segment(data = df_links, aes(x = x, y = y, xend = xend, yend = yend),
                   color = "#555b6a", alpha = 0.3, linewidth = 0.2) +
      geom_point(data = data.frame(x = rv$codex_coords[m_sub$codex_idx, 1],
                                   y = rv$codex_coords[m_sub$codex_idx, 2]),
                 aes(x, y), color = "#00c2ff", size = 0.4, alpha = 0.6) +
      geom_point(data = data.frame(x = rv$merfish_coords[m_sub$merfish_idx, 1],
                                   y = rv$merfish_coords[m_sub$merfish_idx, 2]),
                 aes(x, y), color = "#e040fb", size = 0.4, alpha = 0.6) +
      coord_fixed() + scale_y_reverse() +
      labs(title = "Cell Matching Links (CODEX=cyan, MERFISH=magenta)",
           subtitle = paste("Showing", nrow(m_sub), "of", nrow(m), "pairs"),
           x = "X", y = "Y") +
      theme_phenomenalist()
  })

  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  # TAB 3: Integration & Visualization
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  # ── Build integrated SPE + joint UMAP ──
  observeEvent(input$run_integration, {
    req(rv$codex_spe, rv$merfish_spe, rv$match_idx)
    tracker$capture_parameters(input)
    if (is.null(tracker$analysis_start)) tracker$analysis_started()
    m <- rv$match_idx
    if (nrow(m) < 10) {
      showNotification("Too few matched cells (<10). Adjust threshold.", type = "error")
      return()
    }

    withProgress(message = "Computing joint embedding...", value = 0.1, {

      # Extract expression for matched cells
      codex_expr <- as.matrix(
        assay(rv$codex_spe, ifelse("exprs" %in% assayNames(rv$codex_spe), "exprs",
                                    assayNames(rv$codex_spe)[1]))
      )[, m$codex_idx]
      merfish_expr <- as.matrix(
        assay(rv$merfish_spe, ifelse("exprs" %in% assayNames(rv$merfish_spe), "exprs",
                                      assayNames(rv$merfish_spe)[1]))
      )[, m$merfish_idx]

      setProgress(0.3, detail = "Scaling features")

      # Scale each modality to unit variance
      codex_scaled <- t(scale(t(codex_expr)))
      merfish_scaled <- t(scale(t(merfish_expr)))
      # Replace NAs from zero-variance features
      codex_scaled[is.na(codex_scaled)] <- 0
      merfish_scaled[is.na(merfish_scaled)] <- 0

      # Weight modalities
      w_codex <- input$weight_codex
      w_merfish <- 1 - w_codex
      codex_weighted <- codex_scaled * w_codex
      merfish_weighted <- merfish_scaled * w_merfish

      # Concatenate features: genes × cells → cells × features for UMAP
      combined <- rbind(codex_weighted, merfish_weighted)
      combined_t <- t(combined)  # cells × features

      setProgress(0.5, detail = "Running UMAP")

      # Joint UMAP
      umap_res <- umap(combined_t,
                        n_neighbors = input$n_neighbors,
                        min_dist    = input$min_dist,
                        n_components = 2,
                        metric      = "cosine",
                        ret_model   = FALSE)
      colnames(umap_res) <- c("UMAP_1", "UMAP_2")
      rv$joint_umap <- umap_res

      setProgress(0.8, detail = "Building integrated SPE")

      # Build integrated SPE using CODEX spatial coordinates
      int_coords <- rv$codex_coords[m$codex_idx, , drop = FALSE]

      # Combine colData
      cd_codex <- as.data.frame(colData(rv$codex_spe))[m$codex_idx, , drop = FALSE]
      cd_merfish <- as.data.frame(colData(rv$merfish_spe))[m$merfish_idx, , drop = FALSE]
      # Prefix columns to avoid collisions
      names(cd_codex) <- paste0("codex_", names(cd_codex))
      names(cd_merfish) <- paste0("merfish_", names(cd_merfish))
      cd_combined <- cbind(cd_codex, cd_merfish)
      cd_combined$modality <- "integrated"

      # Assays: protein features as main, transcript as altExp
      int_spe <- SpatialExperiment(
        assays       = list(exprs = codex_expr, logcounts = codex_expr),
        colData      = DataFrame(cd_combined),
        spatialCoords = int_coords
      )
      # Store MERFISH as alternative experiment
      altExp(int_spe, "merfish") <- SummarizedExperiment(
        assays = list(exprs = merfish_expr, logcounts = merfish_expr)
      )
      reducedDim(int_spe, "UMAP_joint") <- umap_res

      rv$integrated_spe <- int_spe

      # Update feature selector
      all_feats <- c(rownames(codex_expr), rownames(merfish_expr))
      rv$all_features <- all_feats
      updateSelectInput(session, "plot_feature", choices = all_feats,
                        selected = all_feats[1])

      setProgress(1, detail = "Done")
      rv$status <- paste(nrow(m), "cells integrated")
      showNotification(paste("Joint embedding computed for", nrow(m), "cells"),
                       type = "message")
    })
  })

  # ── Joint clustering ──
  observeEvent(input$run_clustering, {
    req(rv$joint_umap, rv$integrated_spe)

    withProgress(message = "Clustering...", value = 0.3, {
      # Build kNN graph from UMAP
      nn <- nn2(rv$joint_umap, k = input$k_snn + 1)
      nn_idx <- nn$nn.idx[, -1, drop = FALSE]  # remove self
      n_cells <- nrow(rv$joint_umap)
      k_neighbors <- ncol(nn_idx)

      setProgress(0.5, detail = "Building graph")

      # Vectorized edge list: rep each source cell k times, flatten neighbors
      edge_mat <- cbind(
        rep(seq_len(n_cells), each = k_neighbors),
        as.integer(t(nn_idx))
      )
      g <- graph_from_edgelist(edge_mat, directed = FALSE)
      g <- simplify(g)

      setProgress(0.7, detail = "Leiden clustering")

      # Community detection (Louvain as fallback since Leiden requires leidenAlg)
      cl <- tryCatch(
        cluster_leiden(g, resolution_parameter = input$resolution),
        error = function(e) cluster_louvain(g, resolution = input$resolution)
      )
      clusters <- as.character(membership(cl))
      rv$joint_clusters <- clusters

      # Store in SPE
      colData(rv$integrated_spe)$joint_cluster <- clusters
      colData(rv$integrated_spe)$cluster <- clusters
      colData(rv$integrated_spe)$Phenotype <- clusters  # PCF compatibility

      # Update color-by choices
      cd_cols <- names(colData(rv$integrated_spe))
      updateSelectInput(session, "color_by",
                        choices = c("joint_cluster", "modality", cd_cols),
                        selected = "joint_cluster")

      setProgress(1, detail = paste(length(unique(clusters)), "clusters"))
      rv$status <- paste(length(unique(clusters)), "joint clusters")
      showNotification(paste("Found", length(unique(clusters)), "joint clusters"),
                       type = "message")
    })
  })

  # ── Integration KPIs ──
  output$integration_kpis <- renderUI({
    n_cells  <- if (!is.null(rv$integrated_spe)) ncol(rv$integrated_spe) else "—"
    n_prot   <- if (!is.null(rv$integrated_spe)) nrow(rv$integrated_spe) else "—"
    n_genes  <- if (!is.null(rv$integrated_spe) && "merfish" %in% altExpNames(rv$integrated_spe))
                  nrow(altExp(rv$integrated_spe, "merfish")) else "—"
    n_clust  <- if (!is.null(rv$joint_clusters)) length(unique(rv$joint_clusters)) else "—"
    div(class = "kpi-row",
      div(class = "kpi-item",
        div(class = "kpi-label", "Integrated Cells"),
        div(class = "kpi-value green", format(n_cells, big.mark = ","))
      ),
      div(class = "kpi-item",
        div(class = "kpi-label", "Proteins (CODEX)"),
        div(class = "kpi-value cyan", n_prot)
      ),
      div(class = "kpi-item",
        div(class = "kpi-label", "Genes (MERFISH)"),
        div(class = "kpi-value magenta", n_genes)
      ),
      div(class = "kpi-item",
        div(class = "kpi-label", "Joint Clusters"),
        div(class = "kpi-value amber", n_clust)
      )
    )
  })

  # ── Visualization plots ──
  output$joint_umap <- renderPlot({
    req(rv$joint_umap, rv$integrated_spe)
    color_col <- input$color_by
    cd <- as.data.frame(colData(rv$integrated_spe))
    df <- data.frame(UMAP_1 = rv$joint_umap[, 1], UMAP_2 = rv$joint_umap[, 2])

    if (color_col %in% names(cd)) {
      df$color <- factor(cd[[color_col]])
      n_col <- length(unique(df$color))
      ggplot(df, aes(UMAP_1, UMAP_2, color = color)) +
        geom_point(size = 0.5, alpha = 0.7) +
        scale_color_manual(values = phen_palette_discrete(n_col), name = color_col) +
        labs(title = paste("Joint UMAP —", color_col)) +
        theme_phenomenalist() + coord_fixed() +
        guides(color = guide_legend(override.aes = list(size = 3)))
    } else {
      ggplot(df, aes(UMAP_1, UMAP_2)) +
        geom_point(size = 0.5, alpha = 0.5, color = "#00c2ff") +
        labs(title = "Joint UMAP") +
        theme_phenomenalist() + coord_fixed()
    }
  })

  output$joint_spatial <- renderPlot({
    req(rv$integrated_spe, rv$joint_clusters)
    coords <- as.matrix(spatialCoords(rv$integrated_spe))
    df <- data.frame(x = coords[, 1], y = coords[, 2],
                     cluster = factor(rv$joint_clusters))
    n_cl <- length(unique(rv$joint_clusters))
    ggplot(df, aes(x, y, color = cluster)) +
      geom_point(size = 0.5, alpha = 0.7) +
      scale_color_manual(values = phen_palette_discrete(n_cl)) +
      coord_fixed() + scale_y_reverse() +
      labs(title = "Joint Spatial Clusters", x = "X", y = "Y") +
      theme_phenomenalist() +
      guides(color = guide_legend(override.aes = list(size = 3)))
  })

  output$feature_umap <- renderPlot({
    req(rv$joint_umap, rv$integrated_spe, input$plot_feature)
    feat <- input$plot_feature
    df <- data.frame(UMAP_1 = rv$joint_umap[, 1], UMAP_2 = rv$joint_umap[, 2])

    # Check if feature is in main assay (CODEX) or altExp (MERFISH)
    if (feat %in% rownames(rv$integrated_spe)) {
      df$expr <- as.numeric(assay(rv$integrated_spe, "exprs")[feat, ])
      source_label <- "(CODEX protein)"
    } else if ("merfish" %in% altExpNames(rv$integrated_spe) &&
               feat %in% rownames(altExp(rv$integrated_spe, "merfish"))) {
      df$expr <- as.numeric(assay(altExp(rv$integrated_spe, "merfish"), "exprs")[feat, ])
      source_label <- "(MERFISH transcript)"
    } else {
      return(NULL)
    }

    ggplot(df, aes(UMAP_1, UMAP_2, color = expr)) +
      geom_point(size = 0.5, alpha = 0.7) +
      scale_color_viridis_c(option = "D", name = feat) +
      labs(title = paste(feat, source_label), subtitle = "Joint UMAP") +
      theme_phenomenalist() + coord_fixed()
  })

  output$feature_spatial <- renderPlot({
    req(rv$integrated_spe, input$plot_feature)
    feat <- input$plot_feature
    coords <- as.matrix(spatialCoords(rv$integrated_spe))
    df <- data.frame(x = coords[, 1], y = coords[, 2])

    if (feat %in% rownames(rv$integrated_spe)) {
      df$expr <- as.numeric(assay(rv$integrated_spe, "exprs")[feat, ])
      source_label <- "(CODEX protein)"
    } else if ("merfish" %in% altExpNames(rv$integrated_spe) &&
               feat %in% rownames(altExp(rv$integrated_spe, "merfish"))) {
      df$expr <- as.numeric(assay(altExp(rv$integrated_spe, "merfish"), "exprs")[feat, ])
      source_label <- "(MERFISH transcript)"
    } else {
      return(NULL)
    }

    ggplot(df, aes(x, y, color = expr)) +
      geom_point(size = 0.5, alpha = 0.7) +
      scale_color_viridis_c(option = "D", name = feat) +
      coord_fixed() + scale_y_reverse() +
      labs(title = paste(feat, source_label), subtitle = "Spatial",
           x = "X", y = "Y") +
      theme_phenomenalist()
  })

  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  # TAB 4: Export
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  # Provenance bundle (provenance.json + replay.R)
  output$dl_provenance <- downloadHandler(
    filename = function() paste0("phenomenalist_multimodal_provenance_", Sys.Date(), ".zip"),
    content = function(file) {
      tracker$capture_parameters(input)
      if (is.null(tracker$analysis_start)) tracker$analysis_started()
      tracker$analysis_completed()
      sidecar_files <- c("provenance.json", "replay.R")
      sidecar_files <- sidecar_files[file.exists(file.path(prov_dir, sidecar_files))]
      if (length(sidecar_files) == 0) {
        writeLines("Provenance sidecar not available.", file)
        return()
      }
      zip::zip(zipfile = file, files = sidecar_files, root = prov_dir)
    },
    contentType = "application/zip"
  )

  output$dl_integrated_spe <- downloadHandler(
    filename = function() paste0("phenomenalist_integrated_spe_", Sys.Date(), ".rds"),
    content = function(file) {
      req(rv$integrated_spe)
      saveRDS(rv$integrated_spe, file)
    }
  )

  output$dl_metadata_csv <- downloadHandler(
    filename = function() paste0("phenomenalist_integrated_metadata_", Sys.Date(), ".csv"),
    content = function(file) {
      req(rv$integrated_spe)
      cd <- as.data.frame(colData(rv$integrated_spe))
      coords <- as.data.frame(spatialCoords(rv$integrated_spe))
      out <- cbind(coords, cd)
      if (!is.null(rv$joint_umap)) {
        out$UMAP_1 <- rv$joint_umap[, 1]
        out$UMAP_2 <- rv$joint_umap[, 2]
      }
      write.csv(out, file, row.names = FALSE)
    }
  )

  output$dl_matching_csv <- downloadHandler(
    filename = function() paste0("phenomenalist_cell_matching_", Sys.Date(), ".csv"),
    content = function(file) {
      req(rv$match_idx)
      write.csv(rv$match_idx, file, row.names = FALSE)
    }
  )

  # ── Figure export ──
  make_figure <- function(which_fig) {
    switch(which_fig,
      "Spatial Overlay" = {
        req(rv$codex_coords, rv$merfish_coords)
        df <- rbind(
          data.frame(x = rv$codex_coords[, 1], y = rv$codex_coords[, 2], modality = "CODEX"),
          data.frame(x = rv$merfish_coords[, 1], y = rv$merfish_coords[, 2], modality = "MERFISH")
        )
        ggplot(df, aes(x, y, color = modality)) +
          geom_point(size = 0.3, alpha = 0.4) +
          scale_color_manual(values = c("CODEX" = "#00c2ff", "MERFISH" = "#e040fb")) +
          coord_fixed() + scale_y_reverse() +
          labs(title = "Spatial Overlay — CODEX + MERFISH", x = "X", y = "Y") +
          theme_phenomenalist() +
          guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
      },
      "Joint UMAP" = {
        req(rv$joint_umap, rv$joint_clusters)
        df <- data.frame(UMAP_1 = rv$joint_umap[, 1], UMAP_2 = rv$joint_umap[, 2],
                         cluster = factor(rv$joint_clusters))
        n_cl <- length(unique(rv$joint_clusters))
        ggplot(df, aes(UMAP_1, UMAP_2, color = cluster)) +
          geom_point(size = 0.5, alpha = 0.7) +
          scale_color_manual(values = phen_palette_discrete(n_cl)) +
          labs(title = "Joint UMAP — Multi-Modal Clusters") +
          theme_phenomenalist() + coord_fixed() +
          guides(color = guide_legend(override.aes = list(size = 3)))
      },
      "Joint Spatial Clusters" = {
        req(rv$integrated_spe, rv$joint_clusters)
        coords <- as.matrix(spatialCoords(rv$integrated_spe))
        df <- data.frame(x = coords[, 1], y = coords[, 2],
                         cluster = factor(rv$joint_clusters))
        n_cl <- length(unique(rv$joint_clusters))
        ggplot(df, aes(x, y, color = cluster)) +
          geom_point(size = 0.5, alpha = 0.7) +
          scale_color_manual(values = phen_palette_discrete(n_cl)) +
          coord_fixed() + scale_y_reverse() +
          labs(title = "Joint Spatial Clusters", x = "X", y = "Y") +
          theme_phenomenalist() +
          guides(color = guide_legend(override.aes = list(size = 3)))
      },
      "Matching Links" = {
        req(rv$match_idx, rv$codex_coords, rv$merfish_coords)
        m <- rv$match_idx
        idx <- if (nrow(m) > 2000) sample(nrow(m), 2000) else seq_len(nrow(m))
        m_sub <- m[idx, ]
        df_links <- data.frame(
          x = rv$codex_coords[m_sub$codex_idx, 1], y = rv$codex_coords[m_sub$codex_idx, 2],
          xend = rv$merfish_coords[m_sub$merfish_idx, 1], yend = rv$merfish_coords[m_sub$merfish_idx, 2]
        )
        ggplot() +
          geom_segment(data = df_links, aes(x = x, y = y, xend = xend, yend = yend),
                       color = "#555b6a", alpha = 0.3, linewidth = 0.2) +
          geom_point(data = data.frame(x = rv$codex_coords[m_sub$codex_idx, 1],
                                       y = rv$codex_coords[m_sub$codex_idx, 2]),
                     aes(x, y), color = "#00c2ff", size = 0.4) +
          geom_point(data = data.frame(x = rv$merfish_coords[m_sub$merfish_idx, 1],
                                       y = rv$merfish_coords[m_sub$merfish_idx, 2]),
                     aes(x, y), color = "#e040fb", size = 0.4) +
          coord_fixed() + scale_y_reverse() +
          labs(title = "Cell Matching Links", x = "X", y = "Y") +
          theme_phenomenalist()
      },
      "Distance Histogram" = {
        req(rv$match_idx)
        ggplot(rv$match_idx, aes(x = distance)) +
          geom_histogram(bins = 50, fill = "#00c2ff", color = NA, alpha = 0.8) +
          labs(title = "Matching Distance Distribution",
               x = "Distance", y = "Count") +
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
}

# ═══════════════════════════════════════════════════════════════════════════
shinyApp(ui = ui, server = server)
