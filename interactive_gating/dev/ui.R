library(shiny)
library(shinydashboard)
library(DT)

options(shiny.maxRequestSize=4000*1024^2)

dashboardPage(
  dashboardHeader(title = tags$div(
    tags$img(src = "juggling_squid_cute_.jpeg", height = "40px", style = "margin-right:10px;"),'SquIG')
  ),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Upload Data", tabName = "upload", icon = icon("upload")),
      menuItem("Gating", tabName = "gating", icon = icon("project-diagram")),
      menuItem("Results", tabName = "results", icon = icon("table"))
    )
  ),
  
  dashboardBody(
    tabItems(
      # Upload tab
      tabItem(
        tabName = "upload",
        fluidRow(
          box(
            title = "Upload CSV File",
            width = 12,
            fileInput("csv_file", "Choose CSV File",
                      accept = c("text/csv", 
                                 "text/comma-separated-values,text/plain", 
                                 ".csv",
                                 ".csv.gz",
                                 "application/gzip"))
          )
        ),
        fluidRow(
          box(
            title = "Import Gating Schema",
            width = 12,
            selectInput("import_gates_yn", 
                        "Import a previous gating schema?",
                        choices = c("No" = "no", "Yes" = "yes"),
                        selected = "no"),
            conditionalPanel(
              condition = "input.import_gates_yn == 'yes'",
              fileInput("gate_metadata_file", 
                        "Upload gate_metadata.json file:",
                        accept = c("application/json", ".json")),
              actionButton("restore_gates", 
                           "Apply Gating Schema", 
                           class = "btn-success"),
              hr(),
              uiOutput("schema_validation_msg")
            )
          )
        ),
        fluidRow(
          uiOutput("data_preview_box")
        )
      ),
      
      # Gating tab
      tabItem(
        tabName = "gating",
        fluidRow(
          box(
            title = "Gate Controls",
            width = 3,
            uiOutput("x_axis_select"),
            uiOutput("y_axis_select"),
  
            selectInput("x_scale", "X Variable Scale:",
                        choices = c("Original" = "none",
                                   "log2" = "log2",
                                   "Natural log (ln)" = "ln",
                                   "log10" = "log10",
                                   "Z-score" = "zscore"),
                        selected = "none"),
  
            selectInput("y_scale", "Y Variable Scale:",
                        choices = c("Original" = "none",
                                   "log2" = "log2",
                                   "Natural log (ln)" = "ln",
                                   "log10" = "log10",
                                   "Z-score" = "zscore"),
                        selected = "none"),

            selectInput("numeric_gate_yn", "Input Numeric Values?",
                        choices = c("No" = "no", "Yes" = "yes"),
                        selected = "no"),

            # --- Histogram/slider mode ---
            conditionalPanel(
              condition = "input.numeric_gate_yn == 'no'",
              h5("X Axis Range"),
              plotOutput("x_histogram", height = "80px"),
              uiOutput("x_slider_ui"),
              h5("Y Axis Range"),
              plotOutput("y_histogram", height = "80px"),
              uiOutput("y_slider_ui")
            ),

            # --- Numeric input mode ---
            conditionalPanel(
              condition = "input.numeric_gate_yn == 'yes'",
              h5("X Gate Threshold"),
              fluidRow(
                column(5, selectInput("x_operator", NULL,
                                      choices = c(">" = "gt", ">=" = "gte", 
                                                 "<" = "lt", "<=" = "lte", 
                                                 "between" = "between"),
                                      selected = "gte")),
                column(7, numericInput("x_threshold", NULL, value = 0, step = 0.01))
              ),
              conditionalPanel(
                condition = "input.x_operator == 'between'",
                fluidRow(
                  column(5, tags$label("and <=")),
                  column(7, numericInput("x_threshold_max", NULL, value = 1, step = 0.01))
                )
              ),
              h5("Y Gate Threshold"),
              fluidRow(
                column(5, selectInput("y_operator", NULL,
                                      choices = c(">" = "gt", ">=" = "gte", 
                                                 "<" = "lt", "<=" = "lte", 
                                                 "between" = "between"),
                                      selected = "gte")),
                column(7, numericInput("y_threshold", NULL, value = 0, step = 0.01))
              ),
              conditionalPanel(
                condition = "input.y_operator == 'between'",
                fluidRow(
                  column(5, tags$label("and <=")),
                  column(7, numericInput("y_threshold_max", NULL, value = 1, step = 0.01))
                )
              ),
              br(),
              uiOutput("numeric_gate_preview_stats")
            ),

            hr(),
  
            textInput("gate_name", "Gate Name:", value = "Gate_1"),
            helpText("Names will automatically have ' Gate' appended unless they already contain 'gate'"),
            actionButton("apply_gate", "Apply Gate", class = "btn-primary"),
            helpText("Brush mode: draw a selection on the plot first. Numeric mode: set thresholds above."),
            br(),
            actionButton("reset_gate", "Clear Selection", class = "btn-warning"),
            br(), br(),
            actionButton("reset_all", "Reset All Gates", class = "btn-danger"),
            br(), br(),
            hr(),
            h4("Gate History Tree"),
            helpText("Click a gate to view/branch from that population"),
            uiOutput("gate_tree"),
            hr(),
            h4("Display Options"),
            uiOutput("gate_display_checkboxes"),
            checkboxInput("show_overlaps_only",
                          "Show overlaps only",
                          value = FALSE),
            checkboxInput("show_counts_in_legend",
                          "Show counts in legend",
                          value = TRUE),
            downloadButton("download_current_view", "Save Current View", class = "btn-info btn-sm"),
            hr(),
            h4("Current Path Statistics"),
            verbatimTextOutput("gate_stats")
          ),
          box(width = 9,
            uiOutput("current_gate_info"),
            plotOutput("gate_plot", height = 600, 
                       brush = brushOpts(
                         id = "plot_brush",
                         resetOnNew = TRUE,
                         direction = "xy",
                         stroke = "#ff0000",
                         fill = "#ff000040",
                         delay = 300,
                         delayType = "debounce"
                       ))
          )
        ),
        fluidRow(
          box(
            title = "All Gates Summary",
            width = 12,
            verbatimTextOutput("all_gates_summary")
          )
        )
      ),
      
      # Results tab
      tabItem(
        tabName = "results",
        fluidRow(
          box(
            title = "Download Options",
            width = 12,
            downloadButton("download_plots", "Download Plots (ZIP)"),
            downloadButton("download_csv", "Download Gated Data with Gate Columns (CSV)"),
            downloadButton("download_gate_metadata", "Download Gate Metadata (JSON)", class = "btn-info"),
            downloadButton("download_all", "Download Complete Package (ZIP)"),
            helpText("Complete package includes: all plots (PDF), final gated data (CSV), original data with gate columns (CSV), and gate metadata (JSON)")
          )
        ),
        fluidRow(
          box(
            title = "Final Gated Data Preview",
            width = 12,
            helpText("This shows only the cells that passed all gates. The downloadable CSV includes gate membership columns."),
            DTOutput("final_data")
          )
        )
      )
    )
  ),
  title='SquIG'
)
