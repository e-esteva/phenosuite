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
                                 "application/gzip")),
	      uiOutput("previous_gates_detected")
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
  
  		# NEW: Transformation controls
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
  
  		# NEW: X axis histogram and slider
  		h5("X Axis Range"),
  		plotOutput("x_histogram", height = "80px"),
  		uiOutput("x_slider_ui"),
  
  		# NEW: Y axis histogram and slider
  		h5("Y Axis Range"),
  		plotOutput("y_histogram", height = "80px"),
  		uiOutput("y_slider_ui"),
  
  		hr(),
  
  		textInput("gate_name", "Gate Name:", value = "Gate_1"),
  		helpText("Names will automatically have ' Gate' appended unless they already contain 'gate'"),
  		actionButton("apply_gate", "Apply Gate with Lasso", class = "btn-primary"),
  		helpText("Draw a lasso selection on the plot, then click 'Apply Gate'"),
  		br(),
  		actionButton("reset_gate", "Clear Lasso", class = "btn-warning"),
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
	    ) ,
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
