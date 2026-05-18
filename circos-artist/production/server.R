require(shiny)
require(stringr)
require(glue)
require(shinyjs)
require(shinyFiles)
require(reticulate)
require(circlize)
require(ComplexHeatmap)
require(dplyr)
require(reshape2)
require('shinycssloaders')
source('/srv/shiny-server/phenomenalist/utils/provenance.R')

## Define server logic required to parse args and launch program:

server <- function(input, output,session) {
   
  # return spatial coords file:
  mydata0 <- reactive({

    inFile <- input$file1

    if (is.null(inFile))
      return(NULL)

    tracker$register_input(inFile, input_id = "file1")
    tbl = read.csv(inFile$datapath,row.names=1)

    colnames(tbl)=row.names(tbl)

    return(tbl)
  })
  
  output$clusters <- renderUI({

    log_odds = mydata0()
    req(log_odds)

    print(dim(log_odds))
    print(colnames(log_odds))
    print(row.names(log_odds))

    celltype_selection=input$celltype_selection
    # Wait for the checkbox observer to populate the selection before
    # attempting to match; otherwise we subset to a 0-row frame and render
    # zero textInputs, which silently hides the cluster-rename UI.
    req(length(celltype_selection) > 0)

    # Historical note: the gsubs below were an attempt to reverse Shiny's
    # make.names-style encoding of checkbox values. In practice the
    # checkboxGroupInput preserves the raw celltype strings we hand it, so
    # the transforms mostly create mismatches. Try matching the raw
    # selection first and only fall back to the transformed form when the
    # raw match returns nothing.
    ct_raw   <- celltype_selection
    ct_alt   <- gsub('[.]{1}', ' ', gsub('[.]{2}', '+', celltype_selection))

    pick_rows <- function(sel) na.omit(match(sel, row.names(log_odds)))
    pick_cols <- function(sel) na.omit(match(sel, colnames(log_odds)))

    row_idx <- pick_rows(ct_raw)
    col_idx <- pick_cols(ct_raw)
    if (length(row_idx) == 0 && length(col_idx) == 0) {
      row_idx <- pick_rows(ct_alt)
      col_idx <- pick_cols(ct_alt)
    }
    req(length(row_idx) > 0, length(col_idx) > 0)

    print(colnames(log_odds))
    print(row.names(log_odds))

    log_odds = log_odds[row_idx, col_idx, drop = FALSE]

    numClusters <- as.integer(length(colnames(log_odds)))
    lapply(seq_len(numClusters), function(i) {
      textInput(glue('cluster_{i}'), colnames(log_odds)[i],
                value = colnames(log_odds)[i])
    })

  })
  # return spatial coords file path:
  mydata1 <- reactive({
    
    inFile <- input$file1
    
    if (is.null(inFile))
      return(NULL)
    
    
    return(inFile$datapath)
  })
  
  # returns celltypes
  mydata1 <- reactive({
    
    inFile <- input$file1
    
    if (is.null(inFile))
      return(NULL)
    
    
    tbl = read.csv(inFile$datapath,row.names = 1)
    colnames(tbl)=row.names(tbl)
    
    return(colnames(tbl))
  })
  
  mydata2 <- reactive({
    return(input$ref_selection)
  })
  observeEvent(input$reset_button, {js$resetClick()})
  
  observe({
    celltypes = mydata1()
    updateCheckboxGroupInput(inputId = 'celltype_selection',choices = celltypes,selected = celltypes)
    
    ref = mydata1()
    updateSelectInput(inputId = 'ref_selection',choices=ref,selected="")
  })
  
  
  tempdir=file.path(tempdir(), as.integer(Sys.time()))
  dir.create(tempdir)
  print(as.character(tempdir))
  tempdir0 = as.character(tempdir)
  tracker <- ProvenanceTracker$new("circos_artist", session, tempdir0)
  print(tempdir0)

  output$circos <- renderPlot({
    log_odds = mydata0()
    # Bail early on the very first reactive pass before mydata0() has returned —
    # without this guard the rest of the function runs on NULL and the
    # chordDiagram call at the bottom explodes with a cryptic dimnames error.
    req(log_odds)

    unimodal=mydata2()
    print('Unimodal:')
    print(unimodal)

    source('/srv/shiny-server/phenomenalist/utils/spatial-shiny/renderCircos-shiny_v2.R')

    #source('/Users/ee699/working/TRIC/Phenoptics-shiny-v1/spatial_interactions/utils/spatial-shiny/renderCircos-shiny.R')



    print(log_odds)
    # Only zero-out non-unimodal columns when unimodal actually matches a column,
    # otherwise match() returns NA and `log_odds[,-NA] = 0` throws.
    if (!is.null(unimodal) && nzchar(unimodal)) {
      uni_idx <- match(unimodal, colnames(log_odds))
      if (!is.na(uni_idx)) {
        log_odds[, -uni_idx] = 0
      }
    }
    print(log_odds)
    celltype_selection=input$celltype_selection

    # RACE GUARD: on the first reactive pass after a file upload, the
    # observe() that populates input$celltype_selection via
    # updateCheckboxGroupInput has not fired yet, so celltype_selection is
    # NULL. Subsetting log_odds by na.omit(match(NULL, ...)) yields a 0x0
    # matrix whose dimnames confuse chordDiagramFromMatrix. Bail out here
    # and wait for the next reactive cycle where the selection is populated.
    req(length(celltype_selection) > 0)

    log_odds=log_odds[na.omit(match(celltype_selection,row.names(log_odds))),na.omit(match(celltype_selection,colnames(log_odds)))]

    # Second guard: if none of the selected celltypes are actually present
    # in the log-odds matrix (e.g. the user deselected everything, or the
    # label encoding between the checkbox and the CSV diverged), log_odds
    # is now 0x0 — skip drawing rather than crashing.
    req(nrow(log_odds) > 0, ncol(log_odds) > 0)

    numClusters <- as.integer(length(colnames(log_odds)))
    clusters_ = sapply(seq(numClusters),function(i) glue('cluster_{i}'))
    new_clusters_ = as.character(colnames(log_odds))
    print(table(new_clusters_))

    if(input$edit_clusters > 0){
      for(i in seq(length(clusters_))){
        print(input[[glue('cluster_{i}')]])
        print(colnames(log_odds)[i])

        new_clusters_[new_clusters_==colnames(log_odds)[i]]=input[[glue('cluster_{i}')]]
      }
      print(table(new_clusters_))

      print(length(new_clusters_))
      print(dim(log_odds))

      colnames(log_odds)=new_clusters_
      row.names(log_odds)=new_clusters_
    }


    if(!is.null(input$file1)){
      log_odds=as.matrix(log_odds)

      label=input$run_label
      was_transformed <- FALSE

      # Stabilize -Inf / Inf unconditionally. Previously this was gated on
      # `label != ""` — meaning a user who loaded a log-odds CSV without
      # typing a run label got raw -Inf values passed into colorRamp2, which
      # emits garbage hex strings ("NAFF") that crash chordDiagramFromMatrix
      # deep in col2rgb with a misleading "undefined columns selected" error.
      # The transform now always runs; only the filename suffix is gated on
      # there being a user-supplied label to decorate.
      if (sum(is.infinite(log_odds)) > 0) {
        log_odds        <- log(exp(log_odds) + 1)
        was_transformed <- TRUE
        if (nzchar(label)) {
          label <- glue("{label}_transformed")
        }
      }

      renderCircos(log_odds,label = label,p1=NULL,p2=NULL,out_dir=NULL,continuous_color_scheme = ifelse(input$color_scheme=='Continuous',T,F),scale=input$scale,label_size.cex=input$label_size,transformed = was_transformed,self_interactions = (input$self_interactions == 'Yes'))
    }





    if(input$Run > 0){

      tracker$capture_parameters(input)
      tracker$analysis_started()

      was_transformed <- FALSE
      label=input$run_label
      if(sum(is.infinite(log_odds)) > 0){
        log_odds = log(exp(log_odds)+1)
        was_transformed <- TRUE
        if (nzchar(label)) {
          label <- glue('{label}_transformed')
        }
      }
      renderCircos(log_odds,label = label,p1=NULL,p2=NULL,out_dir=tempdir0,continuous_color_scheme = ifelse(input$color_scheme=='Continuous',T,F),scale=input$scale,label_size.cex=input$label_size,transformed = was_transformed,self_interactions = (input$self_interactions == 'Yes'))
      tracker$analysis_completed()
    }


  })
  output$circos_download <- downloadHandler(
    filename = function(){
      glue("{input$run_label}-circos_artist-output-{Sys.Date()}.zip")
    },
    content = function(file){
      
      
      print(tempdir0) 
      zip::zip(
        zipfile = file,
        files = dir(tempdir0),
        root = tempdir0
      )
    },
    contentType='application/zip'
  )
  
  
}  
