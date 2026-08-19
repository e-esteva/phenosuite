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
require(shinycssloaders)
source('/srv/shiny-server/phenomenalist/utils/provenance.R')

server <- function(input, output,session) {
   
  # return spatial coords file:
  mydata0 <- reactive({

    inFile <- input$file1

    if (is.null(inFile))
      return(NULL)

    tracker$register_input(inFile, input_id = "file1")
    tbl = read.csv(inFile$datapath)

    return(tbl)
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
    
    
    tbl = read.csv(inFile$datapath)
    
    return(unique(tbl$cluster))
  })
  
  observeEvent(input$reset_button, {js$resetClick()})
  
  observe({
    celltypes = mydata1()
    updateCheckboxGroupInput(inputId = 'celltype_selection',choices = celltypes,selected = celltypes)
  })
  
  tempdir=file.path(tempdir(), as.integer(Sys.time()))
  dir.create(tempdir)
  print(as.character(tempdir))
  tempdir0 = as.character(tempdir)
  tracker <- ProvenanceTracker$new("spatial_interactions", session, tempdir0)
  print(tempdir0)



  output$plots <- renderPlot({
    spatial_obj=mydata0()
    # action button:
    run=input$Run
    print(run)
    if(run > 0){
      tracker$capture_parameters(input)
      tracker$analysis_started()
      #text input:
      label=input$run_label
      # instrument resolution:
      resolution=input$instr_res
      # inner radius:
      p1 = input$radius_i
      # outer radius:
      p2 = input$radius_o
      
      
      
      
      
      source('/srv/shiny-server/phenomenalist/spatial_interactions/utils/spatial-shiny/renderCircos-shiny.R')
      source('/srv/shiny-server/phenomenalist/spatial_interactions/utils/spatial-shiny/genHeatmapVolcanoPlot.R')
      
      source_python('/srv/shiny-server/phenomenalist/spatial_interactions/utils/spatial-shiny/pwlo-es-pt.py')
      source("/srv/shiny-server/phenomenalist/spatial_interactions/utils/spatial-shiny/gen-Heatmaps.R")
      
      print(resolution)
      print(head(spatial_obj))
      pd=import('pandas')
      out_dir=glue('{tempdir0}/{label}')
      dir.create(out_dir)
      
      
      withProgress(message = "Running Spatial Interactions Analysis", value = 0, {
        incProgress(1 / 3, detail = "Computing log-odds")
        
        compute_significance=ifelse(input$compute_significance=='Yes',T,F)
        
        metadata=pairwise_logOdds(spatial_obj=pd$DataFrame(spatial_obj),out_dir=out_dir,label=label,resolution=resolution,p1=p1,p2=p2,compute_effect_size = compute_significance)
        
        print(metadata)
        print(class(metadata))
        
        # subset to selection for visualization:
        celltype_selection=input$celltype_selection
        
        print(celltype_selection)
        
        
        if(length(celltype_selection) > 1 && compute_significance){
          incProgress(1 / 3, detail = "Generating Heatmap and Volcano plots")
          
          metadata=metadata[celltype_selection,celltype_selection]
          
          # add heatmap / volcano plot code
          
          # effect size / incidence probability heatmaps:
          genHM_v3(celltype_selection=celltype_selection,out_dir=out_dir)
          
          
        }else{
          # unidirectional case
          metadata[,-match(celltype_selection,colnames(metadata))]=0
        }
        incProgress(1 / 3, detail = "Generating Circos")
        
        renderCircos(metadata,label = label,p1=p1,p2=p2,out_dir=out_dir)
        incProgress(1 / 3, detail = "Done")

      })
      tracker$analysis_completed()
    }

  })
  
  output$spatial_dynamics_download <- downloadHandler(
    filename = function(){
      paste("spatial-dynamics-output-", Sys.Date(), ".zip", sep = "")
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
