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
    
    print(dim(log_odds))
    print(colnames(log_odds))
    print(row.names(log_odds))
    
    celltype_selection=input$celltype_selection
    #celltype_selection=gsub('[+]','pos',celltype_selection)
    #celltype_selection=gsub('[.]','pos',celltype_selection)
    celltype_selection=gsub('[.]{2}','+',celltype_selection)
    celltype_selection=gsub('[.]{1}',' ',celltype_selection)
    print(celltype_selection)
    
    # colnames(log_odds)=gsub('[+]','.',colnames(log_odds))
    # row.names(log_odds)=gsub('[+]','.',row.names(log_odds))
    # 
    # colnames(log_odds)=gsub('[ ]','.',colnames(log_odds))
    # row.names(log_odds)=gsub('[ ]','.',row.names(log_odds))
    
    print(colnames(log_odds))
    print(row.names(log_odds))
    
    print(log_odds[na.omit(match(celltype_selection,row.names(log_odds))),na.omit(match(celltype_selection,colnames(log_odds)))])
    log_odds=log_odds[na.omit(match(celltype_selection,row.names(log_odds))),na.omit(match(celltype_selection,colnames(log_odds)))]
    
    
    
    numClusters <- as.integer(length(colnames(log_odds)))
    lapply(1:numClusters, function(i) {
      textInput(glue('cluster_{i}'),colnames(log_odds)[i],value = colnames(log_odds)[i])
      
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
    unimodal=mydata2()
    print('Unimodal:')
    print(unimodal)
    
    source('/srv/shiny-server/phenomenalist/utils/spatial-shiny/renderCircos-shiny_v2.R')

    #source('/Users/ee699/working/TRIC/Phenoptics-shiny-v1/spatial_interactions/utils/spatial-shiny/renderCircos-shiny.R')
    
    
    
    print(log_odds)
    if(unimodal != ""){
      log_odds[,-match(unimodal,colnames(log_odds))]=0
    }
    print(log_odds)
    celltype_selection=input$celltype_selection
    
    
    
    
    
    log_odds=log_odds[na.omit(match(celltype_selection,row.names(log_odds))),na.omit(match(celltype_selection,colnames(log_odds)))]
    

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
      if(label != ""){
        
        #if(sum(is.infinite(log_odds)) > 0){
        #  log_odds = log(exp(log_odds)+1)
        #  label=glue('{label}_transformed')
        #}
        
	if (sum(is.infinite(log_odds)) > 0) {
  		log_odds        <- log(exp(log_odds) + 1)
  		label           <- glue("{label}_transformed")
  		was_transformed <- TRUE
	}  
        
      }
      renderCircos(log_odds,label = label,p1=NULL,p2=NULL,out_dir=NULL,continuous_color_scheme = ifelse(input$color_scheme=='Continuous',T,F),scale=input$scale,label_size.cex=input$label_size,transformed = was_transformed)
    }
    
    
    
    
    
    if(input$Run > 0){

      tracker$capture_parameters(input)
      tracker$analysis_started()

      was_transformed <- FALSE
      label=input$run_label
      if(sum(is.infinite(log_odds)) > 0){
        log_odds = log(exp(log_odds)+1)
        label=glue('{label}_transformed')
	was_transformed <- TRUE
      }
      renderCircos(log_odds,label = label,p1=NULL,p2=NULL,out_dir=tempdir0,continuous_color_scheme = ifelse(input$color_scheme=='Continuous',T,F),scale=input$scale,label_size.cex=input$label_size,transformed = was_transformed)
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
