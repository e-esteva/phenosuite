require(shiny)
require(stringr)
require(glue)
require(shinyjs)
require(cowplot)
require(ggplot2)
require(shinyFiles)
require(phenomenalist)
require(ggpubr)
source('/srv/shiny-server/phenomenalist/utils/provenance.R')

server=shinyServer( function(input, output, session) {
  
  # returns global cell type set:
  mydata <- reactive({

    inFile=input$PCFs

    if (is.null(inFile))
      return(NULL)

    tracker$register_input(inFile, input_id = "PCFs")

    # to account for multiple files:
    if (length(inFile) > 0) {

      print(inFile)
      
      # to account for multiple files:
      tbls <- lapply(seq(length(inFile$datapath)),function(x){
        tmp=read.csv(glue("{inFile$datapath[x]}"),row.names = 1)
        #names(tmp)=gsub('[.]{2}','+',names(tmp))
        names(tmp)=gsub('[.]{1,2}',' ',names(tmp))
        return(tmp)
      } )
      print('loaded')
      
      tbls_phen <- lapply(tbls,function(x) names(x)[names(x)!='Sample'])
      set_ = unique(unlist(tbls_phen))
      message(set_)
      mapping_=do.call('cbind',lapply(seq(length(tbls_phen)),function(x){
        tmp=rep(0,length(set_))
        tmp[match(tbls_phen[[x]],set_)]=1
        return(tmp)
      }))
      return(set_[rowSums(mapping_)==length(tbls)])
      
    }
  })
  
  # returns rbinded dataset
  mydata0 <- reactive({
    celltypes=mydata()
    inFile=input$PCFs
    
    if (is.null(inFile))
      return(NULL)
    
    
    
    
    # to account for multiple files:
    if (length(inFile) > 0) {
      
      print(inFile)
      message('inFile: ')
      message(inFile)
      message('celltypes: ')
      message(paste(celltypes,sep=','))

      # to account for multiple files:
      tbls <- do.call('rbind',lapply(seq(length(inFile$datapath)),function(x){
        tmp=read.csv(glue("{inFile$datapath[x]}"),row.names = 1)
        #names(tmp)=gsub('[.]{2}','+',names(tmp))
        names(tmp)=gsub('[.]{1,2}',' ',names(tmp))
	message(names(tmp))
        return(tmp[,c(celltypes,'Sample')])
	#return(tmp)
      } ))
      
      return(tbls)
    }
    
  })
  
  # return sample names
  mydata1 <- reactive({
    
    inFile <- input$PCFs
    if(!is.null(inFile)){
      tbls <- lapply(seq(length(inFile$datapath)),function(x){
        tmp=read.csv(glue("{inFile$datapath[x]}"),row.names = 1)
        #names(tmp)=gsub('[.]{2}','+',names(tmp))
        names(tmp)=gsub('[.]{1,2}',' ',names(tmp))
        return(tmp)
      } )
      print('loaded sample names')
      tbls_phen <- lapply(tbls,function(x) unique(x[['Sample']]))
      
      set_ = unique(unlist(tbls_phen))
      message(set_)
      return(set_)
      
      new_samples=new_group_names()
      if(!is.null(new_samples)){
        return(new_samples)
      }
    }else{
      return(NULL)
    }
    
    
  })
  
  
  
  
  
  observeEvent(input$reset_button, {js$resetClick()})
  
  
  
  
  
  
  output$samples <- renderUI({
    if(!is.null(input$PCFs)){
      samples = mydata1()
      
      message(samples)
      
      
      numSamples <- as.integer(length(unique(samples)))
      lapply(1:numSamples, function(i) {
        textInput(glue('sample_{i}'),unique(samples)[i],value = unique(samples)[i])
        
      })
    }
    
    
  })
  new_group_names=reactive({
    
      
      samples = mydata1()
      
      
      numClusters <- as.integer(length(samples))
      groups_ = sapply(seq(numClusters),function(i) glue('sample_{i}'))
      new_groups_ = sapply(groups_,function(x) input[[x]])
      message('new_groups_:')
      message(new_groups_)
      
      for(i in seq(length(groups_))){
        group=input[[glue('sample_{i}')]]
        message(group)
        message(samples[i])
        
        samples[new_groups_==group]=group
      }
      print(table(samples))
      
      return(samples)
      
      
      
    
  })
  
  observe({
    samples = new_group_names()
    updateSelectInput(inputId = 'ref_selection',choices = samples,selected = "")
    
    celltypes = mydata()
    updateSelectInput(inputId = 'celltype_to_analyze',choices = celltypes,selected = "")
    
  })
  
  tempdir=file.path(tempdir(), as.integer(Sys.time()))
  dir.create(tempdir)
  print(as.character(tempdir))
  tempdir0 = as.character(tempdir)
  tracker <- ProvenanceTracker$new("pcf_builder", session, tempdir0)
  print(tempdir0)

  output$plot=renderPlot({
    
    samples=mydata1()
    groups=new_group_names()
    
    data=mydata0()
    print(head(data))
    
    ref = input$ref_selection
    celltype=input$celltype_to_analyze
    print(ref)
    print(celltype)
    
    if(!is.null(data) & input$confirm_pcf > 0){
      tracker$capture_parameters(input)
      tracker$analysis_started()
      global_pcf=data

      for(i in seq(length(groups))){
        group.tmp=groups[i]
        global_pcf$Sample=gsub(samples[i],group.tmp,as.character(global_pcf$Sample))
      }
      print(table(global_pcf$Sample))
      p=ggviolin(global_pcf,x='Sample',y=celltype,color = 'Sample',add = 'boxplot')+geom_hline(yintercept = mean(global_pcf[,match(celltype,names(global_pcf))][global_pcf$Sample==ref]))+theme(legend.position = "none")+xlab('')+ylab('norm PCF')+ggtitle(glue('{celltype} Interactions | all versus {ref}'))+stat_compare_means(ref.group = ref)+coord_flip()+stat_summary(fun = "mean",geom = "point",color = "red")
      print(p)

      if(input$confirm_pcf > 0){
        ggsave(glue('{tempdir0}/${input$run_label}.pdf'),p)
      }
      tracker$analysis_completed()
    }
    
    
    
  })
  
  output$pcf_download <- downloadHandler(
    filename = function(){
      glue("{input$run_label}-pcf_builder-output-{Sys.Date()}.zip")
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
  
  
  
})

