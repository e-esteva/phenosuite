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
  
  
  # return cell type set:
  mydata0 <- reactive({
    	global_intersect=function(list){
  		init=list[[1]]
  		for(i in seq(2,length(list))){
    			init=intersect(init,list[[i]])
  		}
  		return(init)
	}

    
    inFile=input$log_odds
    
    if (is.null(inFile))
      return(NULL)
    
    
    
    
    # to account for multiple files:
    if (length(inFile) > 0) {
      
      
      # to account for multiple files:
      tbls <- lapply(seq(length(inFile$datapath)),function(x) read.csv(glue("{inFile$datapath[x]}"),row.names = 1))
      message('loaded')
      message(inFile$datapath)
      message(length(tbls))
      tbls_phen <- lapply(tbls,function(x) unique(row.names(x)))
      set_ = unique(unlist(tbls_phen))
      for(y in seq_along(tbls_phen)){
	      message(tbls_phen[[y]])
      }
      message(set_)
      mapping_=do.call('cbind',lapply(seq(length(tbls_phen)),function(x){
        tmp=rep(0,length(set_))
        tmp[match(tbls_phen[[x]],set_)]=1
        return(tmp)
      }))
      message(set_[rowSums(mapping_)==length(tbls)])
      }
      return(global_intersect(tbls_phen))
      #return(set_[rowSums(mapping_)==length(tbls)])
      
    })
  
  # log-odds list:
  mydata1 <- reactive({
    inFile=input$log_odds
    if(is.null(inFile)){
      return(NULL)
    }else{
      tracker$register_input(inFile, input_id = "log_odds")
      tbls <- lapply(seq(length(inFile$datapath)),function(x){
        tmp=read.csv(glue("{inFile$datapath[x]}"),row.names = 1)
        colnames(tmp)=row.names(tmp)
        return(tmp)
      } )
      return(tbls)
    }

  })
  
  
  
  
  # ref selection for unimodal view
  mydata2 <- reactive({
    return(input$ref_selection)
  })
  
  observeEvent(input$reset_button, {js$resetClick()})
  
  observe({
    celltypes = mydata0()
    updateCheckboxGroupInput(inputId = 'celltype_selection',choices = celltypes,selected = celltypes)
    
    #ref = mydata2()
    updateSelectInput(inputId = 'ref_selection',choices=celltypes,selected="")
  })
  
  output$clusters <- renderUI({
    
    celltypes=input$celltype_selection
    
    
    
    
    numClusters <- as.integer(length(celltypes))
    lapply(1:numClusters, function(i) {
      textInput(glue('cluster_{i}'),celltypes[i],value = celltypes[i])
      
    })
    
  })
  
  transform_status <- reactive({
    log_odds_ = mydata1()
    celltypes = input$celltype_selection
    req(log_odds_, length(celltypes) > 0)
    per_matrix <- lapply(log_odds_, function(x) {
      ct = intersect(celltypes, row.names(x))
      tmp = as.matrix(x[ct, ct])
      diag(tmp) = 0
      list(
        has_inf      = sum(is.infinite(tmp)) > 0,
        has_negative = any(tmp[is.finite(tmp)] < 0),
        is_nonneg    = all(tmp[is.finite(tmp)] >= 0)
      )
    })
    list(
      has_inf      = sapply(per_matrix, `[[`, "has_inf"),
      has_negative = sapply(per_matrix, `[[`, "has_negative"),
      is_nonneg    = sapply(per_matrix, `[[`, "is_nonneg")
    )
  })

  output$transform_warning <- renderUI({
    status <- transform_status()
    req(status)
    has_inf    <- status$has_inf
    is_nonneg  <- status$is_nonneg
    has_neg    <- status$has_negative

    warnings <- tagList()

    if (any(has_inf) && !all(has_inf)) {
      n_inf <- sum(has_inf)
      n_no_inf <- sum(!has_inf)
      warnings <- tagList(warnings,
        tags$div(
          style = "background-color: #fff3cd; border: 1px solid #ffc107; padding: 10px; border-radius: 5px; margin-bottom: 10px;",
          tags$strong("Scale mismatch (Inf): "),
          tags$span(glue("{n_inf} matrix(es) contain Inf values (will be transformed to 0+ range) while {n_no_inf} matrix(es) have +/- values.")),
          checkboxInput('harmonize_transform', 'Transform all matrices to 0+ range for scale harmonization', value = FALSE)
        )
      )
    }

    if (!any(has_inf) && any(is_nonneg) && any(has_neg)) {
      n_nonneg <- sum(is_nonneg)
      n_neg <- sum(has_neg)
      warnings <- tagList(warnings,
        tags$div(
          style = "background-color: #fff3cd; border: 1px solid #ffc107; padding: 10px; border-radius: 5px; margin-bottom: 10px;",
          tags$strong("Scale mismatch (0+ vs +/-): "),
          tags$span(glue("{n_nonneg} matrix(es) are 0+ (likely previously transformed) while {n_neg} matrix(es) contain negative values.")),
          checkboxInput('harmonize_transform_neg', 'Transform +/- matrices to 0+ range for scale harmonization', value = FALSE)
        )
      )
    }

    if (length(warnings) > 0) warnings else NULL
  })

  tempdir=file.path(tempdir(), as.integer(Sys.time()))
  dir.create(tempdir)
  print(as.character(tempdir))
  tempdir0 = as.character(tempdir)
  tracker <- ProvenanceTracker$new("circos_builder", session, tempdir0)
  print(tempdir0)

  output$circos <- renderPlot({
    
    log_odds_ = mydata1()
    unimodal = mydata2()
    label=input$run_label
    celltypes=input$celltype_selection
    message(paste0(celltypes,collapse=','))
    
    numClusters <- as.integer(length(celltypes))
    clusters_ = sapply(seq(numClusters),function(i) glue('cluster_{i}'))
    new_clusters_ = celltypes
    
    if(input$edit_clusters > 0){
      for(i in seq(length(clusters_))){
        message(input[[glue('cluster_{i}')]])
        message(celltypes[i])
        
        new_clusters_[new_clusters_==celltypes[i]]=input[[glue('cluster_{i}')]]
      }
      print(table(new_clusters_))
      
      print(length(new_clusters_))
      
      message(glue('new cluster names: {paste0(new_clusters_,collapse=",")}'))
      
      
    }
    
    
    if(!is.null(input$log_odds)){
      
      
      
      
      
      print('Unimodal:')
      print(unimodal)
      
      source("/srv/shiny-server/phenomenalist/utils/spatial-shiny/renderCircos-shiny_v2.R")
      
      print('log_odds_ :')
      print(log_odds_)
      discontinuity.check=sum(unlist(lapply(log_odds_,function(x){
        tmp=x[celltypes,celltypes]
	if(unimodal != ""){
		tmp=tmp[,-match(unimodal,celltypes)]=0
	}
        hit=ifelse(sum(is.infinite(as.matrix(tmp))) > 0,T,F)
        return(hit)
      } )))
      if(label != ""){
        if(discontinuity.check > 0 || isTRUE(input$harmonize_transform) || isTRUE(input$harmonize_transform_neg)){
          label=glue('{label}_transformed')
        }
      }
      
      if(input$Run > 0){

        tracker$capture_parameters(input)
        tracker$analysis_started()

        celltypes=input$celltype_selection
        message(paste0(celltypes,collapse=","))

        if(input$action == 'Integrate'){
          label=glue('{label}-integrated')
          for(i in seq(length(log_odds_))){
            log_odds=log_odds_[[i]]
            
            print(celltypes)
            print(row.names(log_odds))
            print(colnames(log_odds))
            
            log_odds=log_odds[celltypes,celltypes]
            print(log_odds)
            colnames(log_odds)=new_clusters_
            row.names(log_odds)=new_clusters_
            log_odds_[[i]]=log_odds
          }
          
          harmonize_all = isTRUE(input$harmonize_transform)
          harmonize_neg = isTRUE(input$harmonize_transform_neg)
          status = transform_status()

          log_odds.tmp = lapply(seq_along(log_odds_),function(i){
            tmp=log_odds_[[i]]
            if(status$has_inf[i] || harmonize_all){
              tmp = log(exp(tmp)+1)
            } else if(harmonize_neg && status$has_negative[i]){
              tmp = log(exp(tmp)+1)
            }
            return(tmp)
          } )
          log_odds=log_odds.tmp[[1]]
          for(i in seq(2,length(log_odds.tmp))){
            log_odds = log_odds + log_odds.tmp[[i]]
          }
          log_odds = log_odds/length(log_odds_)

          print(log_odds)
          if(unimodal != ""){
            label=glue('{label}-{unimodal}_unimodal')
            log_odds[,-match(unimodal,colnames(log_odds))]=0

          }
          write.csv(log_odds,glue('{tempdir0}/{label}.csv'))

          print(log_odds)

          was_transformed = (discontinuity.check > 0) || harmonize_all || harmonize_neg
          renderCircos(log_odds,label = label,p1=NULL,p2=NULL,out_dir=tempdir0,continuous_color_scheme = ifelse(input$color_scheme=='Continuous',T,F),scale=input$scale,discontinuity=was_transformed,label_size.cex=input$label_size,transformed=was_transformed)
          
        }
        if(input$action == 'Harmonize'){
        
	  #source("/srv/shiny-server/phenomenalist/utils/spatial-shiny/renderCircos-shiny-col_fun.R")

          inFile=input$log_odds
          
          message(glue('file name: {inFile$name}'))
          labels=do.call('rbind',strsplit(inFile$name,'[.]'))[,1]
          
          for(i in seq(length(log_odds_))){
            log_odds=log_odds_[[i]]
            log_odds=log_odds[celltypes,celltypes]
            print(log_odds)
            
            colnames(log_odds)=new_clusters_
            row.names(log_odds)=new_clusters_
            log_odds_[[i]]=log_odds
          }
          
          harmonize_all = isTRUE(input$harmonize_transform)
          harmonize_neg = isTRUE(input$harmonize_transform_neg)
          status = transform_status()

          log_odds = lapply(seq_along(log_odds_),function(i){
            tmp=log_odds_[[i]][new_clusters_,new_clusters_]
            message(unimodal)
	    message(colnames(tmp))
	    if(unimodal != ""){
              tmp[,-match(unimodal,colnames(tmp))]=0
            }
            if(status$has_inf[i] || harmonize_all){
	      tmp = log(exp(tmp)+1)
            } else if(harmonize_neg && status$has_negative[i]){
              tmp = log(exp(tmp)+1)
            }
            diag(tmp)=0
            return(tmp)
          })
          
          logOdds=do.call('rbind',log_odds)
          
	  if(input$color_scheme == 'Continuous'){
            col_fun.h = colorRamp2(c(min(logOdds[,abs(colSums(logOdds)) >0]),ifelse(discontinuity.check == 0,median(logOdds[,abs(colSums(logOdds)) >0]),0),max(logOdds[,abs(colSums(logOdds)) >0])), c( "white",'yellow' ,"red"))
          }else{
            col_fun.h = colorRamp2(c(min(logOdds[,abs(colSums(logOdds)) >0]),0,max(logOdds[,abs(colSums(logOdds)) >0])), c("blue",'white' ,"red"))
          }
          
          
          for(i in seq(length(log_odds))){
            
            label.tmp=labels[i]
            label.tmp=glue('{label.tmp}-harmonized')
            
            if(unimodal != ""){
              label.tmp=glue('{label.tmp}-{unimodal}_unimodal')
            }
            
            message(label.tmp)
            message(log_odds[[i]])
            #renderCircos(log_odds[[i]],label = label.tmp,p1=NULL,p2=NULL,out_dir=tempdir0,col_fun=col_fun.h,scale=input$scale)
	    # to adjust discontinuous values automatically to median instead of 0 for median:
	    was_transformed = (discontinuity.check > 0) || harmonize_all || harmonize_neg
	    renderCircos(log_odds[[i]],label = label.tmp,p1=NULL,p2=NULL,out_dir=tempdir0,continuous_color_scheme = ifelse(input$color_scheme=='Continuous',T,F),scale=input$scale,discontinuity=was_transformed,col.fun=col_fun.h,label_size.cex=input$label_size,transformed=was_transformed)
          }


        }
        tracker$analysis_completed()
      }
      
      
      
      
    }
    
    
    
    
    
    
    
    
  })
  output$circos_download <- downloadHandler(
    filename = function(){
      glue("{input$run_label}-circos_builder-output-{Sys.Date()}.zip")
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
