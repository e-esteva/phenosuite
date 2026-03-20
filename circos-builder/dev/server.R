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

	global_union=function(list){
  		init=list[[1]]
  		for(i in seq(2,length(list))){
    			init=union(init,list[[i]])
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
      message(input$celltype_operation)
      if(input$celltype_operation=='Intersect'){
      	return(global_intersect(tbls_phen))
      }else{
	return(global_union(tbls_phen))
      }
      #return(set_[rowSums(mapping_)==length(tbls)])
      
    })
  
  # log-odds list:
  mydata1 <- reactive({
    inFile=input$log_odds
    if(is.null(inFile)){
      return(NULL)
    }else{
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
  
  tempdir=file.path(tempdir(), as.integer(Sys.time()))
  dir.create(tempdir)
  print(as.character(tempdir))
  tempdir0 = as.character(tempdir)
  print(tempdir0)
  
  output$circos <- renderPlot({
    
    log_odds_ = mydata1()
    unimodal = mydata2()
    label=input$run_label
    celltypes=input$celltype_selection
    message(paste0(celltypes,collapse=','))
    # if union selected, augument log-odds matrices to account for new celltypes
    if(input$celltype_operation == 'Union'){
	log_odds_ = lapply(log_odds_,function(x){
		tmp = x
		celltype.idx = match(celltypes,row.names(tmp))
		missing.idx = attr(na.omit(celltype.idx),'na.action')
		if(!is.null(missing.idx)){
			for(y in seq_along(celltypes[missing.idx])){
				tmp = cbind(tmp,rep(-Inf,dim(tmp)[1]))
				tmp = rbind(tmp,rep(-Inf,dim(tmp)[2]))
				colnames(tmp)[dim(tmp)[2]]=celltypes[missing.idx][y]
				row.names(tmp)[dim(tmp)[1]]=celltypes[missing.idx][y]
			}
		}
		return(tmp)
				   
    })
    }
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
    
    edits=data.frame(cbind(celltypes,new_clusters_))
    names(edits)=c('Names In','Names Out')

    write.csv(edits, glue('{tempdir0}/{label}-edits.csv'))

    if(!is.null(input$log_odds)){
      
      
      
      
      
      print('Unimodal:')
      print(unimodal)
      
      source("/srv/shiny-server/Phenoptics-Menu/utils/spatial-shiny/renderCircos-shiny.R")
      
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
        if(discontinuity.check > 0){
          label=glue('{label}_transformed')
        }
      }
      
      if(input$Run > 0){
        
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
          
          log_odds.tmp = lapply(log_odds_,function(x){
            tmp=x
            if(discontinuity.check > 0){
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
          
          renderCircos(log_odds,label = label,p1=NULL,p2=NULL,out_dir=tempdir0,continuous_color_scheme = ifelse(input$color_scheme=='Continuous',T,F),scale=input$scale,discontinuity=ifelse(discontinuity.check > 0,T,F),label_size.cex=input$label_size)
          
        }
        if(input$action == 'Harmonize'){
        
	  #source("/srv/shiny-server/Phenoptics-Menu/utils/spatial-shiny/renderCircos-shiny-col_fun.R")

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
          
          log_odds = lapply(log_odds_,function(x){
            tmp=x[new_clusters_,new_clusters_]
            message(unimodal)
	    message(colnames(tmp))
	    if(unimodal != ""){
              tmp[,-match(unimodal,colnames(tmp))]=0
            }
            #if(sum(is.infinite(as.matrix(tmp)) > 0)){
            if(discontinuity.check > 0){
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
	    renderCircos(log_odds[[i]],label = label.tmp,p1=NULL,p2=NULL,out_dir=tempdir0,continuous_color_scheme = ifelse(input$color_scheme=='Continuous',T,F),scale=input$scale,discontinuity=ifelse(discontinuity.check > 0,T,F),col.fun=col_fun.h,label_size.cex=input$label_size)
          }
          
          
        }
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
