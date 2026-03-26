require(shiny)
require(stringr)
require(glue)
require(shinyjs)
require(shinyFiles)
require(phenomenalist)
require(openai)
require(cowplot)

source('/srv/shiny-server/phenomenalist/utils//RunPhenomenalist-shiny/RunPhenomenalist-shiny.R')
source('/srv/shiny-server/phenomenalist/utils/RunPhenomenalist-shiny/phenomenalist-utils-shiny.R')

server=shinyServer( function(input, output, session) {
  
  
  
  #volumes <- c(Home = fs::path_home(), "R Installation" = R.home(), getVolumes()())
  
  # by setting `allowDirCreate = FALSE` a user will not be able to create a new directory
  #shinyDirChoose(input, "directory", roots = volumes, session = session, restrictions = system.file(package = "base"), allowDirCreate = TRUE)
  
  # return spe object:
  mydata0 <- reactive({
    
    inFile <- input$file1
    
    if (is.null(inFile))
      return(NULL)
    
    if(!is.null(inFile)){    
    	spe = readRDS(inFile$datapath)
    
    	return(spe)
    }
  })
  # return clustering columns
  mydata1 <- reactive({
    
    inFile <- input$file1
    
    if (is.null(inFile))
      return(NULL)
    
    if(!is.null(inFile)){
    	spe = readRDS(inFile$datapath)
    
    	return(str_subset(names(colData(spe)),'cluster'))
    }
  })
  
  #output$directorypath <- renderPrint({
  #  if (is.integer(input$directory)) {
  #    cat("No directory has been selected")
  #  } else {
  #    parseDirPath(volumes, input$directory)
      
  #  }
  #})
  
  observe({
    clustering_resolutions = mydata1()
    updateCheckboxGroupInput(inputId = 'rb0',choices = clustering_resolutions,selected = "")
    
  })
  
  observeEvent(input$reset_button, {js$resetClick()})
  
  
  
  # heatmap of selected clusters
  output$plots_ai <- renderPlot({
  
    spe = mydata0()
    if(!is.null(spe) & !is.null(input$rb0)){
    	group = input$rb0
    	plot_heatmap.mod(x=spe,group_by = group,out_dir = NULL,size.row = 8,size.col = 8)
    }
  })
  
  
  
  
  
  #tempdir=file.path(tempdir(), as.integer(Sys.time()))
  tempdir=glue('/apps/home/rtmp/{session$token}')
  dir.create(tempdir)
  message(as.character(tempdir))
  tempdir0 = as.character(tempdir)
  message(tempdir0)  
  
  # assigns text inputs, makes new cluster labels, generates new heatmap
  output$plots_ai2 <- renderPlot({
    run = input$run
    if(run > 0){
      spe = mydata0()
      group = input$rb0
      
      
      
      new_clusters_ = as.character(spe[[group]])
      new_clusters_ = gsub('[+]','pos',new_clusters_)
      new_clusters_ = gsub('-','_',new_clusters_)

      assays_spe = assays(spe)
      exprs=assays_spe$exprs
      colnames(exprs)=new_clusters_
      
      explanations=list()
      usage_report=list()
      message(glue('running {input$model}'))
      clusters_ = unique(colnames(exprs))
      message('original clusters: ')
      message(paste(unique(clusters_),sep=','))

      withProgress(message = 'Running Analysis', value = 0,{
        incProgress(1/3, detail = 'Generating new annotations')
        for(i in clusters_){
          
          incProgress(0.001, detail = glue('cluster {i}'))
          
          if(length(grep(i,colnames(exprs))) > 1){
            cluster_i.avgExp=rowMeans(exprs[,colnames(exprs)==i])
          }else{
            cluster_i.avgExp=exprs[,colnames(exprs)==i]
          }
          if(input$prompt_alg == "v1"){
            # v1:
            map_vals=sapply(seq(length(cluster_i.avgExp)),function(x) paste0(c(names(cluster_i.avgExp)[x],cluster_i.avgExp[[x]]),collapse = ":"))
            
            prompt=glue('what celltype is described by {paste0(c(map_vals[cluster_i.avgExp >quantile(cluster_i.avgExp,.95) ],map_vals[cluster_i.avgExp < quantile(cluster_i.avgExp,.05)]),collapse = ",")}?')
            print(prompt)
            
            
            content_ =glue("{prompt} Give me a 3 word response")
            print(content_)
          }else{
            markers=do.call('rbind',strsplit(names(cluster_i.avgExp),'_Cytoplasm|_Nucleus_'))[,1]
            map_vals=sapply(seq(length(cluster_i.avgExp)),function(x) paste0(c(markers[x],ifelse(cluster_i.avgExp[[x]] > 0,'+','-')),collapse = ":"))
            prompt=glue('what celltype is described by {paste0(c(map_vals[cluster_i.avgExp >quantile(cluster_i.avgExp,.9) ],map_vals[cluster_i.avgExp < quantile(cluster_i.avgExp,.05)]),collapse = ",")}?')
            print(prompt)
            content_ = glue("{prompt} This is a {input$tissue}. Give me a 3 word response.Give two choices.Explain. Format response like: Choice X: choice; Explanation X: explanation.")
            print(content_)
            
            annotations_ = function(response){
              nl_delim=unlist(strsplit(response,'\n'))
              choice_idx=grep('Choice',nl_delim)
              explanation_idx=grep("Explanation",nl_delim)
              format_1=ifelse(all.equal(choice_idx,explanation_idx)==TRUE,TRUE,FALSE)
              if(format_1){
                nl_delim=nl_delim[choice_idx]
                if(length(grep('Choice [0-9]:',nl_delim)) != 0){
                  choice_delim=do.call('rbind',strsplit(nl_delim,'Choice [0-9]:'))[,2]
                }else{
                  choice_delim=do.call('rbind',strsplit(nl_delim,'Choice [A-Z]:'))[,2]
                }
                
                annotations.tmp = do.call('rbind',strsplit(choice_delim,'; Explanation'))[,1]
                annotations=do.call('rbind',strsplit(annotations.tmp,'^[ ]'))[,2]
                
                explanations.tmp=do.call('rbind',strsplit(choice_delim,'; Explanation'))[,2]
              }else{
                nl_delim.choice=nl_delim[choice_idx]
                choice_delim=do.call('rbind',strsplit(nl_delim.choice,'Choice [0-9]:'))[,2]
                annotations=do.call('rbind',strsplit(choice_delim,'^[ ]'))[,2]
                
                nl_delim.explanation=nl_delim[explanation_idx]
                explanations.tmp=do.call('rbind',strsplit(nl_delim.explanation,'Explanation'))[,2]
                explanations.tmp=do.call('rbind',strsplit(explanations.tmp,'^[ ]'))[,2]
              }
              
              return(list(annotations,explanations.tmp))
            }
          }
          
          
          
          
          
          
          
          
          
          
          
          
          response=try(create_chat_completion(
            model = input$model,
            messages = list(
              list(
                "role" = "system",
                "content" = "You are an expert immunologist."
              ),
              list(
                "role" = "user",
                "content" = content_
              )
            ),temperature = 0,openai_api_key = input$API_KEY
          ))
          message(response)
          message(glue('{response$choices$message.content}'))
	  message(class(response))
          
          
          if(class(response)!='try-error'){
            usage_report[[length(usage_report)+1]]=data.frame(response$usage)
            if(input$prompt_alg != 'v1'){
              annotation_tmp=paste0(annotations_(response$choices$message.content)[[1]],collapse = ", ")
	      annotation_tmp=gsub('-','_',annotation_tmp)
              new_clusters_[grep(glue('^{i}$'),new_clusters_)]=glue('{i}-{annotation_tmp}')
              explanations[[length(explanations)+1]]=paste0(annotations_(response$choices$message.content)[[2]],collapse = ", ")
            }else{
              # v1:
	      annotation_tmp=response$choices$message.content
	      annotation_tmp=gsub('-','_',annotation_tmp)
              new_clusters_[grep(glue('^{i}$'),new_clusters_)]=glue('{i}-{annotation_tmp}')
            }
            
            
            
          }else{
            response=try(create_chat_completion(
              model = input$model,
              messages = list(
                list(
                  "role" = "system",
                  "content" = "You are an expert immunologist."
                ),
                list(
                  "role" = "user",
                  "content" = content_
                )
              ),temperature = 0,
		openai_api_key = input$API_KEY
            ))
            message(response)
            message(class(response))
            if(class(response)!='try-error'){
              usage_report[[length(usage_report)+1]]=data.frame(response$usage)
              if(input$prompt_alg != 'v1'){
                annotation_tmp=paste0(annotations_(response$choices$message.content)[[1]],collapse = ", ")
                annotation_tmp=gsub('-','_',annotation_tmp)
	        new_clusters_[grep(glue('^{i}$'),new_clusters_)]=glue('{i}-{annotation_tmp}')
                explanations[[length(explanations)+1]]=paste0(annotations_(response$choices$message.content)[[2]],collapse = ", ")
              }else{
                # v1:
		annotation_tmp=response$choices$message.content
              	annotation_tmp=gsub('-','_',annotation_tmp)
              	new_clusters_[grep(glue('^{i}$'),new_clusters_)]=glue('{i}-{annotation_tmp}')
                
              }
              
            }
          }
          
          
          
          
          
        }
	message('new clusters:')
        message(paste(unique(new_clusters_),sep=','))
	message('generating usage report ...')
        usage_report_ = colSums(do.call('rbind',usage_report))
        
        
        incProgress(1/3, detail = 'Downloading new annotations')
        new_clusters_clean = do.call('rbind',strsplit(new_clusters_,'[-]'))[,2]
        message(paste(unique(new_clusters_clean),sep=','))

        annotation.count=length(grep('annotated_clusters',names(colData(spe))))
        if(annotation.count == 0){
          spe[['annotated_clusters_marked']]=new_clusters_
	  spe[['annotated_clusters']]=new_clusters_clean
          group = 'annotated_clusters'
        }else{
          new.anno=as.character(glue('annotated_clusters_v{annotation.count+1}'))
          spe[[new.anno]]=new_clusters_clean
	  new.anno.marked=glue('{new.anno}_marked')
	  spe[[new.anno.marked]]=new_clusters_

          group = new.anno
        }
        
        
        #tempdir=parseDirPath(volumes, input$directory)
        #print(as.character(tempdir))
        #tempdir0 = as.character(tempdir)
        #print(tempdir0)
        
        write.csv(usage_report_,glue('{tempdir0}/usage-report-v{annotation.count+1}.csv'))
        if(length(explanations) > 0){
          
          explanations_df = data.frame(do.call('rbind',explanations))
          names(explanations_df)='Explanations'
          write.csv(explanations_df,glue('{tempdir0}/cluster-explanations-report-v{annotation.count+1}.csv'))
        }
        
        # cluster annotation results:
        cluster_annotations=data.frame('Clusters'=unique(spe[[input$rb0]]),'Annotations'=do.call('rbind',strsplit(unique(new_clusters_),'[-]'))[,2])
        write.csv(cluster_annotations,glue('{tempdir0}/cluster-annotations-v{annotation.count+1}.csv'))
        generate_colors <- function(n) {
  		if(n <= 102) {
    			# Use standard palette
    			return(rainbow(n))
  		} else {
    		# Generate more colors
    		hues <- seq(0, 1, length.out = n + 1)[1:n]
    		colors <- hsv(hues, s = 0.8, v = 0.8)
    		return(colors)
  		}
	}
        plot_heatmap.mod(x=spe,group_by = group,out_dir = tempdir0)
	source('/srv/shiny-server/phenomenalist/utils/plot-scatter.R')
        plot_dr.mod(spe,dr='UMAP',color_by = group,out_dir = tempdir0,h = 20,w = 20)
	message('generating spatial plot')
        plot_spatial.mod(spe,color_by = group,out_dir = tempdir0,h = 20,w = 20,colors=generate_colors(length(unique(spe[[group]]))))
        saveRDS(spe,glue('{tempdir0}/spe.rds'))
        spatial_coords = data.frame(spatialCoords(spe))
        spatial_coords$cluster = new_clusters_
        write.csv(spatial_coords,glue('{tempdir0}/annotated_spatial_coords-v{annotation.count+1}.csv'))
        incProgress(1/3, detail = 'Done')
      })
    }
    
  })
  
  output$plots_ai3 <- renderPlot({
    #tempdir=parseDirPath(volumes, input$directory)
    if(input$run > 0){
      spe=readRDS(glue('{tempdir0}/spe.rds'))
      
      annotation.count=length(grep('annotated_clusters',names(colData(spe))))
      if(annotation.count == 1){
        
        group = 'annotated_clusters'
      }else{
        
        annos_v=max(na.omit(as.numeric(do.call('rbind',strsplit(names(colData(spe))[grep('annotated_clusters',names(colData(spe)))],'_v'))[,2])))
        new.anno=as.character(glue('annotated_clusters_v{annos_v}'))
        
        group = new.anno
      }
      plot_heatmap.mod(x=spe,group_by = group,out_dir = NULL,size.row = 8,size.col = 8)
    }
    
    
  })
  output$phenomenalist_download <- downloadHandler(
    filename = function(){
      paste("phenomenalist-auto_phenotpying-output-", Sys.Date(), ".zip", sep = "")
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
