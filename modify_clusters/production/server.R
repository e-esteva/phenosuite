require(shiny)
require(stringr)
require(glue)
require(shinyjs)
require(cowplot)
require(ggplot2)
require(shinyFiles)
require(phenomenalist)
require(Seurat)
require(tidyverse)

# Modify ID
source('/srv/shiny-server/Phenoptics-Menu/utils//RunPhenomenalist-shiny/RunPhenomenalist-shiny.R')
source('/srv/shiny-server/Phenoptics-Menu/utils/RunPhenomenalist-shiny/phenomenalist-utils-shiny.R')
source('/srv/shiny-server/Phenoptics-Menu/utils/heatmap-by-cluster.R')
source('/srv/shiny-server/phenomenalist/utils/provenance.R')

server=shinyServer( function(input, output, session) {
  
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
        if('Seurat' %in% class(spe)){
		return(str_subset(names(spe@meta.data),'cluster|classification|annotation|celltype'))
	}else{
		return(str_subset(names(colData(spe)),'cluster|classification|annotation|celltype'))
	}
    	
    }
  })
  
  seurat_markers <- reactive({
  	spe = mydata0()
	clustering_resolutions = input$rb0
  	if (!is.null(spe)){
    		if('Seurat' %in% class(spe)){
			assay_classes = lapply(spe@assays, class)
        		if('Assay5' %in% assay_classes){
				message('Assay5 detected (plot 1)')
			}
			Idents(spe)=clustering_resolutions
        		m = FindAllMarkers(spe,only.pos=T)
        		m_ = m %>% group_by(cluster) %>% top_n(wt=avg_log2FC,n=5)
			return(m_)
    }
  }
 })
  observe({
    clustering_resolutions = mydata1()
    updateCheckboxGroupInput(inputId = 'rb0',choices = clustering_resolutions,selected = "")
    
  })
  
  #output$directorypath <- renderPrint({
  #  if (is.integer(input$directory)) {
  #    cat("No directory has been selected")
  #  } else {
  #    parseDirPath(volumes, input$directory)
      
  #  }
  #})
  
  observeEvent(input$reset_button, {js$resetClick()})
  
  output$plots <- renderPlot({
    spe = mydata0()
    group = input$rb0
    if('Seurat' %in% class(spe)){
	if(!is.null(spe) & !is.null(group)){
		#obj.sc = as.SingleCellExperiment(spe)
		#spe=SpatialExperiment(assays=assays(obj.sc),rowData=rowData(obj.sc),colData=colData(obj.sc),reducedDims=reducedDims(obj.sc))
                #message(glue('spe nrows: {nrow(spe)} | spe ncols: {ncol(spe)}'))
		if(nrow(spe) > 500){
                        m_ = seurat_markers()
			message('seurat obj detected, using getheatmap by cluster')
			getHeatmapByCluster(spe,celltypes.df=NULL,genes=m_$gene,out=NULL,stat='mean',assay='RNA',cluster_col=group)     
			#plot_heatmap.mod(x=spe,group_by = group,out_dir = NULL,features=m_$gene)

	      	}else{
		 	plot_heatmap.mod(x=spe,group_by = group,out_dir = NULL)
 		}
	 }
	
    }else{
    	if(!is.null(spe) & !is.null(group)){
    		plot_heatmap.mod(x=spe,group_by = group,out_dir = NULL)
    	}
    }
    
  })
  
  
  output$clusters <- renderUI({
    clustering_resolutions = input$rb0
    spe = mydata0()   
    #print(clustering_resolutions)
    #print(unique(spe[[clustering_resolutions]])[length(unique(spe[[clustering_resolutions]]))])
    #print(do.call('rbind',strsplit(as.character(unique(spe[[clustering_resolutions]])[4]),'C[0-9]{1,2}-')))
    
    if(!is.null(spe) & !is.null(clustering_resolutions)){
	if('Seurat' %in% class(spe)){
		cluster_set=unique(spe[[clustering_resolutions]][[1]])
    		numClusters <- as.integer(length(cluster_set))
    		lapply(1:numClusters, function(i) {
      		textInput(glue('cluster_{i}'),as.character(cluster_set[i]),value = ifelse(dim(do.call('rbind',strsplit(as.character(cluster_set[i]),'C[0-9]{1,2}-')))[2] > 1,unlist(strsplit(do.call('rbind',strsplit(as.character(cluster_set[i]),'C[0-9]{1,2}-'))[,2],'[.]$')),as.character(cluster_set[i])))
      
    	})
  	}else{
		if(!is.null(spe) & !is.null(clustering_resolutions)){
        		numClusters <- as.integer(length(unique(spe[[clustering_resolutions]])))
        		lapply(1:numClusters, function(i) {
        		textInput(glue('cluster_{i}'),unique(spe[[clustering_resolutions]])[i],value = ifelse(dim(do.call('rbind',strsplit(as.character(unique(spe[[clustering_resolutions]])[i]),'C[0-9]{1,2}-')))[2] > 1,unlist(strsplit(do.call('rbind',strsplit(as.character(unique(spe[[clustering_resolutions]])[i]),'C[0-9]{1,2}-'))[,2],'[.]$')),unique(spe[[clustering_resolutions]])[i]))

    			})

  		}
  	}	
    }
  })
  
  
  output$plots_new2 <- renderPlot({
    
    spe = mydata0()
    group = input$rb0
        
    if(!is.null(spe) & !is.null(group)){
	if('Seurat' %in% class(spe)){
		numClusters <- as.integer(length(unique(spe[[group]][[1]])))
        	clusters_ = sapply(seq(numClusters),function(i) glue('cluster_{i}'))
        	new_clusters_ = as.character(spe[[group]][[1]])
        	print(table(new_clusters_))

        	for(i in seq(length(clusters_))){
                	message(input[[glue('cluster_{i}')]])
                	message(unique(spe[[group]][[1]])[i])
                	#new_clusters_[grep(glue('^{unique(spe[[group]])[i]}$'),new_clusters_)]=input[[glue('cluster_{i}')]]
                	new_clusters_[new_clusters_==unique(spe[[group]][[1]])[i]]=input[[glue('cluster_{i}')]]
        	}
        	print(table(new_clusters_))

		annotation.count=length(grep('annotated_clusters',names(spe@meta.data)))
	}else{
		numClusters <- as.integer(length(unique(spe[[group]])))
        	clusters_ = sapply(seq(numClusters),function(i) glue('cluster_{i}'))
        	new_clusters_ = as.character(spe[[group]])
        	print(table(new_clusters_))

        	for(i in seq(length(clusters_))){
                	print(input[[glue('cluster_{i}')]])
                	print(unique(spe[[group]])[i])
                	#new_clusters_[grep(glue('^{unique(spe[[group]])[i]}$'),new_clusters_)]=input[[glue('cluster_{i}')]]
                	new_clusters_[new_clusters_==unique(spe[[group]])[i]]=input[[glue('cluster_{i}')]]
        	}
        	print(table(new_clusters_))    

		annotation.count=length(grep('annotated_clusters',names(colData(spe))))
	}
    	if(annotation.count == 0){
      		spe[['annotated_clusters']]=new_clusters_
      		group = 'annotated_clusters'
    	}else{
      		new.anno=as.character(glue('annotated_clusters_v{annotation.count+1}'))
      		spe[[new.anno]]=new_clusters_
      		group = new.anno
    	}
    	#print(group)
    	print(table(spe[[group]]))
	if('Seurat' %in% class(spe)){
		#obj.sc = as.SingleCellExperiment(spe)
		#spe=SpatialExperiment(assays=assays(obj.sc),rowData=rowData(obj.sc),colData=colData(obj.sc),reducedDims=reducedDims(obj.sc))
   		#message(glue('spe nrows: {nrow(spe)} | spe ncols: {ncol(spe)}'))
		message(class(spe))
		if(nrow(spe) > 500){
			m_ = seurat_markers()
			message(head(m_$gene))
			getHeatmapByCluster(spe,celltypes.df=NULL,genes=m_$gene,out=NULL,stat='mean',assay='RNA',cluster_col=group) 
			#	plot_heatmap.mod(x=spe,group_by = group,out_dir = NULL,features=m_$gene)
		}else{
			plot_heatmap.mod(x=spe,group_by = group,out_dir = NULL)
		}
	}else{
    		plot_heatmap.mod(x=spe,group_by = group,out_dir = NULL)
    
	}
    }
    
    
  })
  
  
  tempdir=file.path(tempdir(), as.integer(Sys.time()))
  dir.create(tempdir)
  print(as.character(tempdir))
  tempdir0 = as.character(tempdir)
  tracker <- ProvenanceTracker$new("modify_clusters", session, tempdir0)
  print(tempdir0)

  output$plots_new <- renderPlot({
    spe = mydata0()
    group = input$rb0
    
    
    
    
    render = input$render
    if(render==1){
      
      withProgress(message = 'Running Analysis', value = 0,{
        incProgress(1/3, detail = 'Generating new annotations')
        if ('Seurat' %in% class(spe)){
		numClusters <- as.integer(nrow(unique(spe[[group]])))
		clusters_ = sapply(seq(numClusters),function(i) glue('cluster_{i}'))
        	new_clusters_ = as.character(spe[[group]][[1]])
	}else{
        	numClusters <- as.integer(length(unique(spe[[group]])))
		clusters_ = sapply(seq(numClusters),function(i) glue('cluster_{i}'))
        	new_clusters_ = as.character(spe[[group]])
	}
        
        if ('Seurat' %in% class(spe)){

		for(i in seq(length(clusters_))){
                        new_clusters_[new_clusters_==unique(spe[[group]][[1]])[i]]=input[[glue('cluster_{i}')]]
                }
	}else{
		for(i in seq(length(clusters_))){
          		#new_clusters_[grep(glue('^{unique(spe[[group]])[i]}$'),new_clusters_)]=input[[glue('cluster_{i}')]]
          		new_clusters_[new_clusters_==unique(spe[[group]])[i]]=input[[glue('cluster_{i}')]]
        	}

	}
        print(table(new_clusters_))
	if('Seurat' %in% class(spe)){
		annotation.count=length(grep('annotated_clusters',names(spe@meta.data)))
	}else{
		annotation.count=length(grep('annotated_clusters',names(colData(spe))))
	}
        if(annotation.count == 0){
          spe[['annotated_clusters']]=new_clusters_
          group = 'annotated_clusters'
        }else{
          new.anno=as.character(glue('annotated_clusters_v{annotation.count+1}'))
          spe[[new.anno]]=new_clusters_
          group = new.anno
        }
        
        message(table(spe[[group]]))
        
        prepare_pcf_inputs=function(spe,out_dir,group){
          library(glue)
          message(glue('generating pcf-generation inputs for {group}'))
          
          spatial_obj=data.frame(spatialCoords(spe))
          spatial_obj=cbind(row.names(spatial_obj),spatial_obj)
          
          clust.tmp=spe[[group]]
          spatial_obj$Phenotype=clust.tmp
          names(spatial_obj)=c('Sample Name','Cell X Position','Cell Y Position','Phenotype')
          if(length(unique(spe$Classifier_Label)) == 0){
	    if('Analysis_Region' %in% names(colData(spe))){
            	spatial_obj$`Tissue Category`=spe$Analysis_Region
	    }else{
		    spatial_obj$`Tissue Category`=spe$sample_id
		}

          }else{
            spatial_obj$`Tissue Category`=spe$Classifier_Label
          }
          
          pcf_dir = paste0(out_dir, "/pcf-inputs/")
          dir.create(pcf_dir, showWarnings = FALSE)
          
          out_name=glue('{pcf_dir}/{group}-vectra_cell_seg_data.csv')
          
          write.csv(spatial_obj,out_name,row.names = F)
        }
        
        incProgress(1/3, detail = 'Downloading new annotations to tempdir ...')
        print(as.character(tempdir))
        tempdir0 = as.character(tempdir)
        print(tempdir0)
	if('Seurat' %in% class(spe)){
		obj.sc = as.SingleCellExperiment(spe)
		#spe=SpatialExperiment(assays=assays(obj.sc),rowData=rowData(obj.sc),colData=colData(obj.sc),reducedDims=reducedDims(obj.sc))
		spe=SpatialExperiment(assays=assays(obj.sc),rowData=rowData(obj.sc),colData=colData(obj.sc),reducedDims=reducedDims(obj.sc),spatialCoords = as.matrix(data.frame('x'=obj.sc$x,'y'=obj.sc$y)))
		if(nrow(spe) > 100){
			m_ = seurat_markers()
			plot_heatmap.mod(x=spe,group_by = group,out_dir = tempdir0,features=m_$gene)
		}else{
			plot_heatmap.mod(x=spe,group_by = group,out_dir = tempdir0)
		}
                plot_dr.mod(spe,dr='UMAP',color_by = group,out_dir = tempdir0,h=10,w=20)
		plot_spatial.mod(spe,color_by = group,out_dir = tempdir0,h=10,w=20)
	}else{
        	plot_heatmap.mod(x=spe,group_by = group,out_dir = tempdir0)
        	plot_dr.mod(spe,dr='UMAP',color_by = group,out_dir = tempdir0,h=10,w=20)
        	plot_spatial.mod(spe,color_by = group,out_dir = tempdir0,h=10,w=20)
	}
	names(assays(spe))[grep('scale',names(assays(spe)))]='exprs'
        saveRDS(spe,glue('{tempdir0}/spe.rds'))
        spatial_coords = data.frame(spatialCoords(spe))
        spatial_coords$cluster = new_clusters_
        write.csv(spatial_coords,glue('{tempdir0}/annotated_spatial_coords-v{annotation.count+1}.csv'))
        
        prepare_pcf_inputs(spe = spe,out_dir = tempdir0,group=group)
        incProgress(1/3, detail = 'Done')
      })
    }
    
    
  })
  
  output$phenomenalist_download <- downloadHandler(
    filename = function(){
      paste("phenomenalist-modify_clusters-output-", Sys.Date(), ".zip", sep = "")
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
