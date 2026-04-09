require(dplyr)
require(ggrepel)


server <- function(input, output,session) {
  # return spe object:
  mydata0 <- reactive({
    
    inFile <- input$file1
    
    if (is.null(inFile)){
      return(NULL)
    }else{
      spe = readRDS(inFile$datapath)
      
      return(spe)
    }
     
    
    
    
  })
  # return clustering columns
  mydata1 <- reactive({
    
    #inFile <- input$file1
    spe= mydata0()
    if (is.null(spe)){
      return(NULL)
    }else{
      #spe = readRDS(inFile$datapath)
      
      return(str_subset(names(colData(spe)),'cluster'))
    }
      
    
    
    
  })
  
  
  # returns clusters in selected group:
  mydata2 <- reactive({
    
    inFile <- input$file1
    
    if (is.null(inFile) || is.null(input$rb0)){
      return(NULL)
    }else{
      spe = readRDS(inFile$datapath)
      
      group=input$rb0
      
      print(glue('group: {group}'))
      return(unique(spe[[group]]))
      
      
    }
    
    
    
    
  })
  
  
  observe({
    clustering_resolutions = mydata1()
    updateCheckboxGroupInput(inputId = 'rb0',choices = clustering_resolutions,selected = "")
    
    
    
    
  })
  
  
  observe({
    clusters = mydata2()
    updateCheckboxGroupInput(inputId = 'target_cluster',choices = clusters,selected = "")
    
    
  })
  
  observeEvent(input$reset_button, {js$resetClick()})
  
  tempdir=file.path(tempdir(), as.integer(Sys.time()))
  dir.create(tempdir)
  print(as.character(tempdir))
  tempdir0 = as.character(tempdir)
  print(tempdir0)
  
  sub_clustered_obj=reactive({
    source('/srv/shiny-server/Phenoptics-Menu/utils/RunPhenomenalist-shiny/phenomenalist-utils-shiny.R')
    
    
    
    spe = mydata0()
   

    if(!is.null(spe)){
      x_ = names(colData(spe))
      
      group = input$rb0
      print(group)
      
      target_cluster=input$target_cluster
      
      print('target cluster: ')
      print(target_cluster)
      
      sub_clustering_ranges=input$cluster_res_range
      print('desired clustering reslution(s):')
      print(sub_clustering_ranges)
    }else{
      return(NULL)
    }
    
    if(input$run > 0){
      withProgress(message = glue("Sub-Clustering {group}"), value = 0, {
        
        
      for(r in seq(min(sub_clustering_ranges),max(sub_clustering_ranges))){
        
        new_cols = as.character(spe[[group]])
        for(cluster in target_cluster){
          
          incProgress(1 /(length(target_cluster)*length(sub_clustering_ranges)) , detail = glue("sub-clustering cluster: {cluster} @ leiden resolution: {r}"))
          
          spe.s = spe[,spe[[group]]==cluster]
          
          spe.s=cluster.mod(spe.s,resolution = r )
          print(head(colData(spe.s)))
          print(table(spe.s[[group]]))
          
          y=names(colData(spe.s))
          
          print('new col names:')
          print(y)
          print(all.equal(x_,y))
          if(all.equal(x_,y) != T){
            new_group=y[!(union(x_,y) == intersect(x_,y))]
            print(glue('new_group: {new_group}'))
            
            print(table(spe.s[[new_group]]))
            
            
            new_cols[new_cols == cluster]= paste(glue('{cluster}_{as.character(spe.s[[new_group]])}'),sep="_")
            print(table(new_cols))
            
          }else{
            target_cluster=target_cluster[-match(cluster,target_cluster)]
          }
          
          
        }
        new_col_name=glue('{group}-sub_clustered-res{r}-{paste0(target_cluster,collapse=",")}')
        spe[[new_col_name]]=new_cols
        
	new_col_name_bin = glue('{new_col_name}-binary')
	old_col_names=as.character(spe[[group]])
	new_col_bin = as.vector(sapply(as.character(spe[[new_col_name]] ),function(x) ifelse(x %in% old_col_names,'else',x)))
        spe[[new_col_name_bin]]=new_col_bin

	plot_heatmap.mod(spe,group_by = new_col_name_bin,out_dir = tempdir0)
	plot_heatmap.mod(spe,group_by = new_col_name,out_dir = tempdir0)

        
        
        x_ = names(colData(spe))
      }
      
      
      })
      
      saveRDS(spe,glue('{tempdir0}/spe.rds'))
      
      dir.create(glue('{tempdir0}/mask-inputs/'))
      sub_clusters=str_subset(names(colData(spe)),pattern = 'sub_clustered')
      spatial_df=data.frame(spatialCoords(spe))
      
      for(i in sub_clusters){
        spatial_df$cluster=spe[[i]]
        write.csv(spatial_df,glue('{tempdir0}/mask-inputs/{i}-spatial_anno.csv'))
      }
      
      return(list(spe,new_col_name))
      
    }else{
      return(NULL)
    }
    
  })
  output$plots_clusters=renderPlot({
    group = input$rb0
    spe = mydata0()
    spe.sub_clustered=sub_clustered_obj()
    if(!is.null(spe.sub_clustered)){
      spe=spe.sub_clustered[[1]]
      new_col_name=spe.sub_clustered[[2]]

      if(input$view == 'binary'){
                new_col_name_bin = glue('{new_col_name}-binary')
                old_col_names=as.character(spe[[group]])
                new_col_bin = as.vector(sapply(as.character(spe[[new_col_name]]),function(x) ifelse(x %in% old_col_names,'else',x)))
                spe[[new_col_name_bin]]=new_col_bin

                plot_heatmap.mod(spe,group_by = new_col_name_bin,out_dir = NULL)
      }else{
     	 plot_heatmap.mod(spe,group_by = new_col_name,out_dir = NULL)
      }
    }else{
      if(!is.null(spe)){
	      group = input$rb0
     	      entropy = sapply(unique(spe[[group]]),function(x){
      			spe.s = spe[,spe[[group]]==x]
			values_=as.numeric(as.character(rowMeans(assay(x = spe.s,i = 'exprs'))))
			# to account for negative values:
			values_= exp(values_)
			shannons.entropy = -1*sum(log2(values_/sum(values_))*(values_/sum(values_)))
			return(shannons.entropy)
			

      
      	      })
               
	     entropy_df = data.frame('cluster'= as.character(unique(spe[[group]])),'entropy'=as.numeric(as.character(entropy)))
	     entropy_df$z_score = (entropy_df$entropy - mean(entropy_df$entropy))/sd(entropy_df$entropy)
	     
	     entropy_df = entropy_df %>% arrange(desc(z_score))
	     target_clusters=input$target_cluster
	     if(!is.null(target_clusters)){
		p=ggplot(entropy_df,aes(x=cluster,y=z_score,label=cluster,col=z_score)) +geom_point()+coord_flip()+geom_label_repel(data=entropy_df[match(target_clusters,entropy_df$cluster),],aes(label=cluster))+theme_bw() +ylab("Shannon's Entropy Z-Score")
	     }else{

	     	p=ggplot(entropy_df,aes(x=cluster,y=z_score,label=cluster,col=z_score)) +geom_point()+coord_flip()+geom_text(data=subset(entropy_df,z_score > 2),aes(label=cluster))+theme_bw() +ylab("Shannon's Entropy Z-Score")
	     }
	     ggsave(glue('{tempdir0}/Entropy-by-cluster-selected_clusters-laeled.pdf'),p,h=10,w=10)
	     return(p)
      
      }else{


      
      return(NULL)
     }
    }
  })
  output$sub_clustering_download <- downloadHandler(
      filename = function(){
        glue("{input$run_label}-sub_clustering-output-{Sys.Date()}.zip")
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



