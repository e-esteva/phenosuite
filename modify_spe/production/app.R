require(shiny)
require(stringr)
require(glue)
require(shinyjs)
require(shinyFiles)
require(phenomenalist)
require(cowplot)
require(shinycssloaders)
require(dplyr)
require(ggrepel)

options(shiny.maxRequestSize=1000*1024^2) 
jsResetCode <- "shinyjs.resetClick = function() {history.go(0)}" # Define the js method that resets the page

ui <- fluidPage(
  tags$head(tags$style(".shiny-notification {position: fixed; 
                                             opacity: 1 ;
                       top: 80% ;
                       left: 40% ;
                       height: 90px;
                       width: 300px}")),
  
  # Add navigation tabs
  navbarPage("Clustering Analysis Tool",
    
    # Sub-Clustering Panel
    tabPanel("Sub-Clustering",
      sidebarLayout( 
        sidebarPanel(
          fileInput("file1", "Choose RDS File",
                    accept = c('.rds')
          ),
          checkboxGroupInput('rb0','Select Clustering',choices = "",inline = T,selected = NULL),
          conditionalPanel(condition="!is.null(input.rb0)",checkboxGroupInput('target_cluster','Target Cluster(s)',choices='',inline=T,selected = NULL)),
          
          tags$hr(),
          
          sliderInput("cluster_res_range", "Clustering Resolution Range:",
                      min = 1, max = 10,
                      value = c(1,3)),
        ),
        mainPanel(
          plotOutput("plots_clusters") %>% withSpinner(color="#0dc5c1"),
          tags$hr(),
          textInput('run_label','Results Name'),
          selectInput('view','full or binary view of sub-clusters',choices=c('full','binary'),selected='full'),
          tags$hr(),
          conditionalPanel(condition="!is.null(input.run_label)",actionButton('run','Initiate Sub-Clustering')),
          tags$hr(),
          conditionalPanel(condition='input.run>0',downloadButton(
            outputId = "sub_clustering_download",
            label = "Download Results",
            icon = icon("file-download")
          )),
          tags$hr(),
          useShinyjs(), 
          extendShinyjs(text = jsResetCode, functions = "resetClick"),
          actionButton("reset_button", "Reset Page")
        )
      )
    ),
    
    # Re-Clustering Panel
    tabPanel("Re-Clustering",
      sidebarLayout(
        sidebarPanel(
          fileInput("file2", "Choose RDS File",
                    accept = c('.rds')
          ),
          
          tags$hr(),
          verbatimTextOutput("existing_clusters"),
	  tags$hr(),
          sliderInput("recluster_res_range", "Re-clustering Resolution Range:",
                      min = 1, max = 10, step = 1,
                      value = c(1.0, 3.0)),
        ),
        mainPanel(
          plotOutput("recluster_plots") %>% withSpinner(color="#0dc5c1"),
          tags$hr(),
          textInput('recluster_label','Results Name'),
          tags$hr(),
          conditionalPanel(condition="!is.null(input.recluster_label)",actionButton('run_recluster','Initiate Re-Clustering')),
          tags$hr(),
          conditionalPanel(condition='input.run_recluster>0',downloadButton(
            outputId = "recluster_download",
            label = "Download Results",
            icon = icon("file-download")
          )),
          tags$hr(),
          useShinyjs(), 
          extendShinyjs(text = jsResetCode, functions = "resetClick"),
          actionButton("reset_button2", "Reset Page")
        )
      )
    )
  )
)

server <- function(input, output, session) {
  
  # ===== SUB-CLUSTERING FUNCTIONS =====
  tempdir0=glue('/apps/home/rtmp/{session$token}')
  dir.create(tempdir0)

  message(glue('tempdir0: {tempdir0}'))

  # return spe object for sub-clustering:
  mydata0 <- reactive({
    inFile <- input$file1
    
    if (is.null(inFile)){
      return(NULL)
    }else{
      spe = readRDS(inFile$datapath)
      return(spe)
    }
  })
  
  # return clustering columns for sub-clustering
  mydata1 <- reactive({
    spe= mydata0()
    if (is.null(spe)){
      return(NULL)
    }else{
      return(str_subset(names(colData(spe)),'cluster'))
    }
  })
  
  # returns clusters in selected group for sub-clustering:
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
  observeEvent(input$reset_button2, {js$resetClick()})
  
  #tempdir=file.path(tempdir(), as.integer(Sys.time()))
  #dir.create(tempdir)
  #tempdir0 = as.character(tempdir)
  
  sub_clustered_obj=reactive({
    source('/srv/shiny-server/phenomenalist/utils/RunPhenomenalist-shiny/phenomenalist-utils-shiny.R')
    
    spe = mydata0()
   
    if(!is.null(spe)){
      x_ = names(colData(spe))
      group = input$rb0
      message(group)
      target_cluster=input$target_cluster
      message('target cluster: ')
      message(target_cluster)
      message("=== DEBUG: Input target_cluster ===")
      message(paste(target_cluster, collapse = ", "))
      message("=== DEBUG: Length of target_cluster ===")
      message(length(target_cluster))
      message("=== DEBUG: Class of target_cluster ===")
      message(class(target_cluster))
      sub_clustering_ranges=input$cluster_res_range
      message('desired clustering resolution(s):')
      message(sub_clustering_ranges)
    }else{
      return(NULL)
    }
    
    if(input$run > 0){
      withProgress(message = glue("Sub-Clustering {group}"), value = 0, {
        # Function to create safe cluster identifiers
	create_cluster_identifier <- function(cluster_names) {
  		# Create a short hash of all cluster names for uniqueness
  		cluster_hash <- digest::digest(paste(cluster_names, collapse = "|"), 
                                algo = "md5", serialize = FALSE)
  		short_hash <- substr(cluster_hash, 1, 8)
  
  		# Count of clusters
  		cluster_count <- length(cluster_names)
  
  		return(paste0(cluster_count, "clusters_", short_hash))
	}
	# Initialize variables to track results
  	all_new_columns <- c()
  	last_new_col_name <- NULL
  	last_new_col_name_bin <- NULL

        for(r in seq(min(sub_clustering_ranges),max(sub_clustering_ranges))){
          
          new_cols = as.character(spe[[group]])
          # CRITICAL DEBUG: Check target_cluster again here
          message("=== DEBUG: target_cluster at start of loop ===")
          message(paste(target_cluster, collapse = ", "))
	  processed_clusters <- c()

    	  for(cluster in target_cluster){
            message("=== DEBUG: Processing cluster ===")
            message(cluster)

            incProgress(1 /(length(target_cluster)*length(sub_clustering_ranges)) , detail = glue("sub-clustering cluster: {cluster} @ leiden resolution: {r}"))
            
            spe.s = spe[,spe[[group]]==cluster]
            spe.s=cluster.mod(spe.s,resolution = r )
            print(head(colData(spe.s)))
            message(table(spe.s[[group]]))
            
            y=names(colData(spe.s))
            message('new col names:')
            message(y)
            message(all.equal(x_,y))
            if(all.equal(x_,y)!= T){
              new_group=y[!(union(x_,y) == intersect(x_,y))]
              message(glue('new_group: {new_group}'))
              message(table(spe.s[[new_group]]))
              
              new_cols[new_cols == cluster]= paste(glue('{cluster}_{as.character(spe.s[[new_group]])}'),sep="_")
              message(table(new_cols))
	      processed_clusters <- c(processed_clusters, cluster)

            }else{
              target_cluster=target_cluster[-match(cluster,target_cluster)]
            }
          }
	  # CRITICAL DEBUG: Check what target_cluster is when creating filename
          message("=== DEBUG: target_cluster when creating filename ===")
          message(paste(target_cluster, collapse = ", "))
	  # The problem is likely HERE - let's see what's actually in target_cluster
          message("=== DEBUG: About to create filename with these clusters ===")
          for(i in seq_along(target_cluster)) {
          	message(glue("Cluster {i}: '{target_cluster[i]}'"))
          }
	  if(length(processed_clusters) > 0) {
    		# FIXED: Create short, safe column names
    		cluster_id <- create_cluster_identifier(processed_clusters)
    		new_col_name <- glue('{group}_subclustered_res{r}_{cluster_id}')
    
    		message("=== DEBUG: Generated SHORT column name ===")
    		message(new_col_name)
    
    		spe[[new_col_name]] <- new_cols
    
    		# Create binary version
    		new_col_name_bin <- paste0(new_col_name, "_binary")
    		old_col_names <- as.character(spe[[group]])
    		new_col_bin <- as.vector(sapply(as.character(spe[[new_col_name]]), 
                                   function(x) ifelse(x %in% old_col_names, 'else', x)))
    		spe[[new_col_name_bin]] <- new_col_bin
    
    		message("=== DEBUG: Short filenames will be ===")
    		message(paste("Binary:", new_col_name_bin))
    		message(paste("Regular:", new_col_name))
    		# Track the columns we've created
      		all_new_columns <- c(all_new_columns, new_col_name, new_col_name_bin)
      		last_new_col_name <- new_col_name
      		last_new_col_name_bin <- new_col_name_bin

    		# Save metadata about which clusters were actually processed
    		cluster_metadata <- data.frame(
      			column_name = new_col_name,
      			resolution = r,
      			selected_cluster_count = length(processed_clusters),
      			processed_clusters = paste(processed_clusters, collapse = " | "),
      			timestamp = Sys.time(),
      			stringsAsFactors = FALSE
    		)
    
    		# Write metadata file
    		metadata_file <- file.path(tempdir0, paste0("metadata_", new_col_name, ".csv"))
    		write.csv(cluster_metadata, metadata_file, row.names = FALSE)
    
    		# Now the plot functions will use short names
    		tryCatch({
      			plot_heatmap.mod(spe, group_by = new_col_name_bin, out_dir = tempdir0)
      			plot_heatmap.mod(spe, group_by = new_col_name, out_dir = tempdir0)
    		}, error = function(e) {
      		message("=== DEBUG: Error in plot_heatmap.mod ===")
      		message(e$message)
    		})
    
    		x_ = names(colData(spe))
    
  	} else {
    		message("No clusters were successfully processed for sub-clustering")
    		new_col_name <- NULL
  	}
          #new_col_name=glue('{group}-sub_clustered-res{r}-{paste0(target_cluster,collapse=",")}')
          #message("=== DEBUG: Generated column name ===")
          #message(new_col_name)
	  #spe[[new_col_name]]=new_cols
          
          #new_col_name_bin = glue('{new_col_name}-binary')
          #old_col_names=as.character(spe[[group]])
          #new_col_bin = as.vector(sapply(as.character(spe[[new_col_name]] ),function(x) ifelse(x %in% old_col_names,'else',x)))
          #spe[[new_col_name_bin]]=new_col_bin
	  
	  #message("=== DEBUG: About to call plot_heatmap.mod ===")
          #message("Binary column name:", new_col_name_bin)
          #message("Regular column name:", new_col_name)
	  # The filename error happens in these calls:
          #tryCatch({
          #	plot_heatmap.mod(spe, group_by = new_col_name_bin, out_dir = tempdir0)
          #	plot_heatmap.mod(spe, group_by = new_col_name, out_dir = tempdir0)
          #}, error = function(e) {
          #message("=== DEBUG: Error in plot_heatmap.mod ===")
          #message(e$message)
          #})
          #plot_heatmap.mod(spe,group_by = new_col_name_bin,out_dir = tempdir0)
          #plot_heatmap.mod(spe,group_by = new_col_name,out_dir = tempdir0)
          
          x_ = names(colData(spe))
        }
      })
      
      saveRDS(spe,glue('{tempdir0}/spe.rds'))
      
      dir.create(glue('{tempdir0}/mask-inputs/'))
      sub_clusters=str_subset(names(colData(spe)),pattern = 'subclustered')
      spatial_df=data.frame(spatialCoords(spe))
      
      for(i in sub_clusters){
        spatial_df$cluster=spe[[i]]
        write.csv(spatial_df,glue('{tempdir0}/mask-inputs/{i}-spatial_anno.csv'))
      }
      
      return(list(spe,last_new_col_name,last_new_col_name_bin))
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
      new_col_name_bin=spe.sub_clustered[[3]]
      if(input$view == 'binary'){
        #new_col_name_bin = glue('{new_col_name}-binary')
        #old_col_names=as.character(spe[[group]])
        #new_col_bin = as.vector(sapply(as.character(spe[[new_col_name]]),function(x) ifelse(x %in% old_col_names,'else',x)))
        #spe[[new_col_name_bin]]=new_col_bin
	#message('new_col_name_bin')
	#message(new_col_name_bin)
	#message('new_col_name')
	#message(new_col_name)
	message(names(colData(spe)))
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
  
  # ===== RE-CLUSTERING FUNCTIONS =====
  
  # return spe object for re-clustering:
  recluster_data <- reactive({
    inFile <- input$file2
    
    if (is.null(inFile)){
      return(NULL)
    }else{
      spe = readRDS(inFile$datapath)
      return(spe)
    }
  })
  
  
  # Display existing clustering resolutions
  output$existing_clusters <- renderText({
    spe = recluster_data()
    if(!is.null(spe)){
      existing_clusters = str_subset(names(colData(spe)), 'cluster')
      if(length(existing_clusters) > 0){
        return(paste(existing_clusters, collapse = "\n"))
      } else {
        return("No existing clustering columns found")
      }
    }
    return("")
  })
  
  # Control visibility of existing clusters section
  output$show_existing_clusters <- reactive({
    spe = recluster_data()
    if(!is.null(spe)){
      existing_clusters = str_subset(names(colData(spe)), 'cluster')
      return(length(existing_clusters) > 0)
    }
    return(FALSE)
  })
  
  # Make show_existing_clusters available for conditional panel
  outputOptions(output, "show_existing_clusters", suspendWhenHidden = FALSE)
  
  # Create temp directory for re-clustering
  recluster_tempdir = file.path(tempdir(), paste0("recluster_", session$token, "_", as.integer(Sys.time())))
  dir.create(recluster_tempdir, showWarnings = FALSE)
  recluster_tempdir0 = as.character(recluster_tempdir)
  
  reclustered_obj <- reactive({
    source('/srv/shiny-server/phenomenalist/utils/RunPhenomenalist-shiny/phenomenalist-utils-shiny.R')
    
    spe = recluster_data()
    
    if(!is.null(spe) && input$run_recluster > 0){
      
      # Store existing clustering columns before processing
      existing_cluster_cols = str_subset(names(colData(spe)), 'cluster')
      print('Existing clustering columns:')
      print(existing_cluster_cols)
      
      recluster_ranges = input$recluster_res_range
      print('Re-clustering resolution range:')
      print(recluster_ranges)
      
      withProgress(message = "Re-Clustering", value = 0, {
        
        # Iterate through the resolution range
        resolutions_to_process = seq(min(recluster_ranges), max(recluster_ranges), by = 1)
        
        for(i in seq_along(resolutions_to_process)){
          resolution = resolutions_to_process[i]
          
          incProgress(1/length(resolutions_to_process), 
                     detail = glue("Re-clustering @ leiden resolution: {resolution}"))
          
          # Re-cluster at the selected resolution
          spe = cluster.mod(spe, resolution = resolution)
        }
      })
      
      # Identify NEW clustering columns (those added during this process)
      all_cluster_cols = str_subset(names(colData(spe)), 'cluster')
      new_cluster_cols = setdiff(all_cluster_cols, existing_cluster_cols)
      
      print('New clustering columns:')
      print(new_cluster_cols)
      
      # Generate heatmaps ONLY for new clustering columns
      if(length(new_cluster_cols) > 0){
        for(col in new_cluster_cols){
          plot_heatmap.mod(spe, group_by = col, out_dir = recluster_tempdir0)
        }
      }
      
      # Save the updated spe object
      saveRDS(spe, glue('{recluster_tempdir0}/spe_reclustered.rds'))
      
      # Create spatial annotation files ONLY for new clustering columns
      if(length(new_cluster_cols) > 0){
        dir.create(glue('{recluster_tempdir0}/mask-inputs/'), showWarnings = FALSE)
        spatial_df = data.frame(spatialCoords(spe))
        
        for(i in new_cluster_cols){
          spatial_df$cluster = spe[[i]]
          write.csv(spatial_df, glue('{recluster_tempdir0}/mask-inputs/{i}-spatial_anno.csv'))
        }
      }
      
      # Return the spe object and the lowest NEW clustering column
      lowest_new_col = if(length(new_cluster_cols) > 0) new_cluster_cols[1] else NULL
      
      return(list(spe, lowest_new_col, new_cluster_cols))
    }else{
      return(NULL)
    }
  })
  
  output$recluster_plots <- renderPlot({
    spe = recluster_data()
    spe_reclustered = reclustered_obj()
    
    if(!is.null(spe_reclustered)){
      spe = spe_reclustered[[1]]
      lowest_new_clustering_col = spe_reclustered[[2]]
      
      # Generate heatmap with the lowest NEW clustering if available
      if(!is.null(lowest_new_clustering_col)){
        plot_heatmap.mod(spe, group_by = lowest_new_clustering_col, out_dir = NULL)
      } else {
        # If no new clustering was created, show a message
        plot(1, type="n", xlab="", ylab="", main="No new clustering generated")
        text(1, 1, "Try different resolution parameters")
      }
      
    } else if(!is.null(spe)) {
      # Show existing clustering if available
      cluster_cols = str_subset(names(colData(spe)), 'cluster')
      if(length(cluster_cols) > 0){
        plot_heatmap.mod(spe, group_by = cluster_cols[1], out_dir = NULL)
      } else {
        # Create a simple plot if no clustering is available
        plot(1, type="n", xlab="", ylab="", main="No clustering data available")
        text(1, 1, "Upload an RDS file with clustering data")
      }
    } else {
      return(NULL)
    }
  })
  
  output$recluster_download <- downloadHandler(
    filename = function(){
      glue("{input$recluster_label}-reclustering-output-{Sys.Date()}.zip")
    },
    content = function(file){
      print(recluster_tempdir0) 
      zip::zip(
        zipfile = file,
        files = dir(recluster_tempdir0),
        root = recluster_tempdir0
      )
    },
    contentType='application/zip'
  )
}

# Run the application 
shinyApp(ui = ui, server = server)
