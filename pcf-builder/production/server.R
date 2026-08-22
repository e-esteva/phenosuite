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

# Mean-interaction heatmap helpers. Kept in their own file so the plotting
# logic can be tested without Shiny; deploy pcf-heatmap.R next to server.R.
if (file.exists("pcf-heatmap.R")) {
  source("pcf-heatmap.R")
} else {
  source('/srv/shiny-server/phenomenalist/utils/pcf-heatmap.R')
}

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

    group_names = unique(samples)
    updateSelectizeInput(inputId = 'sample_order', choices = group_names, selected = group_names)

  })
  
  tempdir=file.path(tempdir(), as.integer(Sys.time()))
  dir.create(tempdir)
  print(as.character(tempdir))
  tempdir0 = as.character(tempdir)
  tracker <- ProvenanceTracker$new("pcf_builder", session, tempdir0)
  print(tempdir0)

  # Applies the user's sample -> group renaming to the pooled table. Shared by
  # the violin panel and the heatmap so both label groups identically.
  renamed_pcf <- reactive({
    req(mydata0())
    samples = mydata1()
    groups  = new_group_names()
    global_pcf = mydata0()

    for(i in seq(length(groups))){
      group.tmp=groups[i]
      global_pcf$Sample=gsub(samples[i],group.tmp,as.character(global_pcf$Sample))
    }
    global_pcf
  })

  # Builds the current violin plot from live inputs — no side effects (no
  # saving, no provenance tracking). Shared by the live preview and the
  # confirm-triggered save, so the saved plot always matches what's shown.
  build_pcf_plot <- reactive({
    req(mydata0(), input$celltype_to_analyze, input$ref_selection)

    global_pcf = renamed_pcf()
    ref = input$ref_selection
    celltype=input$celltype_to_analyze

    # User-controlled plot order, top-to-bottom. ggviolin()+coord_flip()
    # renders factor level 1 at the bottom of the flipped axis, so the
    # levels must be set in *reverse* of the user's chosen top-to-bottom
    # order for the on-screen order to match what they picked.
    present_groups = unique(global_pcf$Sample)
    order_ = input$sample_order
    if(!is.null(order_) && setequal(order_, present_groups)){
      global_pcf$Sample = factor(global_pcf$Sample, levels = rev(order_))
    }else{
      global_pcf$Sample = factor(global_pcf$Sample, levels = rev(sort(present_groups)))
    }

    # Reserve headroom on the (flipped) value axis and disable clipping —
    # stat_compare_means() places its p-value labels past the data range,
    # which otherwise get cut off at the panel edge after coord_flip().
    ggviolin(global_pcf,x='Sample',y=celltype,color = 'Sample',add = 'boxplot')+geom_hline(yintercept = mean(global_pcf[,match(celltype,names(global_pcf))][global_pcf$Sample==ref]))+theme(legend.position = "none")+xlab('')+ylab('norm PCF')+ggtitle(glue('{celltype} Interactions | all versus {ref}'))+stat_compare_means(ref.group = ref)+scale_y_continuous(expand = expansion(mult = c(0.05, 0.3)))+coord_flip(clip = "off")+stat_summary(fun = "mean",geom = "point",color = "red")
  })

  # ---- Mean interaction heatmap ------------------------------------------
  # One row per available cell type, one column per labelled sample group,
  # each cell the mean normalised PCF across that group's observations.
  # Column order follows the same "Plot order" control as the violins; here
  # it is used left-to-right, with no coord_flip to reverse.
  heatmap_matrix <- reactive({
    global_pcf = renamed_pcf()
    celltypes  = mydata()
    req(celltypes)

    present_groups = unique(as.character(global_pcf$Sample))
    order_ = input$sample_order
    group_order = if(!is.null(order_) && setequal(order_, present_groups)) order_
                  else sort(present_groups)

    m = pcf_heatmap_matrix(global_pcf, celltypes, group_order = group_order)
    m[pcf_heatmap_row_order(m), , drop = FALSE]
  })

  build_heatmap <- reactive({
    m = heatmap_matrix()
    req(nrow(m) > 0, ncol(m) > 0)
    lbl = if(!is.null(input$run_label) && nzchar(input$run_label)) input$run_label else NULL
    pcf_heatmap_plot(
      m,
      title       = if(is.null(lbl)) "Mean PCF interactions"
                    else glue("{lbl} — mean PCF interactions"),
      subtitle    = glue("{ncol(m)} sample group(s)  |  {nrow(m)} cell types"),
      show_values = isTRUE(input$hm_show_values),
      cap_quantile = if(is.null(input$hm_cap)) 0.98 else input$hm_cap
    )
  })

  # Live preview: cheap to recompute (means over the pooled table), so it
  # follows the rename/order controls immediately rather than waiting on a
  # button the way the violin panel does.
  output$heatmap = renderPlot({ build_heatmap() })

  # Explicit export, once per click: writes the figure and the underlying
  # matrix into the same tempdir the violins use, so "Download Results" picks
  # both up in the zip.
  observeEvent(input$export_heatmap, {
    req(mydata0())
    tracker$capture_parameters(input)
    tracker$analysis_started()

    m   = heatmap_matrix()
    p   = build_heatmap()
    sz  = pcf_heatmap_size(m)
    lbl = if(!is.null(input$run_label) && nzchar(input$run_label)) input$run_label else "pcf"
    safe = gsub("[^A-Za-z0-9]+", "_", lbl)

    ggsave(glue('{tempdir0}/{safe}-mean_interaction_heatmap.pdf'), p,
           width = sz$width, height = sz$height, limitsize = FALSE)
    ggsave(glue('{tempdir0}/{safe}-mean_interaction_heatmap.png'), p,
           width = sz$width, height = sz$height, dpi = 200, limitsize = FALSE)
    # The matrix behind the figure, so values can be reused without
    # re-deriving them from the AUC tables.
    write.csv(as.data.frame(m), glue('{tempdir0}/{safe}-mean_interaction_matrix.csv'))

    showNotification(glue("Heatmap exported ({nrow(m)} x {ncol(m)})"), type = "message")
    tracker$analysis_completed()
  }, ignoreInit = TRUE)

  # Live preview — recomputes whenever celltype/reference/order/groups
  # change, but stays blank until Confirm has been clicked at least once
  # (matches prior behaviour). Does not save anything by itself.
  output$plot=renderPlot({
    req(input$confirm_pcf > 0)
    build_pcf_plot()
  })

  # Explicit save, once per Confirm click — lets the user step through
  # multiple celltypes (select celltype, Confirm, select another, Confirm,
  # ...) and accumulate one PDF per celltype rather than overwriting a
  # single file, so "Download Results" later zips up all of them.
  observeEvent(input$confirm_pcf, {
    req(mydata0())
    tracker$capture_parameters(input)
    tracker$analysis_started()

    p        <- build_pcf_plot()
    celltype <- input$celltype_to_analyze
    safe_ct  <- gsub("[^A-Za-z0-9]+", "_", celltype)
    # Wider than the ggsave default (7x7) — the previous size left the
    # stat_compare_means() p-value labels truncated at the right edge.
    ggsave(glue('{tempdir0}/{input$run_label}-{safe_ct}.pdf'), p, width = 10, height = 7)

    tracker$analysis_completed()
  }, ignoreInit = TRUE)
  
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
