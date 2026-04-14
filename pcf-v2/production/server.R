require(shiny)
require(shinyFiles)
require(stringr)
require(glue)
require(shinyjs)
require(cowplot)
require(reticulate)
require(reshape2)
require(ggpubr)
require(ggsignif)
require('shinycssloaders')
require(dplyr)
require(MatrixGenerics)
source('/srv/shiny-server/phenomenalist/utils/provenance.R')

server <- function(input, output, session) {

  
  mydata <- reactive({
    inFile <- input$vectras
    print(inFile)
    
    
    
    if (is.null(inFile)) {
      return(NULL)
    }
    
    
    
    # to account for multiple files:
    if (length(inFile) > 0) {
      
      
      # to account for multiple files:
      #tbls <- lapply(seq(length(inFile)),function(x) read.csv(glue("{inDir}/{inFile[x]}")))
      tbls <- lapply(seq(length(inFile$datapath)),function(x) read.csv(glue("{inFile$datapath[x]}")))
      print('loaded')
      
      tbls_phen <- lapply(tbls,function(x) unique(x[['Phenotype']]))
      set_ = unique(unlist(tbls_phen))
      print(set_)
      mapping_=do.call('cbind',lapply(seq(length(tbls_phen)),function(x){
        tmp=rep(0,length(set_))
        tmp[match(tbls_phen[[x]],set_)]=1
        return(tmp)
      }))
      return(set_[rowSums(mapping_)==length(tbls)])
      
    }
  })
  
  
  
  observeEvent(input$reset_button, {
    js$resetClick()
  })
  
  observe({
    celltypes <- mydata()
    updateCheckboxGroupInput(inputId = "celltypes_", choices = celltypes, selected = celltypes)
    
    updateSelectInput(inputId = 'ref_selection',choices=celltypes,selected="")
    
  })
  
  tempdir=file.path(tempdir(), as.integer(Sys.time()))
  dir.create(tempdir)
  print(as.character(tempdir))
  tempdir0 = as.character(tempdir)
  tracker <- ProvenanceTracker$new("pcf_v2", session, tempdir0)
  print(tempdir0)


  output$plots <- renderPlot({
    # action button:
    run <- input$Run
    print(run)
    if (run >0 ) {
      # text input:
      label <- input$run_label
      # instrument resolution:
      resolution <- input$instr_res
      
      
      library(reticulate)
      use_virtualenv('r-reticulate')
      message(py_config())
      #py_install("scipy")
      #py_install('matplotlib')
      #py_install("seaborn")
      #py_install("rpy2")
      source_python('/srv/shiny-server/phenomenalist/utils/spatial-shiny/vectra_lib_v4.py')      
      
      pd <- import("pandas")
      withProgress(message = "Running PCF Analysis", value = 0, {

	dp=input$vectras$datapath
        df=unlist(strsplit(dp[1],'[0-9]{1,10}.csv'))

	vectra_output <- extract_data(df)
        
        incProgress(1 / 5, detail = "Computing PCF")
        
        # create celltype selector UI
        celltypes <- input$celltypes_
       

        radius=input$radius
        #message(glue('vectra output: {vectra_output}'))
        pcf_df <- pcf(vectra_output, count_threshold = 10,cell_types=celltypes,radius=radius,resolution = resolution)
	message(reticulate::py_last_error())
	message(glue('pcf df'))
	message(dim(pcf_df)[1])
	message(dim(pcf_df)[2])

        incProgress(1 / 5, detail = "Generating PCF Curves")
        
        
        
        cell1 <- input$ref_selection
        
        # creating new sub-dir for output:
        tempdir1=glue('{tempdir0}/{label}')
        
        dir.create(tempdir1)
        
        ppc=plot_pcf_curves(pcf_df = pd$DataFrame(pcf_df), cell_types_ = celltypes, resolution = resolution, label = label, out_path = tempdir1,radius=radius)
        message(reticulate::py_last_error())

        sample.count=length(unique(ppc[[1]]$Patient))
	message(glue('sample count: {sample.count}'))
	message(head(ppc))         
        

        source("/srv/shiny-server/phenomenalist/utils/spatial-shiny/plot-pfcs-dev.R")
        
        plot_pcfs_R(ppc,out_dir = tempdir1,res=resolution,radius=radius)
        
	
	names(ppc)=c('All',celltypes)
        saveRDS(ppc,glue('{tempdir1}/ppc.rds'))
        
        
        
        incProgress(1 / 5, detail = "Generating PCF AUC Violins")
        message('generating violins')

        if(sample.count == 1){
          pcf_auc <- pcf_AUC(pcf_df = pd$DataFrame(pcf_df), cell_types_ = celltypes, cell1 = cell1, label = label, out_path = tempdir1, resolution = resolution)
          message(reticulate::py_last_error())

          incProgress(1 / 5, detail = "Exporting PCF data")
	  pcf_auc$Sample=rep(label,dim(pcf_auc)[1])

	  write.csv(pcf_auc, glue("{tempdir1}/{label}-PCF_AUCs.csv"))
          
          aucs.m=melt(pcf_auc)
          aucs.m$variable=gsub('[.]{2}','+',as.character(aucs.m$variable))
          aucs.m$variable=gsub('[.]{1}',' ',as.character(aucs.m$variable))
          
          
          p=ggviolin(aucs.m,x='variable',y='value',color = 'variable',add = 'boxplot')+geom_hline(yintercept = mean(aucs.m$value[aucs.m$variable=='All']))+geom_signif(map_signif_level = T,test.args=c('alternative'='two.sided'))+theme(legend.position = "none")+xlab('')+ylab('norm PCF')+ggtitle(glue('{cell1} Interactions'))+stat_compare_means(ref.group = 'All')+coord_flip()+stat_summary(fun = "mean",geom = "point",color = "red")
          ggsave(glue("{tempdir1}/{label}-PCF_AUC_violins.pdf"),p,w=16,h=10)
          print(p)
          
          incProgress(1 / 5, detail = "Done")
        }else{
          # use PPC indexed by cell1 to pull relevant series,
          # concatenate all samples into one long vector for each interaction
          # 
          cell2_list=c('All',celltypes[celltypes != cell1])
          cell1_tbl=ppc[[cell1]]
          
          print('cell1_tbl: ')
          print(cell1_tbl)
          # each entry will be a list corresponding to each sample:
          pcf_aucs = data.frame(do.call('cbind',lapply(cell2_list,function(x){
            cell2_data=cell1_tbl$normPCF[cell1_tbl$Cell_one==x & cell1_tbl$Cell_two==cell1]
            if(length(cell2_data) == 0){
              cell2_data=cell1_tbl$normPCF[cell1_tbl$Cell_one==cell1 & cell1_tbl$Cell_two==x]
            }
            return(unlist(cell2_data))
          })))
          print(cell1)
          print(cell2_list)
          
          print('PCF aucs:')
          print(head(pcf_aucs))
          names(pcf_aucs)=cell2_list
          
          aucs.m=melt(pcf_aucs)
          aucs.m$variable=gsub('[.]{2}','+',as.character(aucs.m$variable))
          aucs.m$variable=gsub('[.]{1}',' ',as.character(aucs.m$variable))
          
          
          incProgress(1 / 5, detail = "Exporting PCF data")
          p=ggviolin(aucs.m,x='variable',y='value',color = 'variable',add = 'boxplot')+geom_hline(yintercept = mean(aucs.m$value[aucs.m$variable=='All']))+geom_signif(map_signif_level = T,test.args=c('alternative'='two.sided'))+theme(legend.position = "none")+xlab('')+ylab('norm PCF')+ggtitle(glue('{cell1} Interactions'))+stat_compare_means(ref.group = 'All')+coord_flip()+stat_summary(fun = "mean",geom = "point",color = "red")
          ggsave(glue("{tempdir1}/{label}-PCF_AUC_violins-global.pdf"),p,w=16,h=10)
          print(p)
          
          samples=unique(cell1_tbl$Patient)
          samples_clean=do.call('rbind',strsplit(samples,'-vectra_cell_seg_data'))
          
          pcf_aucs$Sample=unlist(lapply(samples_clean,function(x) rep(x,length(cell1_tbl$normPCF[[1]]))))
          write.csv(pcf_aucs, glue("{tempdir1}/{label}-PCF_AUCs.csv"))
          
          aucs.m=melt(pcf_aucs,id='Sample')
          aucs.m$variable=gsub('[.]{2}','+',as.character(aucs.m$variable))
          aucs.m$variable=gsub('[.]{1}',' ',as.character(aucs.m$variable))
          
          
          p=ggviolin(aucs.m,x='Sample',y='value',color = 'variable',add = 'boxplot')+geom_hline(yintercept = mean(aucs.m$value[aucs.m$variable=='All']))+geom_signif(map_signif_level = T,test.args=c('alternative'='two.sided'))+xlab('')+ylab('norm PCF')+ggtitle(glue('{cell1} Interactions'))+coord_flip()+stat_summary(fun = "mean",geom = "point",color = "red")
          ggsave(glue("{tempdir1}/{label}-PCF_AUC_violins-bySample.pdf"),p,w=16,h=10)
          
          # individual tests:
          dir.create(glue('{tempdir1}/individual-samples/'))
          for(i in unique(pcf_aucs$Sample)){
            
            aucs.m.tmp=aucs.m %>% subset(Sample==i)
            
            p=ggviolin(aucs.m.tmp,x='variable',y='value',color = 'variable',add = 'boxplot')+geom_hline(yintercept = mean(aucs.m.tmp$value[aucs.m.tmp$variable=='All']))+geom_signif(map_signif_level = T,test.args=c('alternative'='two.sided'))+theme(legend.position = "none")+xlab('')+ylab('norm PCF')+ggtitle(glue('{cell1} Interactions'))+stat_compare_means(ref.group = 'All')+coord_flip()+stat_summary(fun = "mean",geom = "point",color = "red")
            ggsave(glue("{tempdir1}/individual-samples/{i}-PCF_AUC_violins.pdf"),p,w=16,h=10)
          }
          
          
          incProgress(1 / 5, detail = "Done")
        }
        
      })
    }
  })
  output$pcf_download <- downloadHandler(
    filename = function(){
      paste("pcf-output-", Sys.Date(), ".zip", sep = "")
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
