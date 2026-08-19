RunPhenomenalist.shiny=function(segmentation_file,label,failed.markers=NULL,nuclear.markers=NULL,HALO=T,mask.only=NULL,out_dir=getwd(),clustering_res=seq(5,7),classifier_label=NULL,min.cells=1000,else_cytoplasm=F,max.cells=5e5){
  
  suppressPackageStartupMessages({
    library(phenomenalist)
    library(tidyverse)
    library(glue)
    library(cowplot)
    library(ggsci)
    library(shiny)
  })
  
  #source('/Users/ee699/TRIC/Phenomenalist-AnalysisEngine/RunPhenomenalist-Shiny/phenomenalist-utils-shiny.R')
  
  # Default cols to skip for HALO:
  basic_skip_cols = "Blank|blank|DAPI"
  
  data_loc=segmentation_file
  
  ## Default cols to skip for HALO:
  # here if using
  if(is.null(nuclear.markers)){
    halo_skip_cols = "Classification|Completeness|Cytoplasm|Nucleus|Membrane"
    if(length(failed.markers)>0){
      halo_skip_cols = "Classification|Completeness|Cytoplasm|Nucleus|Membrane|"
      halo_skip_cols=paste0(halo_skip_cols,paste(failed.markers,collapse = '|'),collapse = '|')
      
    }
  }else{
    halo_skip_cols = "Classification|Completeness|Membrane|"
    if(length(failed.markers)>0){
      halo_skip_cols=paste0(halo_skip_cols,paste(failed.markers,collapse = ' |'),collapse = '|')
      
    }
  }
  
  
  
  
  
  
  
  
  #labels.files=unlist(strsplit(data_loc,'[/]'))[length(unlist(strsplit(data_loc,'[/]')))]
  #label=unlist(strsplit(labels.files,'.csv'))[1]
  
  
  withProgress(message = 'Running Analysis', value = 0,{
    incProgress(1/6, detail = 'Determining expression columns')
    # by default NULL:
    message(failed.markers)
    message(nuclear.markers)
    expression.columns=phenomenalist.preprocess(segmentation_file,failed.markers = failed.markers,nuclear.markers = nuclear.markers,HALO = HALO,else.cytoplasm = else_cytoplasm)
    message(expression.columns)
    
    if(!HALO){
      skip_cols=paste0(c(basic_skip_cols,'Ch'),collapse = '|')
      if(length(failed.markers)>0){
        skip_cols=paste0(skip_cols,'|')
        skip_cols=paste0(skip_cols,paste(failed.markers,collapse = '|'),collapse = '|')
      }
      
      # if(length(data_loc)==1){
      #   data_csv=data_loc
      #   print(data_csv)
      #   x <- data.table::fread(data_csv, stringsAsFactors = FALSE, data.table = FALSE)
      #   expression.columns=colnames(x)[seq(6,dim(x)[2])]
      #   expression.columns=expression.columns[!str_detect(expression.columns,skip_cols)]
      # }
    }else{
      skip_cols=halo_skip_cols
    }
    
    
    # clustering resolutions:
    resolutions = clustering_res
    # which subset to store (by default all)
    idx=seq(length(resolutions))
    
    
    # Where output dirs should be housed:
    workDir=out_dir
    
    #if (!dir.exists(workDir)) { dir.create(workDir) }
    #setwd(workDir)
    print(workDir)
    
    out.label=NULL
    
    data_csv=segmentation_file
    .mem_log(glue("[pipeline] reading file: {data_csv} | R heap: {.mem_mb()} MB"))
    
    
    message(glue('Skipping columns: {skip_cols}'))
    
    out_dir = glue("{out_dir}/out-phenomenalist")
    if (!dir.exists(out_dir)) { dir.create(out_dir) }
    print('creating nested subdir')
    out_dir = paste0(out_dir, "/", label)
    dir.create(out_dir,showWarnings = F)
    message(out_dir)
    
    if('lock.rds' %in% list.files(out_dir)){
      .mem_log(glue("[pipeline] lock.rds found — loading existing spe.rds | R heap: {.mem_mb()} MB"))
      spe = readRDS(glue('{out_dir}/spe.rds'))
      .mem_log(glue("[pipeline] spe.rds loaded ({ncol(spe)} cells, {nrow(spe)} markers) | R heap: {.mem_mb()} MB"))
      spe = cluster.mod(spe, resolution = resolutions, out_dir = out_dir)
      names(colData(spe))
      saveRDS(spe,glue(out_dir,'/spe.rds'))
      
      cluster_cols = stringr::str_subset(names(colData(spe)), "cluster_leiden")
      
      for(i in seq(length(cluster_cols))){
        
        plot_heatmap.mod(spe, group_by = cluster_cols[i],auto = T,out_dir = out_dir,segment=T)
      }
      
      res_found=do.call('rbind',strsplit(cluster_cols,'[cluster_leiden_res]'))
      all_res=unique(as.numeric(res_found[,19]))
      for(res_val in all_res){
        prepare_mask_inputs(spe=spe,out_dir = out_dir,res = res_val,mask.only = mask.only,label=label)
      }
      
      
    }else{
      incProgress(1/6, detail = 'Generating Phenomenalist object')
      .mem_log(glue("[pipeline] create_object.mod starting | R heap: {.mem_mb()} MB"))
      spe <- create_object.mod(data_csv, skip_cols = skip_cols, transformation = "z", out_dir = out_dir,expression_cols = expression.columns,classifier.label = classifier_label,min.cells = min.cells,max.cells = max.cells)
      .mem_log(glue("[pipeline] create_object.mod done ({ncol(spe)} cells, {nrow(spe)} markers) | R heap: {.mem_mb()} MB"))

      if ("tile_num" %in% names(colData(spe))) plot_spatial(spe, color_by = "tile_num",out_dir = out_dir)
      
      if ("Classifier_Label" %in% names(colData(spe))) plot_spatial(spe, color_by = "Classifier_Label",out_dir = out_dir)
      
      incProgress(1/6, detail = 'Generating spatial objects')
      spatial_dir = paste0(out_dir, "/spatial-expression")
      dir.create(spatial_dir, showWarnings = FALSE)
      .mem_log(glue("[pipeline] plot_spatial starting ({length(names(spe))} markers) | R heap: {.mem_mb()} MB"))
      plot_spatial(spe, color_by = names(spe), out_dir = spatial_dir)
      .mem_log(glue("[pipeline] plot_spatial done | R heap: {.mem_mb()} MB"))
      plot_spatial(spe, color_by = names(spe), smooth = TRUE, out_dir = spatial_dir)
      .mem_log(glue("[pipeline] plot_spatial smooth done | R heap: {.mem_mb()} MB"))
      
      incProgress(1/6, detail = 'Generating dimensionality reductions')
      umap_dir = paste0(out_dir, "/UMAP-expression")
      dir.create(umap_dir, showWarnings = FALSE)
      .mem_log(glue("[pipeline] plot_dr UMAP starting | R heap: {.mem_mb()} MB"))
      plot_dr(spe, dr = "UMAP", color_by = names(spe), out_dir = umap_dir)
      .mem_log(glue("[pipeline] plot_dr UMAP done | R heap: {.mem_mb()} MB"))
      plot_dr(spe, dr = "UMAP", color_by = names(spe), smooth = TRUE, out_dir = umap_dir)
      .mem_log(glue("[pipeline] plot_dr UMAP smooth done | R heap: {.mem_mb()} MB"))
      
      
      
      incProgress(1/6, detail = 'Initiating clustering')
      .mem_log(glue("[pipeline] cluster.mod starting | R heap: {.mem_mb()} MB"))
      spe = cluster.mod(spe, resolution = resolutions, out_dir = out_dir)
      names(colData(spe))
      #saveRDS(spe,glue(out_dir,'/spe.rds'))
      write_rds(spe,glue(out_dir,'/spe.rds'))
      cluster_cols = stringr::str_subset(names(colData(spe)), "cluster_leiden")
      
      for(i in seq(length(cluster_cols))){
        
        plot_heatmap.mod(spe, group_by = cluster_cols[i],auto = T,out_dir = out_dir,segment=T)
      }
      
      res_found=do.call('rbind',strsplit(cluster_cols,'[cluster_leiden_res]'))
      all_res=unique(as.numeric(res_found[,19]))
      for(res_val in all_res){
        prepare_mask_inputs(spe=spe,out_dir = out_dir,res = res_val,mask.only = mask.only,label=label)
      }
      
    }
    incProgress(1/6, detail = 'Done')
    #print(getwd())
    #print(list.files())
    #tar(tarfile = "out-phenomenalist.tar.gz",compression = 'gzip')
    # return(glue('{out_dir}/out-phenomenalist.tar.gz'))
    
  })
}