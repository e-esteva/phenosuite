RunPhenomenalist.shiny=function(segmentation_file,label,failed.markers=NULL,nuclear.markers=NULL,mask.only=NULL,out_dir=getwd(),clustering_res=seq(5,7),classifier_label=NULL,min.cells=1000,else_cytoplasm=F,max.cells=5e5){
  
  suppressPackageStartupMessages({
    library(phenomenalist)
    library(tidyverse)
    library(glue)
    library(cowplot)
    library(ggsci)
    library(shiny)
  })
  
  # Default cols to always skip (non-expression utility columns)
  basic_skip_cols = "Blank|blank"
  # DAPI is typically a nuclear stain used for segmentation only,
  # but in some panels it IS an expression marker. We'll add it to skip
  # only if it's not in the selected expression columns (checked below).
  
  data_loc=segmentation_file
  
  # Build skip_cols generically based on what preprocessing selects
  # These columns should be excluded from the object but are not expression data
  skip_cols_parts <- c("Classification", "Completeness", "Membrane")
  
  # If no nuclear markers specified, also skip Nucleus/Cytoplasm compartment cols
  # (they won't be selected as expression columns anyway)
  if(is.null(nuclear.markers)){
    skip_cols_parts <- c(skip_cols_parts, "Cytoplasm", "Nucleus")
  }
  
  # Add failed markers to skip list
  if(length(failed.markers) > 0){
    skip_cols_parts <- c(skip_cols_parts, failed.markers)
  }
  
  skip_cols <- paste0(c(basic_skip_cols, paste(skip_cols_parts, collapse = "|")), collapse = "|")
  
  withProgress(message = 'Running Analysis', value = 0,{
    incProgress(1/6, detail = 'Determining expression columns')
    # Use the generalized select_intensity_columns (matches intensity OR mean)
    message(failed.markers)
    message(nuclear.markers)
    expression.columns=select_intensity_columns(
      filepath = segmentation_file,
      failed_markers = failed.markers,
      nuclear_markers = nuclear.markers,
      prefer_cytoplasm = else_cytoplasm
    )
    message(expression.columns)
    
    # Add DAPI to skip_cols only if no selected expression column contains DAPI
    if(!any(grepl("DAPI", expression.columns, ignore.case = TRUE))) {
      skip_cols <- paste0(skip_cols, "|DAPI")
      message("DAPI not in expression columns — adding to skip list")
    } else {
      message("DAPI found in expression columns — NOT skipping DAPI")
    }
    
    # clustering resolutions:
    resolutions = clustering_res
    # which subset to store (by default all)
    idx=seq(length(resolutions))
    
    # Where output dirs should be housed:
    workDir=out_dir
    print(workDir)
    
    out.label=NULL
    
    data_csv=segmentation_file
    message(glue('Reading file: {data_csv}'))
    message(glue('Skipping columns: {skip_cols}'))
    
    out_dir = glue("{out_dir}/out-phenomenalist")
    if (!dir.exists(out_dir)) { dir.create(out_dir) }
    print('creating nested subdir')
    out_dir = paste0(out_dir, "/", label)
    dir.create(out_dir,showWarnings = F)
    message(out_dir)
    
    if('lock.rds' %in% list.files(out_dir)){
      spe = readRDS(glue('{out_dir}/spe.rds'))
      spe = cluster.mod(spe, resolution = resolutions, out_dir = out_dir)
      names(colData(spe))
      saveRDS(spe,glue(out_dir,'/spe.rds'))
      
      cluster_cols = stringr::str_subset(names(colData(spe)), "cluster_leiden")
      
      for(i in seq(length(cluster_cols))){
        plot_heatmap.mod(spe, group_by = cluster_cols[i],auto = T,out_dir = out_dir,segment=T)
      }
      
      res_found=do.call('rbind',strsplit(cluster_cols,'[cluster_leiden_res]'))
      prepare_mask_inputs(spe=spe,out_dir = out_dir,res = max(as.numeric(res_found[,19])),mask.only = mask.only,label=label)
      
    }else{
      incProgress(1/6, detail = 'Generating Phenomenalist object')
      spe <- create_object.mod(data_csv, skip_cols = skip_cols, transformation = "z", out_dir = out_dir,expression_cols = expression.columns,classifier.label = classifier_label,min.cells = min.cells,max.cells = max.cells)

      if ("tile_num" %in% names(colData(spe))) plot_spatial(spe, color_by = "tile_num",out_dir = out_dir)
      if ("Classifier_Label" %in% names(colData(spe))) plot_spatial(spe, color_by = "Classifier_Label",out_dir = out_dir)
      
      incProgress(1/6, detail = 'Generating spatial objects')
      spatial_dir = paste0(out_dir, "/spatial-expression")
      dir.create(spatial_dir, showWarnings = FALSE)
      plot_spatial(spe, color_by = names(spe), out_dir = spatial_dir)
      plot_spatial(spe, color_by = names(spe), smooth = TRUE, out_dir = spatial_dir)
      
      incProgress(1/6, detail = 'Generating dimensionality reductions')
      umap_dir = paste0(out_dir, "/UMAP-expression")
      dir.create(umap_dir, showWarnings = FALSE)
      plot_dr(spe, dr = "UMAP", color_by = names(spe), out_dir = umap_dir)
      plot_dr(spe, dr = "UMAP", color_by = names(spe), smooth = TRUE, out_dir = umap_dir)
      
      incProgress(1/6, detail = 'Initiating clustering')
      spe = cluster.mod(spe, resolution = resolutions, out_dir = out_dir)
      names(colData(spe))
      write_rds(spe,glue(out_dir,'/spe.rds'))
      cluster_cols = stringr::str_subset(names(colData(spe)), "cluster_leiden")
      
      for(i in seq(length(cluster_cols))){
        plot_heatmap.mod(spe, group_by = cluster_cols[i],auto = T,out_dir = out_dir,segment=T)
      }
      
      res_found=do.call('rbind',strsplit(cluster_cols,'[cluster_leiden_res]'))
      prepare_mask_inputs(spe=spe,out_dir = out_dir,res = max(as.numeric(res_found[,19])),mask.only = mask.only,label=label)
      
    }
    incProgress(1/6, detail = 'Done')
    
  })
}
