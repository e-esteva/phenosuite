plot_heatmap.mod=function(x, group_by, assay = "logcounts", out_dir = NULL,size.col=4,size.row=4,auto=F,segment=F,features=NULL){
  #remotes::install_version("matrixStats", version="1.1.0")
  
  if (!is(x, "SpatialExperiment")) {
    stop("input is not a SpatialExperiment object")
  }
  if (!is.character(group_by)) {
    stop("`group_by` is not a character string")
  }
  if (!all(group_by %in% names(colData(x)))) {
    stop("not all `group_by` values are present in the object")
  }
  if (!is.character(assay)) {
    stop("`assay` is not a character string")
  }
  if (!assay %in% assayNames(x)) {
    stop("input SpatialExperiment object does not have a `", 
         assay, "` assay")
  }
  if (!is.null(out_dir)) {
    if (!dir.exists(out_dir)) {
      stop("output directory `", out_dir, "` does not exist")
    }
  }
  gradient_colors <- rev(RColorBrewer::brewer.pal(11, "RdYlBu"))
  for (g in group_by) {
    if(!is.null(features)){
	data_ = assay(x,assay)
	e <- scuttle::summarizeAssayByGroup(list(counts = data_[features[features %in% row.names(data_)],]), ids = colData(x)[[g]], statistics = "median")
    	    

    }else{
    	e <- scuttle::summarizeAssayByGroup(x, ids = colData(x)[[g]], statistics = "median")
    
    }
    e <- SummarizedExperiment::assay(e, i = "median")
    
    if(sum(rowSums(e) == 0) > 0){
      e <- scuttle::summarizeAssayByGroup(x, ids = colData(x)[[g]],statistics = "mean")
      e <- SummarizedExperiment::assay(e, i = "mean")

    }

    e <- scale(t(e))
    e[e > 5] <- 5
    e[e < -5] <- -5
    if(auto){
      size.col=8
      cluster.size=length(unique(x[[g]]))
      if(cluster.size <= 20){
        size.row=8
        
        hm <- ComplexHeatmap::Heatmap(e, name = "Expression", 
                                      row_title = g, col = gradient_colors, cluster_rows = TRUE, 
                                      cluster_columns = TRUE,column_names_gp = grid::gpar(fontsize=size.col),row_names_gp = grid::gpar(fontsize=size.row))
        if (is.null(out_dir)) {
          return(hm)
        }
        else {
          out_base <- glue("{out_dir}/{g}-heatmap")
          png(filename = glue("{out_base}.png"), width = 10, 
              height = 5, units = "in", res = 300)
          ComplexHeatmap::draw(hm)
          dev.off()
        }
        
      }else{
        if(segment){
          clust.increments=c(seq(1,cluster.size,20),cluster.size)
          
          if(cluster.size <= 30){
            size.row=3
          }else{
            size.row =2 
          }
          
          for(i in seq(2,length(clust.increments))){
            e.tmp=e[seq(clust.increments[i-1],clust.increments[i]),]
            
            hm <- ComplexHeatmap::Heatmap(e.tmp, name = "Expression", 
                                          row_title = g, col = gradient_colors, cluster_rows = TRUE, 
                                          cluster_columns = TRUE,column_names_gp = grid::gpar(fontsize=size.col),row_names_gp = grid::gpar(fontsize=size.row))
            
            if (is.null(out_dir)) {
              return(hm)
            }
            else {
              out_base <- glue("{out_dir}/{g}-{i-1}-heatmap")
              png(filename = glue("{out_base}.png"), width = 10, 
                  height = 5, units = "in", res = 300)
              ComplexHeatmap::draw(hm)
              dev.off()
            }
            
          }
        }else{
          
          
          if(cluster.size <= 30){
            size.row=3
          }else{
            size.row =2 
          }
          hm <- ComplexHeatmap::Heatmap(e, name = "Expression",
                                        row_title = g, col = gradient_colors, cluster_rows = TRUE,
                                        cluster_columns = TRUE,column_names_gp = grid::gpar(fontsize=size.col),row_names_gp = grid::gpar(fontsize=size.row))
          
          if (is.null(out_dir)) {
            return(hm)
          }
          else {
            out_base <- glue("{out_dir}/{g}-heatmap")
            png(filename = glue("{out_base}.png"), width = 10,
                height = 5, units = "in", res = 300)
            ComplexHeatmap::draw(hm)
            dev.off()
          }
        }
        
        
      }
      
    }else{
      hm <- ComplexHeatmap::Heatmap(e, name = "Expression",
                                    row_title = g, col = gradient_colors, cluster_rows = TRUE,
                                    cluster_columns = TRUE,column_names_gp = grid::gpar(fontsize=size.col),row_names_gp = grid::gpar(fontsize=size.row))
      
      if (is.null(out_dir)) {
        return(hm)
      }
      else {
        out_base <- glue("{out_dir}/{g}-heatmap")
        png(filename = glue("{out_base}.png"), width = 10,
            height = 5, units = "in", res = 300)
        ComplexHeatmap::draw(hm)
        dev.off()
      }
    }
    
  }
}
create_object.mod.v0=function (x, expression_cols = NULL, metadata_cols = NULL, skip_cols = NULL, 
                            clean_names = TRUE, transformation = NULL, out_dir = NULL,down_sample_idx=NULL,classifier.label=NULL,plot_diagnostic=F,min.cells=1000,max.cells=1e5){
  if (is.character(x)) {
    x <- data.table::fread(x, stringsAsFactors = FALSE, data.table = FALSE)
    n1=x$`Classifier.Label`
    n2=x$`Classifier Label`
    if(length(str_subset(pattern='Classifier',string=names(x))) > 0){
    	if(!is.null(classifier.label)){
      		if(is.null(n2)){
        		hits=rowSums(do.call('cbind',lapply(classifier.label,function(y) str_detect(x$`Classifier.Label`,y))))
      		}else{
        		hits=rowSums(do.call('cbind',lapply(classifier.label,function(y) str_detect(n2,y))))
      		}
      
      		x=x[hits==1,]
    	}
    }
    if(!is.null(down_sample_idx) || dim(x)[1] > max.cells){
      if(is.null(down_sample_idx)){
        down_sample_idx=sample(dim(x)[1],max.cells)
      }
      x=x[down_sample_idx,]
    }
  }
  if (!is.data.frame(x)) {
    stop("input is not a file or a data frame")
  }
  # if (!is.null(out_dir)) {
  #   if (dir.exists(out_dir)) {
  #     stop("output directory already exists")
  #   }
  # }
  if (nrow(x) < min.cells) {
    print(x)
    print(nrow(x))
    stop("data frame has too few rows/cells")
  }
  if (ncol(x) < 10) {
    stop("data frame has too few columns")
  }
  message("number of input table rows: ", nrow(x))
  message("number of input table columns: ", ncol(x))
  message("")
  if (!is.null(skip_cols)) {
    if (length(skip_cols) == 1) {
      skip_cols <- str_subset(names(x), pattern = skip_cols)
    }
    x <- x[, setdiff(names(x), skip_cols)]
    message("skipped columns: ", toString(skip_cols), "\n")
  }
  if (is.null(expression_cols)) {
    expression_cols <- detect_exprs_cols(x)
  }
  
  if (length(expression_cols) == 1) {
    expression_cols <- str_subset(names(x), pattern = expression_cols)
  }else{
    expression_cols=paste0(expression_cols,collapse = "|")
    expression_cols <- str_subset(names(x), pattern = expression_cols)
  }
  if (length(setdiff(expression_cols, names(x))) > 0) {
    stop("missing markers: ", toString(setdiff(expression_cols, 
                                               names(x))))
  }
  exprs <- x[, expression_cols]
  if (is.null(metadata_cols)) {
    metadata_cols <- setdiff(names(x), expression_cols)
  }
  if (length(metadata_cols) == 1) {
    metadata_cols <- str_subset(names(x), pattern = metadata_cols)
  }
  if (length(setdiff(metadata_cols, names(x))) > 0) {
    stop("missing metadata columns: ", toString(setdiff(metadata_cols, 
                                                        names(x))))
  }
  message('metadata cols:')
  message(metadata_cols)
  message(class(x))
  message(dim(x)[1])
  message(dim(x)[2])  
  x <- x[, metadata_cols]
  message(class(x))
  message(dim(x)[1])
  message(dim(x)[2])
  message('exprs:')
  message(expression_cols)
  message(class(exprs))
  message(dim(exprs)[1])
  message(dim(exprs)[2])
  if (clean_names) {
    
    clean_col_names <- function(x) {
      
      # check if the input is valid
      if (!is.data.frame(x)) {
        stop("input is not a data frame")
      }
      
      # adjust MAV column names
      names(x) <- str_remove(names(x), ":Cyc_.*")
      names(x) <- str_remove(names(x), ".*:")
      
      # adjust Visiopharm column names
      names(x) <- str_remove(names(x), "MP.*\\(")
      
      # adjust HALO column names
      names(x) <- str_remove(names(x), " Cell Intensity")
      
      # clean column names
      x <- janitor::clean_names(x, case = "none")
      
      # adjust common column names
      if (!"x" %in% names(x)) {
        names(x)[names(x) == "X"] <- "x"
        names(x)[names(x) == "X_X"] <- "x"
        names(x)[names(x) == "CtrX"] <- "x"
        if (all(c("XMin", "XMax") %in% names(x))) {
          x$x <- (x$XMin + x$XMax) / 2
        }
	if (all(c("x_min", "x_max") %in% names(x))) {
          x$x <- (x$x_min + x$x_max) / 2
       	}
      }
      if (!"y" %in% names(x)) {
        names(x)[names(x) == "Y"] <- "y"
        names(x)[names(x) == "Y_Y"] <- "y"
        names(x)[names(x) == "CtrY"] <- "y"
        if (all(c("YMin", "YMax") %in% names(x))) {
          x$y <- (x$YMin + x$YMax) / 2
        }
	if (all(c("y_min", "y_max") %in% names(x))) {
          x$y <- (x$y_min + x$y_max) / 2
       	}
      }
      if (!"z" %in% names(x)) {
        names(x)[names(x) == "Z"] <- "z"
        names(x)[names(x) == "Z_Z"] <- "z"
      }
      if (!"reg" %in% names(x)) {
        names(x)[names(x) == "Region"] <- "reg"
        names(x)[names(x) == "region"] <- "reg"
      }
      if (!"tile_num" %in% names(x)) {
        names(x)[names(x) == "tile_nr"] <- "tile_num"
        names(x)[names(x) == "tile_number"] <- "tile_num"
      }
      if (!"size" %in% names(x)) {
        names(x)[names(x) == "Cell_Area_mm2"] <- "size"
      }
      
      return(x)
    }
    
    exprs <- clean_col_names(exprs)
    x <- clean_col_names(x)
    
    if (!"cell_id" %in% names(x)) {
      names(x)[names(x) == "CellID"] <- "cell_id"
      names(x)[names(x) == "label"] <- "cell_id"
      names(x)[names(x) == "Object_Id"] <- "cell_id"
    }
    if (!"cell_id" %in% names(x)) {
      x$cell_id <- rownames(x)
    }
    if (is.numeric(x$cell_id)) {
      x$cell_id <- str_pad(as.character(x$cell_id), width = 7, 
                           pad = "0")
      x$cell_id <- str_c("C", x$cell_id)
    }
    x$cell_id <- make.names(x$cell_id, unique = TRUE)
  }
  if (!"cell_id" %in% names(x)) {
    stop("data frame must contain `cell_id` column")
  }
  if (!"x" %in% names(x)) {
    stop("data frame must contain `x` column")
  }
  if (!"y" %in% names(x)) {
    stop("data frame must contain `y` column")
  }
  rownames(x) <- x$cell_id
  exprs <- as.matrix(exprs)
  exprs <- exprs[, sort(colnames(exprs))]
  rownames(exprs) <- rownames(x)
  message("number of expression columns: ", ncol(exprs))
  message("number of metadata columns: ", ncol(x))
  message("")
  message("expression columns: ", toString(colnames(exprs)), 
          "\n")
  message("metadata columns: ", toString(colnames(x)), "\n")
  # if (!is.null(out_dir)) {
  #   dir.create(out_dir)
  # }
  if(plot_diagnostic){
    expression_dir=glue('{out_dir}/expression-diagnostic-plots')
    for(i in seq(dim(exprs)[2])){
      pdf(glue('{out_dir}/{colnames(exprs)[i]}.pdf'))
      plot(density(exprs[,i],main=colnames(exprs)[i]))
      abline(v=quantile(exprs[,i]))
      dev.off()
    }
  }
  
  na.idx=seq(dim(t(exprs))[1])[rowVars(t(exprs),useNames = T)==0]
  genes0var=paste0(colnames(exprs)[na.idx],collapse=',')
  
  if(length(na.idx) > 0){
    exprs=exprs[,-na.idx]
    
    message(glue('removed {length(na.idx)} genes due to 0 variance'))
    message(glue("removed: {genes0var}"))
    
  }
  
  s <- SpatialExperiment::SpatialExperiment(assay = list(counts = t(exprs)), 
                                            colData = x, spatialCoordsNames = c("x", "y"))
  if (!is.null(transformation)) {
    s <- transform(s, method = transformation, out_dir = out_dir)
    s <- run_umap(s, n_threads = 4)
  }
  if (!is.null(out_dir)) {
    message("saving object")
    write_rds(s,file=paste0(out_dir, "/spe.rds"))
    #saveRDS(s, paste0(out_dir, "/spe.rds"))
  }
  return(s)
}
create_object.mod=function (x, expression_cols = NULL, metadata_cols = NULL, skip_cols = NULL, 
                            clean_names = TRUE, transformation = NULL, out_dir = NULL,down_sample_idx=NULL,classifier.label=NULL,plot_diagnostic=F,min.cells=1000,max.cells=1e5){
  if (is.character(x)) {
    x <- data.table::fread(x, stringsAsFactors = FALSE, data.table = FALSE)
    n1=x$`Classifier.Label`
    n2=x$`Classifier Label`
    if(length(str_subset(pattern='Classifier',string=names(x))) > 0){
    	if(!is.null(classifier.label)){
      		if(is.null(n2)){
        		hits=rowSums(do.call('cbind',lapply(classifier.label,function(y) str_detect(x$`Classifier.Label`,y))))
      		}else{
        		hits=rowSums(do.call('cbind',lapply(classifier.label,function(y) str_detect(n2,y))))
      		}
      
      		x=x[hits==1,]
    	}
    }
    if(!is.null(down_sample_idx) || dim(x)[1] > max.cells){
      if(is.null(down_sample_idx)){
        down_sample_idx=sample(dim(x)[1],max.cells)
      }
      x=x[down_sample_idx,]
    }
  }
  if (!is.data.frame(x)) {
    stop("input is not a file or a data frame")
  }
  if (nrow(x) < min.cells) {
    print(x)
    print(nrow(x))
    stop("data frame has too few rows/cells")
  }
  if (ncol(x) < 10) {
    stop("data frame has too few columns")
  }
  message("number of input table rows: ", nrow(x))
  message("number of input table columns: ", ncol(x))
  message("")
  if (!is.null(skip_cols)) {
    if (length(skip_cols) == 1) {
      skip_cols <- str_subset(names(x), pattern = skip_cols)
    }
    x <- x[, setdiff(names(x), skip_cols)]
    message("skipped columns: ", toString(skip_cols), "\n")
  }
  if (is.null(expression_cols)) {
    expression_cols <- detect_exprs_cols(x)
  }
  
  if (length(expression_cols) == 1) {
    expression_cols <- str_subset(names(x), pattern = expression_cols)
  }else{
    expression_cols=paste0(expression_cols,collapse = "|")
    expression_cols <- str_subset(names(x), pattern = expression_cols)
  }
  if (length(setdiff(expression_cols, names(x))) > 0) {
    stop("missing markers: ", toString(setdiff(expression_cols, 
                                               names(x))))
  }
  exprs <- x[, expression_cols]
  if (is.null(metadata_cols)) {
    metadata_cols <- setdiff(names(x), expression_cols)
  }
  if (length(metadata_cols) == 1) {
    metadata_cols <- str_subset(names(x), pattern = metadata_cols)
  }
  if (length(setdiff(metadata_cols, names(x))) > 0) {
    stop("missing metadata columns: ", toString(setdiff(metadata_cols, 
                                                        names(x))))
  }
  message('metadata cols:')
  message(metadata_cols)
  message(class(x))
  message(dim(x)[1])
  message(dim(x)[2])  
  x <- x[, metadata_cols]
  message(class(x))
  message(dim(x)[1])
  message(dim(x)[2])
  message('exprs:')
  message(expression_cols)
  message(class(exprs))
  message(dim(exprs)[1])
  message(dim(exprs)[2])
  if (clean_names) {
    
    clean_col_names <- function(x) {
      
      # check if the input is valid
      if (!is.data.frame(x)) {
        stop("input is not a data frame")
      }
      
      # adjust MAV column names
      names(x) <- str_remove(names(x), ":Cyc_.*")
      names(x) <- str_remove(names(x), ".*:")
      
      # adjust Visiopharm column names
      names(x) <- str_remove(names(x), "MP.*\\(")
      
      # adjust HALO column names
      names(x) <- str_remove(names(x), " Cell Intensity")
      
      # clean column names
      x <- janitor::clean_names(x, case = "none")
      
      # --- Resolve X coordinate ---
      if (!"x" %in% names(x)) {
        names(x)[names(x) == "X"] <- "x"
        names(x)[names(x) == "X_X"] <- "x"
        names(x)[names(x) == "CtrX"] <- "x"
        
        # Centroid detection (case-insensitive): matches e.g.
        # "Centroid X", "Centroid_X", "centroid.x", "Centroid_X_px", "XCentroid"
        # X/Y must be isolated (bounded by separator or string edge) to avoid
        # false matches like "Centroid_Y_px" where 'x' appears inside 'px'
        if (!"x" %in% names(x)) {
          centroid_x <- grep("centroid[_.\\s]*x([_.\\s]|$)|(^|[_.\\s])x[_.\\s]*centroid", 
                             names(x), ignore.case = TRUE, value = TRUE)
          if (length(centroid_x) > 0) {
            names(x)[names(x) == centroid_x[1]] <- "x"
            message("Mapped centroid column '", centroid_x[1], "' -> x")
          }
        }
        
        # Bounding box fallback
        if (!"x" %in% names(x)) {
          if (all(c("XMin", "XMax") %in% names(x))) {
            x$x <- (x$XMin + x$XMax) / 2
          }
          if (all(c("x_min", "x_max") %in% names(x))) {
            x$x <- (x$x_min + x$x_max) / 2
          }
        }
      }
      
      # --- Resolve Y coordinate ---
      if (!"y" %in% names(x)) {
        names(x)[names(x) == "Y"] <- "y"
        names(x)[names(x) == "Y_Y"] <- "y"
        names(x)[names(x) == "CtrY"] <- "y"
        
        # Centroid detection (case-insensitive): matches e.g.
        # "Centroid Y", "Centroid_Y", "centroid.y", "Centroid_Y_px", "YCentroid"
        if (!"y" %in% names(x)) {
          centroid_y <- grep("centroid[_.\\s]*y([_.\\s]|$)|(^|[_.\\s])y[_.\\s]*centroid", 
                             names(x), ignore.case = TRUE, value = TRUE)
          if (length(centroid_y) > 0) {
            names(x)[names(x) == centroid_y[1]] <- "y"
            message("Mapped centroid column '", centroid_y[1], "' -> y")
          }
        }
        
        # Bounding box fallback
        if (!"y" %in% names(x)) {
          if (all(c("YMin", "YMax") %in% names(x))) {
            x$y <- (x$YMin + x$YMax) / 2
          }
          if (all(c("y_min", "y_max") %in% names(x))) {
            x$y <- (x$y_min + x$y_max) / 2
          }
        }
      }
      
      # --- Resolve Z coordinate ---
      if (!"z" %in% names(x)) {
        names(x)[names(x) == "Z"] <- "z"
        names(x)[names(x) == "Z_Z"] <- "z"
      }
      if (!"reg" %in% names(x)) {
        names(x)[names(x) == "Region"] <- "reg"
        names(x)[names(x) == "region"] <- "reg"
      }
      if (!"tile_num" %in% names(x)) {
        names(x)[names(x) == "tile_nr"] <- "tile_num"
        names(x)[names(x) == "tile_number"] <- "tile_num"
      }
      if (!"size" %in% names(x)) {
        names(x)[names(x) == "Cell_Area_mm2"] <- "size"
      }
      
      return(x)
    }
    
    exprs <- clean_col_names(exprs)
    x <- clean_col_names(x)
    
    if (!"cell_id" %in% names(x)) {
      names(x)[names(x) == "CellID"] <- "cell_id"
      names(x)[names(x) == "label"] <- "cell_id"
      names(x)[names(x) == "Object_Id"] <- "cell_id"
    }
    if (!"cell_id" %in% names(x)) {
      x$cell_id <- rownames(x)
    }
    if (is.numeric(x$cell_id)) {
      x$cell_id <- str_pad(as.character(x$cell_id), width = 7, 
                           pad = "0")
      x$cell_id <- str_c("C", x$cell_id)
    }
    x$cell_id <- make.names(x$cell_id, unique = TRUE)
  }
  if (!"cell_id" %in% names(x)) {
    stop("data frame must contain `cell_id` column")
  }
  if (!"x" %in% names(x)) {
    stop("data frame must contain `x` column")
  }
  if (!"y" %in% names(x)) {
    stop("data frame must contain `y` column")
  }
  rownames(x) <- x$cell_id
  exprs <- as.matrix(exprs)
  exprs <- exprs[, sort(colnames(exprs))]
  rownames(exprs) <- rownames(x)
  message("number of expression columns: ", ncol(exprs))
  message("number of metadata columns: ", ncol(x))
  message("")
  message("expression columns: ", toString(colnames(exprs)), 
          "\n")
  message("metadata columns: ", toString(colnames(x)), "\n")
  if(plot_diagnostic){
    expression_dir=glue('{out_dir}/expression-diagnostic-plots')
    for(i in seq(dim(exprs)[2])){
      pdf(glue('{out_dir}/{colnames(exprs)[i]}.pdf'))
      plot(density(exprs[,i],main=colnames(exprs)[i]))
      abline(v=quantile(exprs[,i]))
      dev.off()
    }
  }
  
  na.idx=seq(dim(t(exprs))[1])[rowVars(t(exprs),useNames = T)==0]
  genes0var=paste0(colnames(exprs)[na.idx],collapse=',')
  
  if(length(na.idx) > 0){
    exprs=exprs[,-na.idx]
    
    message(glue('removed {length(na.idx)} genes due to 0 variance'))
    message(glue("removed: {genes0var}"))
    
  }
  
  s <- SpatialExperiment::SpatialExperiment(assay = list(counts = t(exprs)), 
                                            colData = x, spatialCoordsNames = c("x", "y"))
  if (!is.null(transformation)) {
    s <- transform(s, method = transformation, out_dir = out_dir)
    s <- run_umap(s, n_threads = 4)
  }
  if (!is.null(out_dir)) {
    message("saving object")
    write_rds(s,file=paste0(out_dir, "/spe.rds"))
  }
  return(s)
}
cluster.mod=function (x, method = c("leiden"), resolution = 1, n_neighbors = 50, 
                      out_dir = NULL,max_clust=500,label=NULL) 
{
  method <- match.arg(method)
  if (!is(x, "SpatialExperiment")) {
    stop("input is not a SpatialExperiment object")
  }
  if (!is.numeric(resolution)) {
    stop("`resolution` is not a number")
  }
  if (!is.numeric(n_neighbors)) {
    stop("`n_neighbors` is not a number")
  }
  if (!"exprs" %in% assayNames(x)) {
    stop("input SpatialExperiment object does not have a `exprs` assay")
  }
  if (!is.null(out_dir)) {
    if (!dir.exists(out_dir)) {
      stop("output directory `", out_dir, "` does not exist")
    }
    clusters_dir <- glue("{out_dir}/clusters")
    dir.create(clusters_dir)
  }
  exprs_mat <- assay(x, "exprs")
  if (method == "leiden") {
    g <- scran::buildSNNGraph(exprs_mat, transposed = FALSE, 
                              k = n_neighbors)
    n_clust_prev <- 0
    for (res_num in resolution) {

      incProgress((1/6)/length(resolution), detail = glue('clustering @ resolution {res_num}'))
      message(glue("clustering using resolution of {res_num}"))
      set.seed(99)
      clusters <- igraph::cluster_leiden(g, objective_function = "modularity", 
                                         resolution_parameter = res_num, n_iterations = 10)
      clusters <- clusters$membership
      res_str <- format(as.numeric(res_num), nsmall = 1)
      res_str <- stringr::str_pad(res_str, width = 3, side = "left", 
                                  pad = "0")
      res_str <- stringr::str_c("res", res_str)
      n_clust <- length(unique(clusters))
      if(n_clust <= max_clust){
	      if(is.null(label)){
        		clusters_label <- glue("cluster_{method}_{res_str}_clust{n_clust}")
        	}else{
			clusters_label <- glue("cluster_{method}_{res_str}_clust{n_clust}_{label}")
		}
	clusters <- as.character(clusters)
        clusters <- stringr::str_pad(clusters, width = 2, 
                                     side = "left", pad = "0")
        clusters <- stringr::str_c("C", clusters)
        
        if (n_clust > n_clust_prev) {
          x[[clusters_label]] <- factor(clusters)
          n_clust_prev <- n_clust
          if (!is.null(out_dir)) {
            cluster_summary <- janitor::tabyl(x[[clusters_label]])
            colnames(cluster_summary) <- c("cluster", "cells_num", 
                                           "cells_freq")
            readr::write_csv(cluster_summary, glue("{clusters_dir}/{clusters_label}-summary.csv"))
            plot_heatmap(x, group_by = clusters_label, 
                         out_dir = clusters_dir)
            plot_spatial(x, color_by = clusters_label, 
                         out_dir = clusters_dir)
            plot_dr(x, dr = "UMAP", color_by = clusters_label, 
                    out_dir = clusters_dir)
          }
          
        }
      }else{
         
     	message(glue('{n_clust} clusters found, ending clustering'))
        if(!is.null(out_dir)){
          saveRDS(glue('clusters:{res_num}'),glue('{out_dir}/lock.rds'))
	}
	break
      }
    }
  }
  if (!is.null(out_dir)) {
    message("saving object")
    saveRDS(x, paste0(out_dir, "/spe.rds"))
  }
  return(x)
}
phenomenalist.preprocess=function(x,failed.markers=NULL,nuclear.markers=NULL,else.cytoplasm=F,HALO=T){
  if(is.null(failed.markers) & is.null(nuclear.markers)){
    return(NULL)
  }else{
    if(HALO){
      if(is.null(nuclear.markers)){
        
        x=data.table::fread(x)
        
        names_original=names(x)
        names(x)=tolower(names(x))
        failed.markers=tolower(failed.markers)
        
        cell.intensity=names(x)[grep(tolower('Cell.Intensity'),names(x))]
        cell.intensity=cell.intensity[-unlist(lapply(failed.markers,function(marker) grep(marker,cell.intensity)))]
        expression.columns=cell.intensity
        expression.columns = names_original[match(expression.columns,names(x))]
        print(expression.columns)
      }else{
        message("selected job: nucleus else aggregate")
        message('nuclear markers:')
        message(nuclear.markers)
        message('failed markers: ')
        message(failed.markers)
        
        nuclear_markers=nuclear.markers
        
        x=data.table::fread(x)
        names_original=names(x)
        names(x)=tolower(names(x))
        nuclear_markers=tolower(nuclear_markers)
        failed.markers=tolower(failed.markers)
        
        expression.columns=names(x)[unlist(lapply(nuclear_markers,function(m) grep(m,names(x))))]
        
        nuclear.intensity.cols=expression.columns[grep(tolower('Nucleus Intensity'),expression.columns)]
        if(length(nuclear.intensity.cols) == 0){
          nuclear.intensity.cols=expression.columns[grep(tolower('Nucleus.Intensity'),expression.columns)]
        }
        cell.intensity=names(x)[grep(tolower('Cell Intensity'),names(x))]
        if(length(cell.intensity) == 0){
          cell.intensity=names(x)[grep(tolower('Cell.Intensity'),names(x))]
        }
        message('full cell intensity:')
        message(cell.intensity)
        message(c(nuclear_markers,failed.markers))
	dots_detected=F
        nuclear_markers=do.call('rbind',strsplit(nuclear.intensity.cols,' nucleus intensity'))[,1]
        if(all.equal(nuclear_markers,nuclear.intensity.cols)){
		dots_detected=T
		nuclear_markers=do.call('rbind',strsplit(nuclear.intensity.cols,'.nucleus.intensity'))[,1]
	}
    	message('nuclear markers:')
        message(nuclear_markers)
	if(dots_detected){
		cell.intensity=cell.intensity[-unlist(lapply(c(nuclear_markers,failed.markers),function(marker) grep(glue('^{tolower(marker)}.cell.intensity'),cell.intensity)))]
	}else{
        	cell.intensity=cell.intensity[-unlist(lapply(c(nuclear_markers,failed.markers),function(marker) grep(glue('^{tolower(marker)} cell intensity'),cell.intensity)))]
        }
	message('reduced cell intensity:')
        message(cell.intensity)
        job2.selected_expression_cols=c(nuclear.intensity.cols,cell.intensity)
        if(else.cytoplasm){
          cytoplasm.intensity=names(x)[grep(tolower('Cytoplasm Intensity'),names(x))]
          if(length(cytoplasm.intensity) == 0){
            cytoplasm.intensity=names(x)[grep(tolower('Cytoplasm.Intensity'),names(x))]
          }
          cytoplasm.intensity.nuc_failed=lapply(c(nuclear_markers,failed.markers),function(marker) grep(glue('^{marker}$'),cytoplasm.intensity))
          if(length(unlist(cytoplasm.intensity.nuc_failed)) > 0){
            cytoplasm.intensity=cytoplasm.intensity[-unlist(cytoplasm.intensity.nuc_failed)]
          }
          
          job3.selected_expression_cols=c(nuclear.intensity.cols,cytoplasm.intensity)
          expression.columns=job3.selected_expression_cols
        }else{
          expression.columns=job2.selected_expression_cols
          print(expression.columns)
        }
        expression.columns = names_original[match(expression.columns,names(x))]
      }
      
    }else{
      x=data.table::fread(x)
      # if('V1' %in% colnames(x)){
      #   init.idx=6
      # }else{
      #   init.idx=5
      # }
      # ## or just look for first hit of Ch: (Ch1Cy1 likely always first)
      # if(length(grep("^Ch",colnames(x))) > 0){
      #   init.idx=min(grep("^Ch",colnames(x)))
      # }
      init.idx=max(na.omit(match(c('label',     'y',     'x', 'size', 'y_min', 'x_min', 'y_max', 'x_max'),colnames(x))))+1
      message(glue('init.idx: {init.idx}'))
      cell.intensity=colnames(x)[seq(init.idx,dim(x)[2])]
      message(glue('cell.intensity: {cell.intensity}'))

      if(!is.null(failed.markers)){
        cell.intensity=cell.intensity[-unlist(lapply(failed.markers,function(marker) grep(marker,cell.intensity)))]
      }
      
      expression.columns=cell.intensity
    }
    
    return(expression.columns)
    
  }
}
select_intensity_columns <- function(filepath,
                                     failed_markers = NULL,
                                     nuclear_markers = NULL,
                                     prefer_cytoplasm = FALSE) {

  x <- data.table::fread(filepath)
  original_names <- names(x)
  names_lower <- tolower(original_names)

  # Helper: find columns matching a compartment keyword (case-insensitive)
  find_compartment <- function(compartment) {
    which(grepl(compartment, names_lower))
  }

  # Helper: find which column indices match any of a set of marker names
  marker_indices <- function(idx_pool, markers) {
    if (is.null(markers) || length(markers) == 0) return(integer(0))
    markers_lower <- tolower(markers)
    unlist(lapply(markers_lower, function(m) {
      idx_pool[grepl(m, names_lower[idx_pool])]
    }))
  }

  # Identify compartment column pools
  cell_idx     <- find_compartment("cell.*intensity|intensity.*cell")
  nucleus_idx  <- find_compartment("nucleus.*intensity|intensity.*nucleus")
  cyto_idx     <- find_compartment("cytoplasm.*intensity|intensity.*cytoplasm")

  # All intensity columns as a fallback

  all_intensity_idx <- find_compartment("intensity")

  if (is.null(nuclear_markers)) {
    # Job 1: use Cell Intensity for everything, drop failed markers
    selected_idx <- if (length(cell_idx) > 0) cell_idx else all_intensity_idx
    selected_idx <- setdiff(selected_idx, marker_indices(selected_idx, failed_markers))

  } else {
    # Job 2/3: nuclear markers → Nucleus Intensity; others → Cell or Cytoplasm
    nuc_selected <- marker_indices(nucleus_idx, nuclear_markers)

    # Remove nuclear + failed markers from the "other" compartment
    drop_markers <- c(nuclear_markers, failed_markers)

    if (prefer_cytoplasm && length(cyto_idx) > 0) {
      other_selected <- setdiff(cyto_idx, marker_indices(cyto_idx, drop_markers))
    } else {
      other_selected <- setdiff(cell_idx, marker_indices(cell_idx, drop_markers))
    }

    selected_idx <- c(nuc_selected, other_selected)
  }

  # Return original-cased column names
  original_names[selected_idx]
}

prepare_mask_inputs=function(spe,mask.only,out_dir,res,failed.markers,label=NULL){
  library(glue)
  message(glue('generating mask-generation inputs for resolution {res} clusters'))
  
  spatial_obj=data.frame(spatialCoords(spe))
  res.hits=str_detect(names(colData(spe)),glue('res{res}'))
  hits.found=ifelse(sum(res.hits)>0,T,F)
  while(hits.found==F){
    res=res-1
    res.hits=str_detect(names(colData(spe)),glue('res{res}'))
    hits.found=ifelse(sum(res.hits)>0,T,F)
  }
  clust.tmp=colData(spe)[[names(colData(spe))[res.hits]]]
  spatial_obj$cluster=clust.tmp
  
  mask_inputs_dir = paste0(out_dir, "/mask-inputs/")
  dir.create(mask_inputs_dir, showWarnings = FALSE)
  if(!is.null(label)){
    out_name=glue('{mask_inputs_dir}/{label}-spatial-anno-res{res}.csv')
  }else{
    out_name=glue('{mask_inputs_dir}/spatial-anno-res{res}.csv')
  }
  write.csv(spatial_obj,out_name)
  
  message('generating marker metadata')
  all.markers=row.names(spe)
  all.markers=do.call('rbind',strsplit(all.markers,'*_Cytoplasm|_Nucleus'))[,1]
  if(!is.null(mask.only)){
    working.markers=c(all.markers,mask.only)
  }else{
    working.markers=all.markers
  }
  if(!is.null(label)){
    out_name=glue('{mask_inputs_dir}/{label}-marker-metadata.csv')
  }else{
    out_name=glue('{mask_inputs_dir}/marker-metadata.csv')
  }
  write.csv(working.markers,out_name)
  
}
#source('/Users/ee699/working/TRIC/phenomenalist/R/plot-scatter.R')
source('/srv/shiny-server/phenomenalist/utils/plot-scatter.R')
plot_dr.mod=function (x, dr, color_by, assay = "logcounts", smooth = FALSE, 
                      range = c(0.01, 0.99), out_dir = NULL,h=NULL,w=NULL,pdf=F,interactive=F) 
{
  if (!is(x, "SpatialExperiment")) {
    stop("input is not a SpatialExperiment object")
  }
  if (!is.character(dr)) {
    stop("`dr` is not a character string")
  }
  if (!dr %in% reducedDimNames(x)) {
    stop("SpatialExperiment object does not have a dimensionality reduction `", 
         dr, "`")
  }
  if (!is.null(out_dir)) {
    if (!dir.exists(out_dir)) {
      stop("output directory `", out_dir, "` does not exist")
    }
  }
  if (all(color_by %in% names(colData(x)))) {
    vals <- as.data.frame(colData(x), stringsAsFactors = FALSE)[color_by]
  }
  else if (all(color_by %in% rownames(x))) {
    if (!is.character(assay)) {
      stop("`assay` is not a character string")
    }
    if (!assay %in% assayNames(x)) {
      stop("input SpatialExperiment object does not have a `", 
           assay, "` assay")
    }
    vals <- assay(x, i = assay)
    vals <- as.data.frame(t(vals), stringsAsFactors = FALSE)[color_by]
  }
  else {
    stop("not all `color_by` values are present in the object")
  }
  coords <- reducedDim(x, type = dr)
  coords <- coords[, 1:2]
  colnames(coords) <- paste0(dr, 1:2)
  coords <- as.data.frame(coords)
  coords <- cbind(coords, vals[rownames(coords), , drop = FALSE])
  for (val in colnames(vals)) {
    p <- plot_scatter(data = coords, x = names(coords)[1], 
                      y = names(coords)[2], color_by = val, smooth = smooth, 
                      range = range, title = val)
    if (is.null(out_dir)) {
      if(interactive){
        require(plotly)
        return(htmlwidgets::createWidget(ggplotly(p)))
      }else{
        return(p)
      }
      
    }
    else {
      out_base <- glue("{out_dir}/{val}-{dr}")
      if (smooth) {
        out_base <- glue("{out_base}-smooth")
      }
      message(glue("generating {dr} plot for {val}"))
      
      if(interactive){
        require(plotly)
        
        htmlwidgets::saveWidget(ggplotly(p),glue('{out_base}.html'))
        
      }else{
        ggsave(filename = ifelse(pdf,glue("{out_base}.pdf"),glue("{out_base}.png")), plot = p, 
               width = ifelse(is.null(w),8,w), height = ifelse(is.null(h),5,h))
      }
      
    }
  }
  return(p)
}
plot_spatial.mod=function (x, color_by, assay = "logcounts", smooth = FALSE, range = c(0.01, 
                                                                                       0.99), out_dir = NULL,h=NULL,w=NULL,pdf=F,colors=NULL,interactive=F) 
{
  if (!is(x, "SpatialExperiment")) {
    stop("input is not a SpatialExperiment object")
  }
  if (!is.null(out_dir)) {
    if (!dir.exists(out_dir)) {
      stop("output directory `", out_dir, "` does not exist")
    }
  }
  if (all(color_by %in% names(colData(x)))) {
    vals <- as.data.frame(colData(x)[color_by], stringsAsFactors = FALSE)
  }
  else if (all(color_by %in% rownames(x))) {
    if (!is.character(assay)) {
      stop("`assay` is not a character string")
    }
    if (!assay %in% assayNames(x)) {
      stop("input SpatialExperiment object does not have a `", 
           assay, "` assay")
    }
    vals <- assay(x, i = assay)
    vals <- t(vals)[, color_by, drop = FALSE]
    vals <- as.data.frame(vals, stringsAsFactors = FALSE)
  }
  else {
    stop("not all `color_by` values are present in the object")
  }
  # if (!is.null(out_dir)) {
  #   dir.create(out_dir, showWarnings = FALSE)
  # }
  coords <- spatialCoords(x)
  coords <- as.data.frame(coords)
  coords$y <- coords$y * -1
  coords <- cbind(coords, vals[rownames(coords), , drop = FALSE])
  ratio <- max(coords$y)/max(coords$x)
  for (val in colnames(vals)) {
    if(is.null(colors)){
      p <- plot_scatter(data = coords, x = "x", y = "y", color_by = val, 
                        smooth = smooth, range = range, title = val, aspect_ratio = ratio)
    }else{
      p <- plot_scatter(data = coords, x = "x", y = "y", color_by = val, 
                        smooth = smooth, range = range, title = val, aspect_ratio = ratio)+scale_color_manual(values=colors)
    }
    
    if (is.null(out_dir)) {
      if(interactive){
        require(plotly)
        htmlwidgets::createWidget(ggplotly(p))
      }else{
        return(p)
      }
      
    }
    else {
      out_base <- glue("{out_dir}/{val}-spatial")
      if (smooth) {
        out_base <- glue("{out_base}-smooth")
      }
      message(glue("generating spatial plot for {val}"))
      if(interactive){
        htmlwidgets::saveWidget(ggplotly(p),glue('{out_base}.html'))
      }else{
        ggsave(filename = ifelse(pdf,glue("{out_base}.pdf"),glue("{out_base}.png")), plot = p, 
               width = ifelse(is.null(w),8,w), height = ifelse(is.null(h),5,h))
      }
      
    }
  }
}
assign_celltype_with_template=function(obj,phenotyping_template,cluster,mclust=F){

  g=cluster
  multiply_vectors <- function(boolean_gates) {
    Reduce(`*`, boolean_gates)
  }
  assay='logcounts'
  e <- scuttle::summarizeAssayByGroup(obj, ids = colData(obj)[[g]], assay.type = assay, statistics = "median")
  e <- SummarizedExperiment::assay(e, i = "median")
  e <- scale(t(e))

  if(mclust){
    require(mclust)
    scaled_ = assay(obj,'exprs')
    mclust.cols=apply(e,2,function(x){Mclust(x,G=2)})
    # find higher mean group
    group_means=lapply(mclust.cols,function(x) x$parameters$mean)

  }


  markers=lapply(seq(nrow(phenotyping_template)),function(x) unlist(strsplit(phenotyping_template$MARKERS[x],'[,]')))

  message('generating marker gates')
  marker_gates=do.call('cbind',lapply(markers,function(x){
    tmp=x
    boolean_gates=lapply(tmp,function(y){

      # pull last character (sign):
      sign.tmp=substr(y,nchar(y),nchar(y))
      # strip sign, leave only marker:
      marker.tmp=unlist(strsplit(y,'[+]|[-]'))
      # ensure marker has no numbers or letters after: (e.g. CD3 --> CD3 not CD31)
      pattern <- glue("{marker.tmp}(?![a-zA-Z0-9])")
      if(sign.tmp == '+'){
        if(!mclust){
          return(e[,grep(pattern,colnames(e),perl=T)] >0)
        }else{
          marker.idx=grep(pattern,colnames(e),perl=T)
          return(mclust.cols[[marker.idx]]$classification == names(group_means[[marker.idx]])[group_means[[marker.idx]]==max(group_means[[marker.idx]])])
        }

      }else{
        if(!mclust){
          return(e[,grep(pattern,colnames(e),perl=T)] <0)
        }else{
          marker.idx=grep(pattern,colnames(e),perl=T)
          return(mclust.cols[[marker.idx]]$classification == names(group_means[[marker.idx]])[group_means[[marker.idx]]==min(group_means[[marker.idx]])])
        }

      }

    })

    # entry-wise multiplication of boolean decision vectors:
    decision_vector=multiply_vectors(boolean_gates = boolean_gates)
    return(decision_vector)

  }))
  colnames(marker_gates)=phenotyping_template$CELLTYPE
  message('assigning celltypes')
  print(head(marker_gates))
  assign_celltype=sapply(seq(nrow(marker_gates)),function(x){
    tmp=colnames(marker_gates)[marker_gates[x,]==1]
    if(length(tmp) == 0){
      return('Unannotated')
    }else{
      return(paste0(tmp,collapse = ','))
    }
  })
  mapping=cbind(row.names(marker_gates),assign_celltype)
  print(mapping)
  assigned_celltype=as.character(obj[[g]])
  for(i in seq_along(mapping[,1])){
    assigned_celltype[assigned_celltype==mapping[i,1]]=mapping[i,2]
  }
  obj[[glue('{g}_annotations_template')]]=assigned_celltype
  return(list(mapping = mapping, marker_gates = marker_gates))
  #return(obj)

}
