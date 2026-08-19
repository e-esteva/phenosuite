getHeatmapByCluster=function(obj,celltypes.df=NULL,genes,out=NULL,stat='median',assay='RNA',cluster_col=NULL){
  library(SingleCellExperiment)
  #library(iCellR)
  library(Seurat)
  library(ComplexHeatmap)
  library(cowplot)
  library(colorRamp2)
  library(glue)
  
  
  if(class(obj) == "iCellR"){
    x <- SingleCellExperiment(assay = list(counts = obj@main.data[genes[genes %in% row.names(obj@main.data)],]), colData = obj@best.clust$clusters)
    
    e <- scuttle::summarizeAssayByGroup(x, ids = obj@best.clust$clusters, 
                                        assay.type = 'counts', statistics = stat)
  }else{
    if(class(obj)=='Seurat'){
	assay_classes = lapply(obj@assays, class)
	if("Assay5" %in% assay_classes[[assay]]){
		message('assay5 protocol')
		data_ = LayerData(obj, assay = assay, layer = "data")
		x <- SingleCellExperiment(assay = list(counts = data_[genes[genes %in% row.names(data_)],]), colData = obj[[cluster_col]])
	}else{
		message('no assay5 assay detected')
		message(assay_classes[[assay]])
      		x <- SingleCellExperiment(assay = list(counts = obj[[assay]]@data[genes[genes %in% row.names(obj[[assay]]@data)],]), colData = obj[[cluster_col]])
	}
      message(glue('x dims: {dim(x)}'))
      message(glue('length of ids: {length(obj[[cluster_col]][[1]])}'))
      if(!(cluster_col %in% names(obj@meta.data))){
	message(glue('{cluster_col} not found in data'))
      }
      else{
      	message(cat(head(obj[[cluster_col]][[1]]),'\n'))
      }
      e <- scuttle::summarizeAssayByGroup(x, ids = obj[[cluster_col]][[1]], 
                                          assay.type = 'counts', statistics = stat)
    }else{
      stop('Input needs to be either iCellR or Seurat')
    }
  }
  
  
  e <- SummarizedExperiment::assay(e, i = stat)
  e <- scale(t(e))
  e[e > 4] <- 4
  e[e < -4] <- -4
  gradient_colors <- rev(RColorBrewer::brewer.pal(11, "RdYlBu"))
  e_ = do.call('cbind',apply(e,2,na.omit))
  #celltypes.df=data.frame('celltype'=c('HSC','ST-HSC/MPP3','Kras-specific MPP','MPP4','MyP','MEP','Proliferative MPP'))
  
  
  
  h=ComplexHeatmap::Heatmap(e_, name = "Expression", 
                          row_title = 'clusters', col = gradient_colors,cluster_rows = F,cluster_columns = F,row_split = factor(row.names(e_),levels=row.names(e_)))
  
  if(!is.null(out)){
    pdf(glue('{out}/Heatmap_by_cluster-{stat}.pdf'))
    h
    dev.off()
  }else{
    return(h)
  }
}

