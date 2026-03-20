generateHeatMap_VolcanoPlot=function( celltype_selection,spatial_obj,out_dir){
  
  require(ggplot2)
  require(ComplexHeatmap)
  require(dplyr)
  require(glue)
  
  print(celltype_selection)
  # KS test p-values:
  pvals_path=glue("{out_dir}/{list.files(out_dir,pattern='KS')}")
  pvals = read.csv(pvals_path)
  print(pvals)
  pvals= pvals[,-1]
  row.names(pvals)=names(pvals)
  print(pvals)
  
  celltype_selection = gsub(' ','.',celltype_selection)
  celltype_selection = gsub('-','.',celltype_selection)
  celltype_selection = gsub('[+]','pos',celltype_selection)
  celltype_selection = gsub('pos','.',celltype_selection)
  print(celltype_selection)
  pvals = pvals[match(celltype_selection,row.names(pvals)),match(celltype_selection,names(pvals))]
  
  # spatial_obj$celltype=gsub('[ ]','.',spatial_obj$celltype)
  # spatial_obj$celltype=gsub('-','.',spatial_obj$celltype)
  # spatial_obj$celltype=gsub('[+]','pos',spatial_obj$celltype)
  
  
  logOdds_path=glue("{out_dir}/{list.files(out_dir,pattern='logOdds')}")
  spatialIX.raw=read.csv(logOdds_path)
  
  spatialIX.raw=spatialIX.raw[,-1]
  row.names(spatialIX.raw)=names(spatialIX.raw)
  spatialIX.raw=spatialIX.raw[match(celltype_selection,names(spatialIX.raw)),match(celltype_selection,names(spatialIX.raw))]
  print(spatialIX.raw)
  
  print('spatial celltypes: ')
  print(colnames(spatialIX.raw))
  print('pval celltypes: ')
  print(colnames(pvals))
  
  ## pvalues:
  # outgoing:
  colnames(pvals)=rownames(pvals)
  print(pvals)
  transformed_pvalues.o=-log10(as.matrix(pvals))
  upper_bound=quantile(transformed_pvalues.o[!is.infinite(transformed_pvalues.o)],.99)
  if(sum(is.infinite(transformed_pvalues.o)) > 0){
    transformed_pvalues.o[is.infinite(transformed_pvalues.o)]=upper_bound
  }
  
  h=ComplexHeatmap::Heatmap(matrix = transformed_pvalues.o,name = '-log10 pvalues',cluster_rows = F,cluster_columns = F)
  pdf(glue('{out_dir}/interaction_pvalues-Heatmap.pdf'))
  print(h)
  dev.off()
  
  
  colnames(spatialIX.raw)=gsub('[/]','-',colnames(spatialIX.raw))
  rownames(spatialIX.raw)=gsub('[/]','-',rownames(spatialIX.raw))
  # volcano:
  volcano_dir=glue('{out_dir}/volcano-plots')
  dir.create(volcano_dir,showWarnings = F)
  setwd(volcano_dir)
  
  
  transformed_pvalues.df=data.frame(transformed_pvalues.o)
  spatialIX.raw.df=data.frame(spatialIX.raw)
  cFull=cbind(transformed_pvalues.df,spatialIX.raw.df)
  for(i in seq(dim(spatialIX.raw)[2])){
    c.subset=cFull[,c(i,i+dim(spatialIX.raw)[2])]
    names(c.subset)=c('pvalues','logOdds')
    c.subset$clusters=row.names(c.subset)
    p=ggplot(c.subset,aes(x=pvalues,y=logOdds,label=clusters))+geom_point()+ggrepel::geom_label_repel()+geom_hline(yintercept = 0)+xlab('-log10 pval')
    ggsave(glue('volcanoPlot-{colnames(spatialIX.raw)[i]}.pdf'),p)
  }
}
