generateHeatMap_VolcanoPlot.v2=function( celltype_selection,out_dir){
  
  require(ggplot2)
  require(ComplexHeatmap)
  require(dplyr)
  require(glue)
  
  print(celltype_selection)
  # hypergeometric test p-values:
  counts_path=glue("{out_dir}/{list.files(out_dir,pattern='phyper-counts')}")
  counts = read.csv(counts_path)
  print(counts)
  counts= counts[,-1]
  row.names(counts)=names(counts)
  print(counts)
  
  celltype_selection = gsub(' ','.',celltype_selection)
  celltype_selection = gsub('-','.',celltype_selection)
  celltype_selection = gsub('[+]','pos',celltype_selection)
  print(celltype_selection)
  counts = counts[match(celltype_selection,row.names(counts)),match(celltype_selection,names(counts))]
  
  celltype.table=read.csv(glue("{out_dir}/{list.files(out_dir,pattern='celltype_table')}"),row.names = 1)
  
  row.names(celltype.table)=gsub('[ ]','.',row.names(celltype.table))
  row.names(celltype.table)=gsub('-','.',row.names(celltype.table))
  row.names(celltype.table)=gsub('[+]','pos',row.names(celltype.table))
  ## computing hypergeometric test:
  tabl_ = celltype.table
  print(tabl_)
  pvalues_o = do.call('cbind',lapply(seq(dim(counts)[2]),function(z){
    c.tmp=sapply(row.names(counts),function(x) phyper(counts[match(x,row.names(counts)),z]-1,tabl_[x,],sum(celltype.table$celltype)-tabl_[x,],tabl_[row.names(counts)[z],],lower.tail=F) )
    return(c.tmp)
  }))
  
  logOdds_path=glue("{out_dir}/{list.files(out_dir,pattern='logOdds')}")
  spatialIX.raw=read.csv(logOdds_path)
  print(spatialIX.raw)
  spatialIX.raw=spatialIX.raw[,-1]
  row.names(spatialIX.raw)=names(spatialIX.raw)
  spatialIX.raw=spatialIX.raw[match(celltype_selection,names(spatialIX.raw)),match(celltype_selection,names(spatialIX.raw))]
  
  ## pvalues:
  # outgoing:
  colnames(pvalues_o)=rownames(pvalues_o)
  print(pvalues_o)
  transformed_pvalues.o=-log10(pvalues_o)
  upper_bound=quantile(transformed_pvalues.o[!is.infinite(transformed_pvalues.o)],.99)
  transformed_pvalues.o[is.infinite(transformed_pvalues.o)]=upper_bound
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
