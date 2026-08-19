genHM_v3=function(celltype_selection,out_dir){
  
  require(ggplot2)
  require(ComplexHeatmap)
  require(dplyr)
  require(glue)
  
  
  effect_sizes=read.csv(glue("{out_dir}/{list.files(out_dir,pattern='KS-effect_sizes')}"),row.names = 1)
  colnames(effect_sizes)=row.names(effect_sizes)
  
  effect_sizes=effect_sizes[celltype_selection,celltype_selection]
  
  h=ComplexHeatmap::Heatmap(matrix = effect_sizes,name = 'Effect Sizes',cluster_rows = F,cluster_columns = F)
  pdf(glue('{out_dir}/Downsampled-Kolmogorov-Smirnov-Effect_size-Heatmap.pdf'))
  print(h)
  dev.off()
  
  # incidence tables:
  incidence_table=read.csv(glue("{out_dir}/{list.files(out_dir,pattern='probabilities_matrix')}"),row.names = 1)
  colnames(incidence_table)=row.names(incidence_table)
  
  incidence_table=incidence_table[celltype_selection,celltype_selection]
  
  h=ComplexHeatmap::Heatmap(matrix = incidence_table,name = 'Probability',cluster_rows = F,cluster_columns = F)
  pdf(glue('{out_dir}/Neighborhood-Incidence-Probability-Heatmap.pdf'))
  print(h)
  dev.off()
  
  
  
  
}