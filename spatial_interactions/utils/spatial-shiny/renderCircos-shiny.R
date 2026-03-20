renderCircos=function(logOdds,label,p1,p2,out_dir){
  require(colorspace)
  require(ComplexHeatmap)
  require(glue)
  require(circlize)
  require(colorRamp2)
  
  logOdds=as.matrix(logOdds)
  diag(logOdds)=0
  if(sum(is.infinite(logOdds))>0){
    logOdds=log(exp(logOdds)+1)
  }
  col_fun = colorRamp2(c(min(logOdds[,abs(colSums(logOdds)) >0]),median(logOdds[,abs(colSums(logOdds)) >0]),max(logOdds[,abs(colSums(logOdds)) >0])), c( "white",'yellow' ,"red"))
  par(cex = 0.9,mar = c(1.25, 1.25, 1.25, 1.25))
  chordDiagram(logOdds, grid.col = seq(dim(logOdds)[2]), symmetric = F, col = col_fun,directional = 1)
  title(glue(label,' | ',p1,';',p2))

  lgd_links = Legend(at = round(as.vector(quantile(logOdds[,abs(colSums(logOdds)) >0])),4), col_fun = col_fun, 
                   title_position = "topleft", title = "Links")

  lgd_list_vertical = packLegend(lgd_links)
  draw(lgd_list_vertical, x = unit(4, "mm"), y = unit(4, "mm"), just = c("left", "bottom"))
  
  pdf(glue('{out_dir}/{label}-circos.pdf'))
  logOdds=as.matrix(logOdds)
  diag(logOdds)=0
  
  col_fun = colorRamp2(c(min(logOdds[,abs(colSums(logOdds)) >0]),median(logOdds[,abs(colSums(logOdds)) >0]),max(logOdds[,abs(colSums(logOdds)) >0])), c( "white",'yellow' ,"red"))
  par(cex = 0.9,mar = c(1.25, 1.25, 1.25, 1.25))
  chordDiagram(logOdds, grid.col = seq(dim(logOdds)[2]), symmetric = F, col = col_fun,directional = 1)
  title(glue(label,' | ',p1,';',p2))
  
  
  # circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  #   xlim = get.cell.meta.data("xlim")
  #   ylim = get.cell.meta.data("ylim")
  #   sector.name = get.cell.meta.data("sector.index")
  #   circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  #   }, bg.border = NA)
  # 
  
  
  lgd_links = Legend(at = round(as.vector(quantile(logOdds[,abs(colSums(logOdds)) >0])),4), col_fun = col_fun, 
                     title_position = "topleft", title = "Links")
  
  lgd_list_vertical = packLegend(lgd_links)
  draw(lgd_list_vertical, x = unit(4, "mm"), y = unit(4, "mm"), just = c("left", "bottom"))
  dev.off()
}