renderCircos=function(logOdds,label,p1=NULL,p2=NULL,out_dir=NULL,col_fun=NULL,scale=F){
  require(colorspace)
  require(ComplexHeatmap)
  require(glue)
  require(circlize)
  require(colorRamp2)
  
  logOdds=as.matrix(logOdds)
  diag(logOdds)=0
  
  if(is.null(col_fun)){
    col_fun = colorRamp2(c(min(logOdds[,abs(colSums(logOdds)) >0]),0,max(logOdds[,abs(colSums(logOdds)) >0])), c( "white",'yellow' ,"red"))
  }else{
    col_fun = col_fun
  }
  
  par(cex = 0.9,mar = c(1.25, 1.25, 1.25, 1.25))
  chordDiagram(logOdds, grid.col = ggsci::pal_aaas()(dim(logOdds)[2]), symmetric = F, col = col_fun,directional = 1,scale = scale)
  if(!is.null(p1)){
    title(glue(label,' | ',p1,';',p2))
  }else{
    title(label)
  }
  

  lgd_links = Legend(at = round(as.vector(quantile(logOdds[,abs(colSums(logOdds)) >0])),4), col_fun = col_fun, 
                   title_position = "topleft", title = "Links")

  lgd_list_vertical = packLegend(lgd_links)
  draw(lgd_list_vertical, x = unit(4, "mm"), y = unit(4, "mm"), just = c("left", "bottom"))
  
  if(!is.null(out_dir)){
    pdf(glue('{out_dir}/{label}-circos.pdf'))
    logOdds=as.matrix(logOdds)
    diag(logOdds)=0
    
    #col_fun = colorRamp2(c(min(logOdds[,abs(colSums(logOdds)) >0]),median(logOdds[,abs(colSums(logOdds)) >0]),max(logOdds[,abs(colSums(logOdds)) >0])), c( "white",'yellow' ,"red"))
    par(cex = 0.9,mar = c(1.25, 1.25, 1.25, 1.25))
    chordDiagram(logOdds, grid.col = ggsci::pal_aaas()(dim(logOdds)[2]), symmetric = F, col = col_fun,directional = 1,scale=scale)
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
  
}
