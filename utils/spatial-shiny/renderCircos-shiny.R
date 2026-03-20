renderCircos=function(logOdds,label,p1,p2,out_dir,continuous_color_scheme=T,scale=F,discontinuity=F,col.fun=NULL,label_size.cex=0.9){
  require(colorspace)
  require(ComplexHeatmap)
  require(glue)
  require(circlize)
  require(colorRamp2)
  
  logOdds=as.matrix(logOdds)
  message(class(logOdds))
  message(colnames(logOdds))
  message(row.names(logOdds))

  diag(logOdds)=0
  if(sum(is.infinite(logOdds))>0){
    
    discontinuity=T

   }
   if(discontinuity){
    if(continuous_color_scheme){
	    if(is.null(col.fun)){
		col_fun = colorRamp2(c(min(logOdds[,abs(colSums(logOdds)) >0]),median(logOdds[,abs(colSums(logOdds)) >0]),max(logOdds[,abs(colSums(logOdds)) >0])), c( "white",'yellow' ,"red"))
	}else{
		col_fun=col.fun
	     }
    	
    }else{
      if(is.null(col.fun)){
	 col_fun = colorRamp2(c(min(logOdds[,abs(colSums(logOdds)) >0]),median(logOdds[,abs(colSums(logOdds)) >0]),max(logOdds[,abs(colSums(logOdds)) >0])), c( "blue",'white' ,"red"))
    	}else{
		col_fun =col.fun
	}
    }
  }else{
  	if(continuous_color_scheme){
		if(is.null(col.fun)){

    			col_fun = colorRamp2(c(min(logOdds[,abs(colSums(logOdds)) >0]),0,max(logOdds[,abs(colSums(logOdds)) >0])), c( "white",'yellow' ,"red"))
		}else{
			col_fun = col.fun
		}
  	}else{
		if(is.null(col.fun)){
    			col_fun = colorRamp2(c(min(logOdds[,abs(colSums(logOdds)) >0]),0,max(logOdds[,abs(colSums(logOdds)) >0])), c( "blue",'white' ,"red"))
		}else{
			col_fun = col.fun
		}
  	}
  }
  #par(cex = label_size.cex,mar = c(1.25, 1.25, 1.25, 1.25))
  #chordDiagram(logOdds, grid.col = seq(dim(logOdds)[2]), symmetric = F, col = col_fun,directional = 1,scale=scale)
  #title(glue(label,' | ',p1,';',p2))

  #lgd_links = Legend(at = round(as.vector(quantile(logOdds[,abs(colSums(logOdds)) >0])),4), col_fun = col_fun, 
  #                 title_position = "topleft", title = "Links")

  #lgd_list_vertical = packLegend(lgd_links)
  #draw(lgd_list_vertical, x = unit(4, "mm"), y = unit(4, "mm"), just = c("left", "bottom"))

  par(cex=label_size.cex,mar = c(1.25, 1.25, 1.25, 1.25))

  chordDiagram(logOdds, annotationTrack = "grid", preAllocateTracks = list(track.height = 0.1),scale=T,col = col_fun,grid.col = seq(dim(logOdds)[2]))
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
                xlim = get.cell.meta.data("xlim")
                xplot = get.cell.meta.data("xplot")
                ylim = get.cell.meta.data("ylim")
                sector.name = get.cell.meta.data("sector.index")
                if(abs(xplot[2] - xplot[1]) < 20) {
                        circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise",
                                niceFacing = TRUE, adj = c(0, 0.5))
                } else {
                        circos.text(mean(xlim), ylim[1], sector.name, facing = "inside",
                                niceFacing = TRUE, adj = c(0.5, 0))
                }
                }, bg.border = NA)
   lgd_links = Legend(at = round(as.vector(quantile(logOdds[,abs(colSums(logOdds)) >0])),4), col_fun = col_fun,
                   title_position = "topleft", title = "log-odds")

   lgd_list_vertical = packLegend(lgd_links)

   draw(lgd_list_vertical, x = unit(4, "mm"), y = unit(4, "mm"), just = c("left", "bottom"))
  
  if(!is.null(out_dir)){

  	pdf(glue('{out_dir}/{label}-circos.pdf'))
  
  	#par(cex = label_size.cex,mar = c(1.25, 1.25, 1.25, 1.25))
  	#chordDiagram(logOdds, grid.col = seq(dim(logOdds)[2]), symmetric = F, col = col_fun,directional = 1,scale=scale)
  	#title(glue(label,' | ',p1,';',p2))
  
  
  
  
  	#lgd_links = Legend(at = round(as.vector(quantile(logOdds[,abs(colSums(logOdds)) >0])),4), col_fun = col_fun, 
        #             title_position = "topleft", title = "Links")
  
  	#lgd_list_vertical = packLegend(lgd_links)
  	#draw(lgd_list_vertical, x = unit(4, "mm"), y = unit(4, "mm"), just = c("left", "bottom"))

	par(cex=label_size.cex,mar = c(1.25, 1.25, 1.25, 1.25))

	chordDiagram(logOdds, annotationTrack = "grid", preAllocateTracks = list(track.height = 0.1),scale=T,col = col_fun,grid.col = seq(dim(logOdds)[2]))
	title(glue(label,' | ',p1,';',p2))
	circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    		xlim = get.cell.meta.data("xlim")
    		xplot = get.cell.meta.data("xplot")
    		ylim = get.cell.meta.data("ylim")
    		sector.name = get.cell.meta.data("sector.index")
    		if(abs(xplot[2] - xplot[1]) < 20) {
        		circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise",
                    		niceFacing = TRUE, adj = c(0, 0.5))
    		} else {
        		circos.text(mean(xlim), ylim[1], sector.name, facing = "inside",
                    		niceFacing = TRUE, adj = c(0.5, 0))
    		}
		}, bg.border = NA)
	lgd_links = Legend(at = round(as.vector(quantile(logOdds[,abs(colSums(logOdds)) >0])),4), col_fun = col_fun,
                   title_position = "topleft", title = "log-odds")

	lgd_list_vertical = packLegend(lgd_links)
	draw(lgd_list_vertical, x = unit(4, "mm"), y = unit(4, "mm"), just = c("left", "bottom"))
  	dev.off()
   }
}
