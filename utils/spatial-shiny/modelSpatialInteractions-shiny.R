ModelSpatialInteractions.shiny=function(spatial_obj,cluster_col=NULL,downsample_num=NULL,resolution=0.377,p1=3,p2=30,out_dir=NULL,label,return_matrix=F,spe=NULL,draw=F,tolerance=0.001,view=F){
  #######################################
  library(circlize)
  library(ComplexHeatmap)
  library(dplyr)
  library(reshape2)
  library(glue)
  
  if(is.null(out_dir)){out_dir=getwd()}
  dir.create(glue('{out_dir}/{label}'),showWarnings = F)
  out_dir=glue('{out_dir}/{label}')
  withProgress(message = 'Generating Circos Plot', value = 0,{
    #source('maintainingProportions.R')
    # load in phenomenalist object:
    # spe = readRDS('spe.rds')
    # check clustering columns:
    if(is.null(cluster_col)){
      test.arrange=spatial_obj
    }else{
      spe=spatial_obj
      test.arrange=data.frame(spatialCoords(spe))
      
      test.arrange$cluster=spe[[cluster_col]]
      
    }
    incProgress(1/4, detail = 'Pulling spatial coords')
    
    
    
    
    if(!is.null(downsample_num)){
      idx=downsample.maintain.prop(test.arrange$cluster,ds.x = downsample_num,tolerance=tolerance)
      test.arrange=test.arrange[idx,]
    }
    
    
    arranged.coords=test.arrange%>%arrange(cluster)
    incProgress(2/4, detail = 'computing global dist matrix ...')
    
    message('computing global dist matrix ...')
    arranged.dist=dist(arranged.coords[,c('x','y')])
    
    incProgress(3/4, detail = 'converting dist object to operable matrix ...')
    message('converting dist object to operable matrix ...')
    arranged.dist.mat=melt(as.matrix(arranged.dist),varnames=c('row','col'))
    
    rm(arranged.dist)
    
    test3=matrix(arranged.dist.mat$value,dim(test.arrange)[1],dim(test.arrange)[1])
    rm(arranged.dist.mat)
    
    message('binning nodes by cluster ...')
    cluster.bins=table(arranged.coords$cluster)
    cluster.ranges=as.vector(c(0,cumsum(cluster.bins)))
    print(cluster.ranges)  
    cluster.blocks=lapply(seq(2,length(cluster.ranges)),function(x) seq(cluster.ranges[x-1]+1,cluster.ranges[x]))
    
    
    # conversion from pixel to micron--will have to adjust by platform:
    ##physical.pixel=0.3774
    physical.pixel=resolution
    # minimum distance:
    #p1=3
    # maximum distance
    #p2=30
    p2=p2+p1
    
    p1.scaled=p1/physical.pixel
    p2.scaled=p2/physical.pixel
    
    # only change from v1 --> names should be names(cluster.bin)
    counts=matrix(0,length(cluster.blocks),length(cluster.blocks))
    colnames(counts)=names(cluster.bins)
    rownames(counts)=names(cluster.bins)
    logOdds=matrix(0,length(cluster.blocks),length(cluster.blocks))
    colnames(logOdds)=names(cluster.bins)
    rownames(logOdds)=names(cluster.bins)
    message('computing log-odds of interactions ...')
    incProgress(4/4, detail = 'computing log-odds of interactions ...')
    
    for(i in seq(length(cluster.blocks))){
      interaction.space.tmp=test3[,cluster.blocks[[i]]]
      p1.all.interactions=sapply(seq(dim(interaction.space.tmp)[2]),function(x) length(interaction.space.tmp[,x][interaction.space.tmp[,x] > p1.scaled &interaction.space.tmp[,x]< p2.scaled]))
      logOdds.tmp=NULL
      counts.tmp=NULL
      for(j in seq(length(cluster.blocks))){
        q.clust.tmp=interaction.space.tmp[cluster.blocks[[j]],]
        p1.p2.interactions=sapply(seq(dim(q.clust.tmp)[2]),function(x) length(q.clust.tmp[,x][q.clust.tmp[,x] > p1.scaled & q.clust.tmp[,x]< p2.scaled]))
        
        uniq.interactions=length(unique(unlist( lapply(seq(dim(q.clust.tmp)[2]),function(x) match(q.clust.tmp[,x][q.clust.tmp[,x] > p1.scaled & q.clust.tmp[,x]< p2.scaled],q.clust.tmp[,x]) )) ))
        #counts.tmp[j]=ifelse(length(p1.p2.interactions) > 0,sum(p1.p2.interactions),0)
        counts.tmp[j]=uniq.interactions
        
        if(length(p1.p2.interactions) > 0){
          elA=sum(p1.p2.interactions)/sum(p1.all.interactions)
        }else{
          elA=0
        }
        elB=cluster.bins[[j]]/sum(cluster.bins)
        logOdds.tmp[j]=log(elA/elB)
      }
      logOdds[,i]=logOdds.tmp
      counts[,i]=counts.tmp
    }
    logOdds.unadj=logOdds
    logOdds=log(exp(logOdds)+(1-min(exp(logOdds))))
    
    ## computing hypergeometric test:
    tabl_ = table(spatial_obj$cluster)
    pvalues_o = do.call('cbind',lapply(seq(dim(counts)[2]),function(z){
      c.tmp=sapply(row.names(counts),function(x) phyper(counts[match(x,row.names(counts)),z]-1,tabl_[x],dim(spatial_obj)[1]-tabl_[x],tabl_[row.names(counts)[z]],lower.tail=F) )
      return(c.tmp)
    }))
    pvalues_i = do.call('cbind',lapply(seq(dim(counts)[2]),function(z){
      c.tmp=sapply(row.names(counts),function(x) phyper(counts[z,match(x,colnames(counts))]-1,tabl_[colnames(counts)[z] ],dim(spatial_obj)[1]-tabl_[colnames(counts)[z] ],tabl_[x],lower.tail=F) )
      return(c.tmp)
    }))
    if(!is.null(out_dir)){
      metadata=list(logOdds.unadj,logOdds,counts,pvalues_i,pvalues_o)
      if(is.null(downsample_num)){
        # saveRDS(logOdds.unadj,glue('{out_dir}/{label}-logOdds.rds'))
        # saveRDS(logOdds,glue('{out_dir}/{label}-logOdds-adj.rds'))
        # saveRDS(counts,glue('{out_dir}/{label}-counts.rds'))
        # saveRDS(pvalues_o,glue('{out_dir}/{label}-pvalue-matrix-outgoing.rds'))
        # saveRDS(pvalues_i,glue('{out_dir}/{label}-pvalue-matrix-incoming.rds'))
        # 
        saveRDS(metadata,glue('{out_dir}/{label}-metadata.rds'))
        pdf(glue('{out_dir}/{label}_no-self-interactions-circos.pdf'))
        diag(logOdds)=0
        
        col_fun = colorRamp2(c(min(logOdds),median(logOdds),max(logOdds)), c( "white",'yellow' ,"red"))
        par(cex = 0.9,mar = c(1.25, 1.25, 1.25, 1.25))
        chordDiagram(logOdds, grid.col = seq(dim(logOdds)[2]), symmetric = F, col = col_fun,directional = 1)
        title(glue(label,' | ',p1,';',p2))
        
        lgd_links = Legend(at = round(as.vector(quantile(logOdds)),4), col_fun = col_fun,
                           title_position = "topleft", title = "Links")
        
        lgd_list_vertical = packLegend(lgd_links)
        draw(lgd_list_vertical, x = unit(4, "mm"), y = unit(4, "mm"), just = c("left", "bottom"))
        dev.off()
        
      }else{
        # saveRDS(logOdds.unadj,glue('{out_dir}/{label}-{downsample_num}-logOdds.rds'))
        # saveRDS(logOdds,glue('{out_dir}/{label}-{downsample_num}-logOdds-adj.rds'))
        # saveRDS(counts,glue('{out_dir}/{label}-{downsample_num}-counts.rds'))
        # saveRDS(pvalues_o,glue('{out_dir}/{label}-{downsample_num}-pvalue-matrix-outgoing.rds'))
        # saveRDS(pvalues_i,glue('{out_dir}/{label}-{downsample_num}-pvalue-matrix-incoming.rds'))
        saveRDS(metadata,glue('{out_dir}/{label}-{downsample_num}-metadata.rds'))
        
        pdf(glue('{out_dir}/{label}_{downsample_num}-no-self-interactions-circos.pdf'))
        diag(logOdds)=0
        
        col_fun = colorRamp2(c(min(logOdds),median(logOdds),max(logOdds)), c( "white",'yellow' ,"red"))
        par(cex = 0.9,mar = c(1.25, 1.25, 1.25, 1.25))
        chordDiagram(logOdds, grid.col = seq(dim(logOdds)[2]), symmetric = F, col = col_fun,directional = 1)
        title(glue(label,' | ',p1,';',p2))
        
        lgd_links = Legend(at = round(as.vector(quantile(logOdds)),4), col_fun = col_fun,
                           title_position = "topleft", title = "Links")
        
        lgd_list_vertical = packLegend(lgd_links)
        draw(lgd_list_vertical, x = unit(4, "mm"), y = unit(4, "mm"), just = c("left", "bottom"))
        dev.off()
      }
    }
     
    message('generating circos plots ...')
    
    
    if(draw){
      ### no self-interactions:
      diag(logOdds)=0
      
      col_fun = colorRamp2(c(min(logOdds),median(logOdds),max(logOdds)), c( "white",'yellow' ,"red"))
      if (view) {
        
        par(cex = 0.9,mar = c(1.25, 1.25, 1.25, 1.25))
        chordDiagram(logOdds, grid.col = seq(dim(logOdds)[2]), symmetric = F, col = col_fun,directional = 1)
        title(glue(label,' | ',p1,';',p2))
        
        lgd_links = Legend(at = round(as.vector(quantile(logOdds)),4), col_fun = col_fun, 
                           title_position = "topleft", title = "Links")
        
        lgd_list_vertical = packLegend(lgd_links)
        draw(lgd_list_vertical, x = unit(4, "mm"), y = unit(4, "mm"), just = c("left", "bottom"))
        
        	
      }
    }
    
  })
  

  if(return_matrix){
    return(list(logOdds.unadj,logOdds,counts,pvalues_i,pvalues_o))
  }
  
}
