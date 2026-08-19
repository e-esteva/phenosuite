plot_pcfs_R=function(ppc,out_dir=NULL,res=0.377,radius=30){
  
  require(gridExtra)
  require(MatrixGenerics)

  stepsize=(radius/res)/500
  
  sample.count=length(unique(ppc[[1]]$Patient))
  


  grobs=lapply(seq(length(ppc)),function(x){
    
    ppc_ = ppc[[x]] %>% data.frame() %>% select(c(Cell_one,Cell_two,PCF))
    
    pcfs_ = data.frame(do.call('cbind',lapply(ppc_$PCF,function(x) x)))
    
    if(x < length(ppc)){
      colnames(pcfs_)[seq(x*sample.count)]=ppc_$Cell_one[seq(x*sample.count)]
      colnames(pcfs_)[seq(x*sample.count+1,length(colnames(pcfs_)))]=ppc_$Cell_two[seq(x*sample.count+1,length(colnames(pcfs_)))]
      
    }else{
      colnames(pcfs_)=ppc_$Cell_one
    }
    
    
    celltypes=unique(colnames(pcfs_))
    
    if(sample.count > 1){
      require(tidyr)
      
      # compute average PCF over samples for each point:
      pcfs_sample_norm=data.frame(do.call('cbind',lapply(celltypes,function(x){
        means=rowMeans(pcfs_[,colnames(pcfs_)==x])
        return(means)
      })))
      
      names(pcfs_sample_norm)=celltypes
      
      # compute variance across samples for each point:
      pcfs_sample_norm.vars=data.frame(do.call('cbind',lapply(celltypes,function(x){
        vars=rowVars(as.matrix(pcfs_[,colnames(pcfs_)==x]))
        return(vars)
      })))
      names(pcfs_sample_norm.vars)=celltypes
      
      pcfs_ = pcfs_sample_norm
      
      pcfs_$Radius=seq(dim(pcfs_)[1])
      pcfs_sample_norm.vars$Radius=seq(dim(pcfs_)[1])
      
      
      mu.unicode='\u00B5'
      
      #### join mean and variance tables on Radius:
      means_long.pcf <- pivot_longer(pcfs_, -Radius, values_to = "mean", names_to = "variable")
      vars_long.pcf <- pivot_longer(pcfs_sample_norm.vars, -Radius, values_to = "var", names_to = "variable")
      
      df_join.pcf <- means_long.pcf %>% 
        left_join(vars_long.pcf)
      
      p=ggplot(data = df_join.pcf, aes(x = Radius*res*stepsize,y=value, group = variable)) + 
        geom_line(aes(y = mean, color = variable), size = 1) + 
        geom_ribbon(aes(y = mean, ymin = mean - var, ymax = mean + var, fill = variable), alpha = .2) +
        xlab(glue('Radius ({mu.unicode}m)'))+
        ylab('PCF')+
        ggtitle(colnames(pcfs_)[x])+
        theme_bw() +  
        theme(legend.key = element_blank()) + 
        geom_hline(yintercept = 1)+
        theme(legend.title = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ylim(0,2)
      
    }else{
      
      
      
      pcfs_$Radius=seq(dim(pcfs_)[1])
      
      pcfs_melted=melt(pcfs_,id='Radius')
      
      
      mu.unicode='\u00B5'
      p=ggplot(pcfs_melted,aes(x=Radius*res*stepsize,y=value,col=variable))+geom_line()+geom_hline(yintercept = 1)+xlab(glue('Radius ({mu.unicode}m)'))+ylab('PCF')+ggtitle(colnames(pcfs_)[x])+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ylim(0,2)
    }
    return(p)
  })
  if(!is.null(out_dir)){
    g=gridExtra::arrangeGrob(grobs = grobs)
    
    ggsave(glue('{out_dir}/PCF-plots.pdf'),g,w=20,h=20)
  }else{
    return(grobs)
  }
}
