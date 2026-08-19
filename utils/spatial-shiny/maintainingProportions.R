downsample.maintain.prop=function(data,ds.x,tolerance=0.001,seed=NULL){
  if(is.null(seed)){
    set.seed(99)
  }else{
    set.seed(seed)
  }
  
  ref.prop=as.numeric(as.character(table(data)))/sum(as.numeric(as.character(table(data))))
  ds.idx=sample(length(data),ds.x,replace = F)
  data.ds0=data[ds.idx]
  ds.prop=as.numeric(as.character(table(data.ds0)))/sum(as.numeric(as.character(table(data.ds0))))
  mae=mean(abs(ref.prop-ds.prop))
  if(mae<=tolerance){
    return(ds.idx)
  }else{
    correction.vector=round((ref.prop-ds.prop)*ds.x)
    keep.idx=list()
    for(i in seq(length(correction.vector))){
      if(correction.vector[i] < 0){
        # idx of element in downsampled
        var.tmp.idx=ds.idx[data.ds0==names(table(data))[i]]
        # randomly dropping excess cells from this set:
        random.drop=var.tmp.idx[-sample(length(var.tmp.idx),abs(correction.vector[i]),replace = F)]
        keep.idx[[i]]=random.drop
      }else{
        if(correction.vector[i] > 0){
          # idx of element in ref
          var.baseline.idx=seq(length(data))[data==names(table(data))[i]]
          # idx of element in downsampled
          var.tmp.idx=ds.idx[data.ds0==names(table(data))[i]]
          # removing downsampled idx from baseline
          data.tmp=var.baseline.idx[-var.tmp.idx]
          # drawing the necessary amount randomly from remaining idx
          data.draw=data.tmp[sample(length(data.tmp),abs(correction.vector[i]),replace = F)]
          keep.idx[[i]]=c(var.tmp.idx,data.draw)
        }else{
          var.tmp.idx=ds.idx[data.ds0==names(table(data))[i]]
          keep.idx[[i]]=var.tmp.idx
        }
        
      }
    }
    return(unlist(keep.idx))
  }
  
}
