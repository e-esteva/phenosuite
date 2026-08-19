getGPT_models=function(major=NULL,pricing=F,cache=NULL){
  require(curl)
  require(glue)
  require(stringr)
  require(stringi)
  
  if(is.null(cache)){
  	tmp=tempfile()
  	curl::curl_download('https://openai.com/pricing',tmp)
  	src_code =  suppressWarnings(readLines(tmp))
  	gpt_entry=grep('gpt',src_code)
  
  	span_delim=unlist(strsplit(src_code[gpt_entry],'[<span>][</span>]'))
  	t0=span_delim[grep('>gpt-',span_delim)]
  	gpt_models=gsub(pattern = "^[^gpt]*gpt",replacement = "",x = t0)
  	gpt_models_ = gpt_models[grep('^-',gpt_models)]
  	gpt_models_clean = unique(paste('gpt',gpt_models_,sep = ""))
  	# hanging '-i' --> '-instruct'
  	gpt_models_clean=gsub(pattern = '-i$',replacement = '-instruct',gpt_models_clean)
  	if(!pricing){
    		if(!is.null(major)){
      			gpt_models_clean=str_subset(gpt_models_clean,glue('gpt-{major}'))
    		}
    	return(gpt_models_clean)
  	}else{
    		t0=span_delim[grep('[$]',span_delim)]
    		models_pricing=lapply(seq(length(t0)),function(x) span_delim[grep('[$]',span_delim)[x]-seq(11)][11])
    		model_hits=unique(unlist(models_pricing))
    		model_hits.gpt=model_hits[grep('gpt',model_hits)]
    		gpt_hits.gpt.clean=paste('gpt',gsub(pattern = "^[^gpt]*gpt",replacement = "",x = model_hits.gpt),sep="")
    
    		pricing=t0[grep('s=\"f-body-1\">',t0)]
    		pricing.clean=paste('0',gsub(pattern = "^[^0.]*[0.]",replacement = "",x = pricing),sep = "")
    		pricing.clean.inputs=pricing.clean[seq(1,length(pricing.clean),2)]
    		pricing.clean.outputs=pricing.clean[seq(2,length(pricing.clean),2)]
    		# only first 12 rows correspond to the 6 gpts (make this selection a function of the number of gpts)
    		gpt_size=length(model_hits.gpt)
    		# pricing concatenated:
    		pricing_df=data.frame(cbind(pricing.clean.inputs,pricing.clean.outputs))
    		names(pricing_df)=c('Input','Output')
    		# each GPT has duplicate entries
    		pricing_df = pricing_df[seq(1,dim(pricing_df)[1],3),]
    		pricing_df.gpts = pricing_df[seq(gpt_size),]
    		pricing_df.gpts$Model=gpt_hits.gpt.clean
    
  	}
  }else{
	date=Sys.Date()
	gpts=list.files(cache)
	cache.dates=do.call('rbind',strsplit(gpts,'[_]'))[,1]
	if(date-max(as.Date(cache.dates)) > 90){
		tmp=tempfile()
        	curl::curl_download('https://openai.com/pricing',tmp)
        	src_code =  suppressWarnings(readLines(tmp))
        	gpt_entry=grep('gpt',src_code)

        	span_delim=unlist(strsplit(src_code[gpt_entry],'[<span>][</span>]'))
        	t0=span_delim[grep('>gpt-',span_delim)]
        	gpt_models=gsub(pattern = "^[^gpt]*gpt",replacement = "",x = t0)
       		gpt_models_ = gpt_models[grep('^-',gpt_models)]
        	gpt_models_clean = unique(paste('gpt',gpt_models_,sep = ""))
        	# hanging '-i' --> '-instruct'
        	gpt_models_clean=gsub(pattern = '-i$',replacement = '-instruct',gpt_models_clean)
        	if(!pricing){
                	if(!is.null(major)){
                        	gpt_models_clean=str_subset(gpt_models_clean,glue('gpt-{major}'))
                	}
        		return(gpt_models_clean)
        	}else{
                	t0=span_delim[grep('[$]',span_delim)]
                	models_pricing=lapply(seq(length(t0)),function(x) span_delim[grep('[$]',span_delim)[x]-seq(11)][11])
                	model_hits=unique(unlist(models_pricing))
                	model_hits.gpt=model_hits[grep('gpt',model_hits)]
                	gpt_hits.gpt.clean=paste('gpt',gsub(pattern = "^[^gpt]*gpt",replacement = "",x = model_hits.gpt),sep="")

                	pricing=t0[grep('s=\"f-body-1\">',t0)]
                	pricing.clean=paste('0',gsub(pattern = "^[^0.]*[0.]",replacement = "",x = pricing),sep = "")
                	pricing.clean.inputs=pricing.clean[seq(1,length(pricing.clean),2)]
                	pricing.clean.outputs=pricing.clean[seq(2,length(pricing.clean),2)]
                	# only first 12 rows correspond to the 6 gpts (make this selection a function of the number of gpts)
                	gpt_size=length(model_hits.gpt)
                	# pricing concatenated:
                	pricing_df=data.frame(cbind(pricing.clean.inputs,pricing.clean.outputs))
                	names(pricing_df)=c('Input','Output')
                	# each GPT has duplicate entries
                	pricing_df = pricing_df[seq(1,dim(pricing_df)[1],3),]
                	pricing_df.gpts = pricing_df[seq(gpt_size),]
                	pricing_df.gpts$Model=gpt_hits.gpt.clean

        	}
		write.csv(gpt_models_clean,glue('{cache}/{Sys.Date()}_GPT_models.csv'))
	}else{
		gpts_latest=read.csv(gpts,header=NULL)
		return(gpts_latest$V1)
	}
  }
}
