def model_spatial_interactions(spatial_obj,out_dir,label,resolution,p1,p2,min_count):
  import os 
  import numpy as np
  import pandas as pd
  import scipy
  from scipy import spatial, io, sparse
  import re
  from itertools import chain
  
  

  os.makedirs(out_dir,exist_ok=True)
  
  print('outdir')
  print(out_dir)
  min_count = min_count
  
  
  dist_mat = scipy.spatial.distance.pdist(spatial_obj[['x','y']])
  dist_square = scipy.spatial.distance.squareform(dist_mat)

  # pull all class idx:
  celltypes=[re.sub(pattern='[+]',repl='pos',string=x) for x in spatial_obj['celltype'] ]
  classes = list(set(celltypes))



  class_idx = [[x for x in range(len(celltypes)) if len(re.findall(pattern='^'+str(y)+'$', string = celltypes[x])) > 0] for y in classes]
  
  class_sizes = np.array([len(x) for x in class_idx])

  resolution=resolution
  # minimum distance:
  p1=p1
  # maximum distance
  p2=p2
  p2=p2+p1
  
  p1_scaled=p1/resolution
  p2_scaled=p2/resolution

  phyper=[]
  logOdds=[]
  for i in range(len(class_idx)):
    # only consider celltypes with > min_count:
    if class_sizes[i] > min_count:
      interaction_space = dist_square[:,class_idx[i]]
      interaction_space_df = pd.DataFrame(interaction_space)
      print('computing '+str(classes[i])+' global interactions within '+str(p1)+':'+str(p2)+' microns ..')
      r_all_interactions = np.sum(np.array([len(interaction_space_df[(interaction_space_df[interaction_space_df.columns[x]] > p1_scaled) & (interaction_space_df[interaction_space_df.columns[x]] < p2_scaled)][x]) for x in range(len(interaction_space_df.columns))]))
  
      q_logOdds = []
      phyper_metadata = []
      for j in range(len(class_idx)):
        print('computing '+str(classes[i])+' | '+str(classes[j])+' interactions within '+str(p1)+':'+str(p2)+' microns ..')
      
        q_interaction_subset = interaction_space[class_idx[j],:]
        q_interaction_subset_df = pd.DataFrame(q_interaction_subset)
    
        r_q_interactions = np.sum(np.array([len(q_interaction_subset_df[(q_interaction_subset_df[q_interaction_subset_df.columns[x]] > p1_scaled) & (q_interaction_subset_df[q_interaction_subset_df.columns[x]] < p2_scaled)][x]) for x in range(len(q_interaction_subset_df.columns))]))
        r_q_interactions_idx = [(q_interaction_subset_df[(q_interaction_subset_df[q_interaction_subset_df.columns[x]] > p1_scaled) & (q_interaction_subset_df[q_interaction_subset_df.columns[x]] < p2_scaled)][x]).index for x in range(len(q_interaction_subset_df.columns))]
    
        elA = r_q_interactions / r_all_interactions
        elB = len(class_idx[j]) / len(celltypes) 
    
      
        q_logOdds.append(np.log(elA / elB))
      
        # unlist indices:
        # for hypergeometric test: need uniq 
        r_q_interactions_idx_ = list(chain.from_iterable(r_q_interactions_idx))
        r_q_interactions_uniq_idx_ = len(set(r_q_interactions_idx_))
        phyper_metadata.append(r_q_interactions_uniq_idx_)
      
      logOdds.append(q_logOdds)
      phyper.append(phyper_metadata)
    
  
  logOdds_df = pd.DataFrame(np.vstack(logOdds).T)
  logOdds_df.columns = classes
  logOdds_df.index = classes
  

  logOdds_df.to_csv(str(out_dir)+'/logOdds-matrix-'+str(label)+'.csv')
  
  phyper_counts=pd.DataFrame(phyper).T
  phyper_counts.columns = classes
  phyper_counts.index = classes
  phyper_counts.to_csv(str(out_dir)+'/phyper-counts-matrix-'+str(label)+'.csv')
  
  table_ = spatial_obj.celltype.value_counts()
  table_.to_csv(str(out_dir)+'/celltype_table.csv')
  
  return logOdds_df
