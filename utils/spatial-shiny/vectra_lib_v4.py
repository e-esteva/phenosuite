import glob
import re
import os

import numpy as np
import pandas as pd
from scipy import stats, spatial

import matplotlib as mpl
import matplotlib.pyplot as plt



import pickle

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', 500)
pd.options.mode.chained_assignment = None  # default='warn'

#%matplotlib inline
mpl.rcParams['font.size'] = 14
import seaborn as sns






def read_csv_tsv(filename):
    file = pd.read_csv(filename, delimiter='\t') #try tsv
    if 'Path' not in file.columns: #try comma-delim
        file = pd.read_csv(filename)        
    #if 'Path' not in file.columns:
    #    raise Exception(filename+" did not open properly!")
        
    #detect_legacy_vectra
    if 'Sample_Name' in file.columns: #underscores bad
        replace_parens = lambda x: '(' + x.group(0) + ')'
        file.columns = (file.columns.str.replace('_', ' ')
                        .str.replace('Opal [\S]*', replace_parens)
                        .str.replace('Normalized Counts Total Weighting',
                                     '(Normalized Counts, Total Weighting)')
                        .str.replace('HLA DR', 'HLA-DR')
                       )
    return file


def extract_data(directory, classification = None, verbose = True, 
                 drop_nan = True, drop_duplicates = False, debug = False):
    """
    Extracts cell information from cell_seg_data files and pairs it with 
    corresponding score files from score_data.

    Args:
        directory: string
            string of parent directory containing all files
        classification: function
            must take a row of the dataframe and output a string 
        verbose: bool
            output all quality-control checking
        drop_nan: bool
            remove rows with NaN phenotype
        drop_duplicates: bool
            whether to remove duplicates with the same file name. Most recently
            modified will be kept.
        debug: bool
            if True, only compile score info and return score_files

    Returns:
        output: list of dicts
            each corresponding to an image file
            
        unless debug is True, in which case:
        score_files: pandas df
            each extracted score file for debugging purposes
    """
    
    #file_list = glob.glob(directory+'/**/*_cell_seg_data.csv', recursive=True)
    file_list = glob.glob(directory+'/**/*.csv', recursive=True)
    print('file_list: '+str(file_list))
    #pattern = re.compile('.*/([^/]*)_cell_seg_data.csv')
    pattern = re.compile('.*/([^/]*).csv')
    output = [{'File Path': name, 'File Name': pattern.search(name).group(1)} 
             for name in file_list]
    
    if drop_duplicates:
        files_df = pd.DataFrame(output)
        duplicated_files = files_df[
            files_df.duplicated('File Name', keep=False)
            ].copy()
        duplicated_files['Time'] = duplicated_files['File Path'].apply(
            os.path.getmtime)
        duplicated_files.sort_values('Time', inplace=True)
        
        to_remove = duplicated_files[duplicated_files
                                     .duplicated('File Name', keep='last')]
        output_filtered = [item for item in output 
                           if to_remove['File Path']
                           .isin([item['File Path']]).sum() == 0]
        output = output_filtered
        
    score_file_list = glob.glob(directory+'/**/*_score_data.txt', 
                                recursive=True)
    
    if len(score_file_list) == 0: #no score files
        if verbose:
            print("No score files detected!\n")
        for sample in output:
            file = read_csv_tsv(sample['File Path'])
            sample['Sample Name'] = file.loc[0,'Sample Name']

            df = file[['Cell X Position', 'Cell Y Position', 'Tissue Category',
                       'Phenotype']].copy()
            df.columns = ['Cell X Position', 'Cell Y Position', 'Tissue Category',
                          'Phenotype']

            if drop_nan:
                df = df[~df['Phenotype'].isna()]

            if classification is not None:
                df['Classification'] = df.apply(classification, axis=1)

            sample['Data'] = df
        return output

    # if score files are present:
    score_files = pd.concat([
        pd.concat([
            read_csv_tsv(score_file) if not drop_duplicates
            else read_csv_tsv(score_file).assign(
                Time = os.path.getmtime(score_file))
            for score_file in score_file_list 
            if sample['File Name'] in score_file
            ]) 
        for sample in output], sort=False)
    
    if drop_duplicates:
        score_files.sort_values('Time', inplace=True)
        score_files = score_files[~score_files.duplicated(
            subset=['Sample Name', 'First Stain Component', 
                    'Second Stain Component'], keep='last')]
        
    score_files.set_index('Sample Name', inplace=True)
    
    stain_properties = pd.concat(
        (score_files[[
            'First Cell Compartment', 
            'First Stain Component'
            ]].T.reset_index(drop=True),
         score_files[[
            'Second Cell Compartment', 
            'Second Stain Component'
             ]].T.reset_index(drop=True)
        ), axis=1).T
        
    stain_properties.columns = ['compartment', 'component']
    
    threshold_names = list(stain_properties['component'].drop_duplicates() 
                           + ' Threshold')

    score_verification = score_files.groupby('Sample Name')[
        threshold_names + ['Number of Cells']].std()>0
    
    if verbose:
        print('Stains and compartments:\n', 
              stain_properties.drop_duplicates().reset_index(drop=True),
             '\nDiscrepancies: ', 
              score_verification[score_verification.sum(axis=1)>0]
              if (score_verification.sum().sum()>0) else 'None'
             )
        
    if debug:
        return score_files
    
    thresholds = score_files.groupby('Sample Name').mean()
        
    for sample in output:
        file = read_csv_tsv(sample['File Path'])
        sample['Sample Name'] = file.loc[0,'Sample Name']
        
        stain_column_names = [
            row['compartment'] + ' ' + row['component'] 
            + ' Mean (Normalized Counts, Total Weighting)' 
            for idx, row in 
            stain_properties.loc[sample['Sample Name']].iterrows()
            ]
        stain_names = [name.split(' ')[1] for name in stain_column_names]

        df = file[['Cell X Position', 'Cell Y Position', 'Tissue Category',
                   'Phenotype'] + stain_column_names].copy()
        df.columns = ['Cell X Position', 'Cell Y Position', 'Tissue Category',
                      'Phenotype'] + stain_names 

        for stain in threshold_names:
            current_stain = stain.split(' ')[0]
            df[current_stain] = np.sign(df[current_stain].fillna(0) 
                                - thresholds.loc[sample['Sample Name'], stain] 
                                + 0.000001)  # round up to avoid 0
        
        if drop_nan:
            df = df[~df['Phenotype'].isna()]
        
        if classification is not None:
            df['Classification'] = df.apply(classification, axis=1)

        sample['Data'] = df
    return output

def count_cells(output, grouping = 'Phenotype', density = True):
    """
    Counts number or density of cells for each image.

    Args:
        output: list of dicts
            from extract_data
        grouping: 'Phenotype' or 'Classification' or list of variables
            whether to use VECTRA Phenotype or self-generated classification
        density: bool
            whether to density-normalize each image

    Returns:
        pandas df, with columns as image names and rows as cell types
    """
    
    counts_table = pd.concat(
        [sample['Data'].groupby(grouping).size().rename(sample['Sample Name']) 
         for sample in output], axis=1, sort=True).fillna(0)

    errors = counts_table.columns[counts_table.sum()==0]
    if len(errors) > 0:
        counts_table.drop(errors, axis=1, inplace=True)
        print("Warning: samples were omitted due to missing values: ", errors)
    
    if density:
        counts_table = (counts_table/counts_table.sum())       
        
    return counts_table

def pcf(output, cell_types = None, phenotype = 'Phenotype', 
        count_threshold = 1,radius=30,resolution=0.377):
    """
    Calculates pair correlation function using spatstat from R.
    Extra dependencies: rpy2 (3.2.2+), R (3.6+), spatstat (1.62-2+)

    Args:
        output: list of dicts
            from extract_data
        cell_types: list
            Names of pertinent cell types, which will be prepended by 'All'. 
            Default is all recognized phenotypes
        phenotype: string, 'Phenotype' or 'Classification'
            whether to use VECTRA Phenotype or self-generated classification
        count_threshold: int
            threshold to not calculate a specific pcf

    Returns:
        pandas df, with each row a sample
    """

    #only needs rpy2 if calculating pcf
    import rpy2
    import rpy2.robjects as robjects
    from rpy2.robjects.packages import importr
    from rpy2.robjects import pandas2ri
    
    pandas2ri.activate()
    r = robjects.r
    spatstat = importr("spatstat.geom")
    spatstat_ = importr("spatstat.explore")
    
    if cell_types is None:
        cell_types = list(pd.concat([sample['Data'][phenotype] 
                                    for sample in output]
                                    ).value_counts().index)       
    cell_types = ['All'] + cell_types
        
    pcf_output = []
    
    for sample in output:
        selection = sample['Data'] 
        patient_name = sample['File Name']
        sample_name = sample['Sample Name']

        for cell_one in cell_types:
            
            starting_index_two = cell_types.index(cell_one)
            for cell_two in cell_types[starting_index_two:]:            

                if ((sum(selection[phenotype] == cell_one) 
                    < count_threshold and cell_one != "All")
                    or (sum(selection[phenotype] == cell_two) 
                    < count_threshold and cell_two != "All")):
                    
                    cell_count = [len(selection[phenotype]) if x == 'All' 
                                  else sum(selection[phenotype] == x) 
                                  for x in [cell_one, cell_two]]
                    
                    pcf_output.append([patient_name, sample_name,cell_one,
                                       cell_two, "NA", "NA", min(cell_count), 
                                       cell_count[0], cell_count[1]])

                else:
                    markers = selection[['Tissue Category', phenotype]];
                    x_coords = np.round(
                        selection['Cell X Position']).astype(int);
                    y_coords = np.round(
                        selection['Cell Y Position']).astype(int);
                    point_cloud = spatstat.ppp(
                        x_coords, y_coords, 
                        robjects.IntVector(
                            [x_coords.min()-1, x_coords.max()+1]),
                        robjects.IntVector(
                            [y_coords.min()-1, y_coords.max()+1]),
                        marks=markers.astype(str));
                    point_cloud = r.cut(point_cloud, z=phenotype);
                    print('point-cloud diagnostics')
                    print(point_cloud)
                    celltypes=list(sample['Data']["Phenotype"])
                    marks_ =list(point_cloud.rx2('marks'))
                    mapping=dict(zip(celltypes,marks_))
                    #very close to default. constant for stability
                    #x_values = np.arange(0, 125, .25) 
                    rmax=radius/resolution
                    # coercing to 500 to approximate default of 513
                    step=rmax/500
                    x_values = np.arange(0,rmax,step)
                    
                    if cell_two == "All": #regular pcf
                        cell_count = r.c(point_cloud.rx2("n"), 
                                         point_cloud.rx2("n"));
                        lambda1 = r.density(point_cloud, at="points",
                                            leaveoneout=False);
                        lambdasums  = r.sum(1/lambda1)**2;
                        pcfplot = spatstat_.pcfinhom(
                            point_cloud, robjects.FloatVector(lambda1), 
                            correction = "isotropic", r = x_values, 
                            renormalise = False
                            );
                        
                        
                    elif cell_one == "All": #dot pcf
                        
                        #c2 = spatstat.subset_ppp(
                        #    point_cloud, point_cloud.rx2('marks') == mapping[cell_two]
                        #    );
                        c2 = spatstat.subset_ppp(
                            point_cloud, np.array(tuple(point_cloud.rx2('marks'))) == mapping[cell_two]
                        );

                        cell_count = r.c(point_cloud.rx2("n"), c2.rx2("n"));
                        lambda1 = r.density(point_cloud, at="points", 
                                            leaveoneout=False);
                        lambda2 = r.density(c2, at="points",
                                            leaveoneout=False);
                        lambdasums = r.sum(1/lambda1)*r.sum(1/lambda2);
                        pcfplot = spatstat_.pcfdot_inhom(
                            point_cloud, cell_two,
                            robjects.FloatVector(lambda2),
                            robjects.FloatVector(lambda1), r=x_values,
                            correction = "isotropic"
                            );
                        
                    else: #cross pcf
                        #c1 = spatstat.subset_ppp(
                        #    point_cloud, point_cloud.rx2('marks') == mapping[cell_one]
                        #    );
                        #c2 = spatstat.subset_ppp(
                        #    point_cloud, point_cloud.rx2('marks') == mapping[cell_two]
                        #    );
                        c1 = spatstat.subset_ppp(
                            point_cloud, np.array(tuple(point_cloud.rx2('marks'))) == mapping[cell_one]
                            );
                        c2 = spatstat.subset_ppp(
                            point_cloud, np.array(tuple(point_cloud.rx2('marks'))) == mapping[cell_two]
                            );

                        cell_count = r.c(c1.rx2("n"), c2.rx2("n"));
                        lambda1 = r.density(c1, at="points", 
                                            leaveoneout=False);
                        lambda2 = r.density(c2, at="points",
                                            leaveoneout=False);
                        lambdasums = r.sum(1/lambda1)*r.sum(1/lambda2);
                        pcfplot = spatstat_.pcfcross_inhom(
                            point_cloud, cell_one, cell_two, 
                            robjects.FloatVector(lambda1), 
                            robjects.FloatVector(lambda2),
                            r = x_values, correction = "isotropic"
                            );

                    #pcftable = pcfplot["iso"];
                    # has to be int index:
                    pcftable=pcfplot[list(pcfplot.names).index('iso')]
                    # first val is Inf so set to 0:
                    #pcftable.iloc[0] = 0
                    pcftable[0]=0
                    normalization = (lambdasums/(
                        (x_coords.max() - x_coords.min()) 
                        * (y_coords.max() - y_coords.min())
                        )**2)[0]
                    pcf_output.append([
                        patient_name, sample_name, cell_one, cell_two, 
                        pcftable, pcftable/normalization,
                        sum(pcftable/normalization),
                        normalization, min(cell_count), 
                        cell_count[0], cell_count[1]
                        ]) 

    pcf_output = pd.DataFrame(pcf_output, columns = [
        'Patient', 'Sample Name', 'Cell_one', 'Cell_two', 'PCF', 'normPCF',
        'PCFsum', 'normalization', 'min_count', 'count_one', 'count_two'
        ])
    
    return pcf_output

def pcf_subset(df, cell1, cell2, min_count=20, max_radius = 200):
    """
    Analysis of pcfs of two specified cell types

    Args:
        df: pandas dataframe
            pcf dataframe from vl.pcf
        cell1, cell2: str
            Names of cell types of interest (or 'All') 
        min_count: int or list
            exclude samples with fewer than this many cells of each type
        max_radius: int
            steps to sum pcf to. steps are approx. 0.25 microns

    Returns:
        pandas df, with each row a sample
    """
    cell_mask = (((df['Cell_one']==cell1) & (df['Cell_two']==cell2)) 
                 | ((df['Cell_one']==cell2) & (df['Cell_two']==cell1)))
    count_mask = ((df['min_count'] > min_count) 
                  if type(min_count) == int 
                  else ((df['count_one'] > min_count[0]) 
                        & (df['count_two'] > min_count[1])))
    pcf_table = df[cell_mask & count_mask].copy()
    pcf_table['normalization'] = pd.to_numeric(pcf_table['normalization'])
    pcf_table['PCF'] = pcf_table['PCF'].apply(
        lambda x: np.insert(x[1:],0,0)
        ) / pcf_table['normalization']
    pcf_table = pcf_table.assign(
        pcf_sum = pcf_table['PCF'].apply(lambda x: sum(x[:max_radius]))
        )
    return pcf_table

def error_median(x):
    """
    Calculate shaded plot regions by bootstrapping for error bars. 
    
    Example of use for plotting:
    pcf_median = data.groupby(variable)['PCF'].apply(
               lambda x: np.median(np.vstack(x), axis=0))
    pcf_error = data.groupby(variable)['PCF'].apply(error_median)

    """
    bootstrap = np.vstack([np.median(
        np.vstack(x)[np.random.choice(len(x), len(x), replace=True),:]
        , axis=0) for q in range(100)])
    
    return [np.percentile(bootstrap, 2.5, axis = 0), 
            np.percentile(bootstrap, 97.5, axis = 0)]
    
def nearest_neighbor(output, phenotype = 'Phenotype', k = 1):
    """
    Args:
        output: list of dicts
            from extract_data
        phenotype: string, 'Phenotype' or 'Classification'
            whether to use VECTRA Phenotype or self-generated classification
        k: int
            number of neighbors to average
    Returns:
        pandas df, 
        Mean distance from each cell type (first index level) to the 
        closest k neighbors of a given cell type (second index level), for each 
        sample (column)
    """
    
    nn_output = []
      
    knn_mean_function = lambda x: (x.apply(pd.Series.nsmallest, axis=1, n=k)
                                    .mean(axis = 1))
    knn_group_function = lambda group: (group.groupby(group.columns, axis=1)
                                             .apply(knn_mean_function)
                                             .mean())
        
    for sample in output:
        selection = sample['Data'] 
        sample_name = sample['Sample Name']
        
        labels = selection[phenotype]
        nn = spatial.distance.pdist(selection[['Cell X Position', 
                                               'Cell Y Position']])
        nn = spatial.distance.squareform(nn)
        nn = pd.DataFrame(nn, index=labels, columns=labels)
        
        missing_counts = nn.index.value_counts() < k
        if missing_counts.sum() > 0:
            print("Warning: some cell types in ", sample_name,
                  " have fewer than ", str(k), " cells.\n", 
                  missing_counts[missing_counts].index.values)

        np.fill_diagonal(nn.values, np.nan) #remove self-distances

        nn_output.append(nn.groupby(nn.index)
                         .apply(knn_group_function)
                         .stack().rename(sample_name))

    nn_output = pd.concat(nn_output, axis=1)
    nn_output.index.set_names(['From', 'To'], inplace=True)
    
    return nn_output

def run_all(output, folder_location):
    """
    Args:
        output: list of dicts
            from extract_data
        folder_location: string
            where to save files
    Returns:
        None
    """
    xy = pd.concat([x['Data'].assign(Sample = x['Sample Name']) 
                    for x in output])
    xy.rename(columns={'Sample':'Sample Name'}, inplace=True)
    xy.to_csv(folder_location+"/xy.csv")
    
    pcf_df = pcf(output, phenotype = 'Classification', count_threshold=10)
    pcf_df.to_csv(folder_location+"/pcf.csv")
     
    nn_output = []
    for k in [1,5]:
        nn_output.append(nearest_neighbor(output, phenotype = 'Classification',
                                          k = k).stack().rename(k))
    nn = pd.concat(nn_output, axis=1)
    nn.index.set_names('Sample Name', level=2, inplace=True)
    nn.to_csv(folder_location+"/nearest_neighbor.csv")
    
    return None







def interaction_subset(df, cell1, cell2, min_count=20):
    pcf_table = df[(((df['Cell_one']==cell1) & (df['Cell_two']==cell2)) | ((df['Cell_one']==cell2) & (df['Cell_two']==cell1)) ) & (df['min_count'] > min_count)]
    #pcf_table = pcf_table[pcf_table['PCF'] != 'NA']
    return pcf_table



def plot_difference(data, ax = False):
    pcf_dm = np.median(np.vstack(data['normPCF']),axis=0)
    pcf_dsem = error_median(data['normPCF'])
    
    if ax == False:
        fig, ax = plt.subplots()
        fig.set_size_inches(8,5)

    ax.plot(x_data,  pcf_dm)
    ax.fill_between(x_data, pcf_dsem[0], pcf_dsem[1],  alpha = .4, lw=0)

    ax.set_title(data['Cell_one'].iloc[0] + ' vs. ' + data['Cell_two'].iloc[0], fontsize = 18);
    ax.plot(x_data, np.ones(len(x_data)), 'k--');
    ax.set_xlabel('Radius ($\mu$m)', fontsize = 16);
    ax.set_ylabel('PCF', fontsize = 16);
    ax.tick_params(axis='both', which='major', labelsize=14);





x_data = np.arange(0, 125, .25) 






def error_median(x):
    bootstrap = np.vstack([np.median(np.vstack(x)[np.random.choice(len(x), len(x), replace=True),:],axis=0) for q in range(100)])
    return [np.percentile(bootstrap, 2.5, axis = 0), np.percentile(bootstrap, 97.5, axis = 0)]

def plot_difference(data, ax = False):
    pcf_dm = np.median(np.vstack(data['normPCF']),axis=0)
    pcf_dsem = error_median(data['normPCF'])
    
    if ax == False:
        fig, ax = plt.subplots()
        fig.set_size_inches(8,5)

    ax.plot(x_data,  pcf_dm)
    ax.fill_between(x_data, pcf_dsem[0], pcf_dsem[1],  alpha = .4, lw=0)

    ax.set_title(data['Cell_one'].iloc[0] + ' vs. ' + data['Cell_two'].iloc[0], fontsize = 18);
    ax.plot(x_data, np.ones(len(x_data)), 'k--');
    ax.set_xlabel('Radius ($\mu$m)', fontsize = 16);
    ax.set_ylabel('PCF', fontsize = 16);
    ax.tick_params(axis='both', which='major', labelsize=14);

def interaction_subset(df, cell1, cell2, min_count=0):
    pcf_table = df[(((df['Cell_one']==cell1) & (df['Cell_two']==cell2)) | ((df['Cell_one']==cell2) & (df['Cell_two']==cell1)) ) & (df['min_count'] > min_count)]
    #pcf_table = pcf_table[pcf_table['PCF'] != 'NA']
    return pcf_table

#pcf_df_tmp=pcf_df[SAMPLE_NUMBER]
#CELL_TYPES_ALIAS = ['All', 'pDCs', 'T cells','Medullary macrophages', 'Uninfected CD169+ SCS macrophages']
def plot_pcf_curves(pcf_df,cell_types_,resolution,label,out_path,radius):
  CELL_TYPES_ALIAS =['All']
  CELL_TYPES_ALIAS.extend(cell_types_)
  
  CELL_TYPES = CELL_TYPES_ALIAS
  COLOR_LIST = sns.color_palette(None, len(CELL_TYPES))
  print(COLOR_LIST)
  
  #x_data = np.arange(0, 125, .25) 
  rmax=radius/resolution
  stepsize=rmax/500
  
  x_data = np.arange(0,rmax,stepsize)
  
  # fig = plt.figure(figsize=(10, 10))
  # #outer = mpl.gridspec.GridSpec(3, 3, wspace=0.2, hspace=0.2)
  # outer = mpl.gridspec.GridSpec(int(np.ceil(np.sqrt(len(CELL_TYPES))))+1, int(np.ceil(np.sqrt(len(CELL_TYPES)))), wspace=0.4, hspace=0.4)
  # print(str(outer.ncols))
  # print(str(outer.nrows))
  # print(len(CELL_TYPES))
  # print(CELL_TYPES)
  
  pcfs_global=[]
  for i in range(len(CELL_TYPES)):
    print('i: '+str(i))
    print(CELL_TYPES[i])
    # inner = mpl.gridspec.GridSpecFromSubplotSpec(2, 1,
    #                 subplot_spec=outer[i], wspace=0.2, hspace=0.2, height_ratios=[1,2])
    # 
    # 
    # ax = plt.Subplot(fig, inner[0])
    # ax2 = plt.Subplot(fig, inner[1])
    # ax.set_title(CELL_TYPES_ALIAS[i], fontsize = 16);
    # ax2.plot(x_data, np.ones(len(x_data)), 'k--');
    # 
    pcfs_i=[]
    for j in range(len(CELL_TYPES)):
        print('j: '+str(j))
        print(CELL_TYPES[j])
        
        data = interaction_subset(pcf_df, CELL_TYPES[i], CELL_TYPES[j])
        pcfs_i.append(data)
        
        print('interaction subset computed')
        #pcf_dm = np.median(np.vstack(data['normPCF']),axis=0)
        #pcf_dsem = error_median(data['normPCF'])
        
        #ax.plot(x_data,  pcf_dm, color = COLOR_LIST[j])
        #ax.fill_between(x_data, pcf_dsem[0], pcf_dsem[1],  alpha = .4, lw=0, color = COLOR_LIST[j])    
        #ax2.plot(x_data,  pcf_dm, color = COLOR_LIST[j])
        #ax2.fill_between(x_data, pcf_dsem[0], pcf_dsem[1],  alpha = .4, lw=0, color = COLOR_LIST[j])    
    
    print('compiling DF: ')
    pcfs_i_df=pd.DataFrame(np.vstack(pcfs_i),columns = [
        'Patient', 'Sample Name', 'Cell_one', 'Cell_two', 'PCF', 'normPCF',
        'PCFsum', 'normalization', 'min_count', 'count_one', 'count_two'
        ])
    print(pcfs_i_df.head())
    
    
    pcfs_global.append(pcfs_i_df)
    
    # ax.set_ylim(3.5, 40)  # outliers only
    # ax2.set_ylim(0, 3)  # most of the data
    # ax.set_xlim([0, 50]);
    # ax2.set_xlim([0, 50]);
    # ax.spines['bottom'].set_visible(False)
    # ax2.spines['top'].set_visible(False)
    # ax.xaxis.tick_top()
    # ax.tick_params(labeltop=False)  # don't put tick labels at the top
    # ax2.xaxis.tick_bottom()
    # d = .015  # how big to make the diagonal lines in axes coordinates
    # kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
    # ax.plot((-d, +d), (-2*d, +2*d), **kwargs)        # top-left diagonal
    # ax.plot((1 - d, 1 + d), (-2*d, +2*d), **kwargs)  # top-right diagonal
    # kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
    # ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
    # ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal
    
    # if i%3 > 0:
    #     ax.set_yticks([])
    #     ax2.set_yticks([])
    # if i < 3:
    #     ax.set_xticks([])
    #     ax2.set_xticks([])        
    # 
    # fig.add_subplot(ax)
    # fig.add_subplot(ax2)

  #fig.add_subplot(111, frameon=False)
  #plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
  #plt.grid(False)
  #plt.xlabel("Radius ($\mu$m)", fontsize=16)
  #plt.ylabel("Pair Correlation Function", fontsize=16)


  # supported values are 'best', 'upper right', 'upper left', 'lower left', 'lower right', 'right', 'center left', 'center right', 'lower center', 'upper center', 'center'
  #lgd = ax2.legend(['Poisson'] + CELL_TYPES_ALIAS, bbox_to_anchor=(1.05, .8), fontsize=14);
  #lgd = ax2.legend(['Poisson'] + CELL_TYPES_ALIAS, loc='lower right', bbox_to_anchor=(0.5, 0., 0.5, 0.5), fontsize=14); 
  #lgd = ax2.legend( ['Poisson'] + CELL_TYPES_ALIAS,loc='center', 
  #           bbox_to_anchor=(0.5, 0., 0.5, 0.5),fancybox=False, shadow=False, ncol=outer.ncols)
             
  #leg = ax2.get_legend()
  #print(dir(leg))
  #colors_=['k']
  #colors_hex=[COLOR_LIST.as_hex()[x] for x in range(len(COLOR_LIST))]
  #colors_.extend(colors_hex)
  #print(colors_)
  #for l in colors_:
  #  leg.legendHandles[colors_.index(l)].set_color(colors_[colors_.index(l)])
    
  #fig.tight_layout();
  #fig.show(warn=False);
  #fig.savefig(out_path+'/'+str(label)+'_PCF-plot.pdf')
  return pcfs_global
  #return None

def pcf_AUC(pcf_df,cell_types_,cell1,label,out_path,resolution):
  CELL_TYPES = ['All']
  CELL_TYPES.extend(cell_types_)
  print(CELL_TYPES)
  CELL_TYPES_ALIAS = CELL_TYPES
 
  CELL2_LIST = CELL_TYPES
  print(CELL_TYPES)
  CELL2_LIST.pop(CELL_TYPES.index(cell1))
  print(CELL_TYPES)
  LABELS = CELL_TYPES
  print(cell1)
  print(LABELS)
  
  print('Cell 2 list:')
  print(CELL2_LIST)
  COLORS = sns.color_palette(None, len(CELL_TYPES))
  STEP_TO_UM = resolution

  def interaction_subset(df, cell1, cell2, min_count=0):
    pcf_table = df[(((df['Cell_one']==cell1) & (df['Cell_two']==cell2)) | ((df['Cell_one']==cell2) & (df['Cell_two']==cell1)) ) & (df['min_count'] > min_count)]
    #pcf_table = pcf_table[pcf_table['PCF'] != 'NA']
    return pcf_table


  def pvalue_text(data1, data2, verbose = False):
      pvalue = stats.mannwhitneyu(data1, data2, alternative='two-sided')[1]
      if verbose:
          print(pvalue)
      if pvalue < 0.001:
          return '***'
      elif pvalue < 0.01:
          return '**'
      elif pvalue < 0.05:
          return '*'
      else:
          return '$^{n.s.}$'
  print('computing results')
  results = []
  for cell2 in CELL2_LIST:
      print(cell1)
      print(cell2)
      results.append((interaction_subset(pcf_df, cell1, cell2 )['normPCF']*STEP_TO_UM).values)
      #results.append(np.array(np.array(interaction_subset(pcf_df, cell1, cell2 )['PCF_AUC'])[0])*STEP_TO_UM)

  print("results completed")
  fig, ax = plt.subplots()
  fig.set_size_inches(10,10)


  results_= [x[0] for x in results]
  print('len results_: '+str(len(results)))
  #results_ = results
  medianprops = dict(linestyle='-', linewidth=1, color='black')
  vplot = ax.violinplot(results_, showextrema=False);
  bplot = ax.boxplot(results_, widths = 0.25, whis=[5, 95], medianprops=medianprops, showfliers=False, patch_artist=True);

  for patch, line, violin, color in zip(bplot['boxes'], bplot['medians'], vplot['bodies'], COLORS):
      patch.set_facecolor(color)
      violin.set_facecolor(color)
      line.set_color('lightgray') #
  
  ax.set_xticklabels(LABELS, fontsize=13 );
  ax.get_xaxis().set_tick_params(direction='out')

  
  #ax.set_title(cell1 + ' vs. ' + cell2, fontsize=20);
  ax.set_ylabel(str(cell1)+' AUC', fontsize=20);
  ax.set_title(str(cell1)+' Interaction', fontsize=20);
  #ax.add_patch(mpl.patches.Rectangle((-1, -50), 8.5, 6000, alpha = 0.2, fc = plt.cm.tab10(7) ))
  innit_val=0
  
  
  for l in LABELS:
    
    ax.plot([0.1, 0.38*(1+innit_val)], [0.8, 0.8], 'k', transform=ax.transAxes)
    #print(str(bplot['medians'][0].get_data()))
    #ax.plot([0.1, bplot['medians'][0].get_data()[1][0]], [0.8, 0.8], 'k', transform=ax.transAxes)
    #ax.hlines(y=bplot['medians'][0].get_data()[1][0],xmin=LABELS[0],xmax=LABELS[len(LABELS)])
    ax.text(0.1+0.25*(innit_val)/(len(LABELS)/4), 0.805+0.05*innit_val/(len(LABELS)/4), pvalue_text(results_[0],results_[innit_val], verbose=True), horizontalalignment='center', transform=ax.transAxes)
    innit_val += 1
    
  #ax.plot([0.1, 0.62], [0.9, 0.9], 'k', transform=ax.transAxes)
  #ax.text(0.5, 0.905, pvalue_text(results_[0],results_[2], verbose=True), horizontalalignment='center', transform=ax.transAxes)

  #ax.plot([0.62, 0.72], [0.8, 0.8], 'k', transform=ax.transAxes)
  #ax.text(0.67, 0.805, pvalue_text(results_[2],results_[1], verbose=True), horizontalalignment='center', transform=ax.transAxes)

  #ax.plot([0.1, 0.92], [0.7, 0.7], 'k', transform=ax.transAxes)
  #ax.text(0.4, 0.7, pvalue_text(results_[0],results_[3], verbose=True), horizontalalignment='center', transform=ax.transAxes)

  fig.tight_layout();
  #fig.savefig(out_path+'/'+str(label)+'_AUC_PCF-Violins.pdf')



  results_df = pd.DataFrame(np.vstack(results_).T)
  results_df.columns=CELL2_LIST
  #results_df.to_csv(out_path+'PCFs.csv')
  
  return results_df


  
  



