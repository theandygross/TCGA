#! /usr/bin/env python

import os as os
import sys as sys
import pickle as pickle

import pandas as pd
from matplotlib.pylab import subplots
from numpy import diag

from Data.Containers import Dataset
from Processing.Helpers import frame_svd, extract_pc

def initialize_rppa(run, cancer):
    '''Load in run and data'''
    data = Dataset(cancer, run, 'RPPA')
    
    pc = data.df.groupby(level=0).apply(extract_pc)
    pc = pd.DataFrame(pc.dropna().to_dict()).T
    pc = pc[data.df.groupby(level=0).size() > 1]
    data.phos_loadings = pc.gene_vec
    data.phos = pd.DataFrame(pc.pat_vec.to_dict()).T
    
    '''Normalize data and calculate principal components'''
    df = data.df
    norm = ((df.T - df.mean(1)) / df.std(1)).T
    U,S,vH = frame_svd(norm)
    data.pc1, data.pc2 = vH[0], vH[1]
    data.loading1, data.loading2 = U[0], U[1]
    
    '''Reconstruct data without first pc'''
    S_n = S.copy()
    S_n[0] = 0
    rest = U.dot(pd.DataFrame(diag(S_n)).dot(vH.T))
    
    '''Extract PCs for all gene sets'''
    pc = {p: extract_pc(rest.ix[g]) for p,g, in run.gene_sets.iteritems() 
          if len(set(rest.ix[g].index.get_level_values(0))) > 4}
    pc = pd.DataFrame({p: s for p,s in pc.iteritems() if s}).T
    
    data.pc_genes = pc.gene_vec
    data.pct_var = pc.pct_var
    data.pathways = pd.DataFrame(pc.pat_vec.to_dict()).T
    df = data.df.copy()
    df.index = pd.Index(map(tuple, df.index))

    data.features = pd.concat([df, data.phos, data.pathways], 
                            keys=['protiens','phos_pc','pathways'])
    
    '''Create figures for pathway gene sets'''
    pathway_plot_folder = data.path + '/Figures/PathwayPlots/'
    if not os.path.isdir(pathway_plot_folder):
        os.makedirs(pathway_plot_folder)
        
    for feature in data.pathways.index:
        fig, axs = subplots(1,2, figsize=(10,6))
        data.pc_genes.ix[feature].order().plot(kind='bar', ax=axs[0])
        axs[0].annotate('Explained Variation: %.2f' % data.pct_var[feature], 
                        (.03,.97), xycoords='axes fraction', ha='left', va='top')
        axs[0].set_ylabel('Eigen-Patient Loading')
        data.pathways.ix[feature].hist(ax=axs[1])
        axs[1].set_xlabel('Eigen-Gene Loading')
        try:
            fig.tight_layout()
            fig.savefig(pathway_plot_folder + feature)
        except:
            print feature
            
    for feature in data.phos.index:
        fig, axs = subplots(1,2, figsize=(10,6))
        data.phos_loadings.ix[feature].order().plot(kind='bar', ax=axs[0])
        axs[0].set_ylabel('Eigen-Patient Loading')
        data.phos.ix[feature].hist(ax=axs[1])
        axs[1].set_xlabel('Eigen-Gene Loading')
        try:
            fig.tight_layout()
            fig.savefig(pathway_plot_folder + feature)
        except:
            print feature
        
          
    data.save()
    
if __name__ == "__main__":
    report_path = sys.argv[1]
    cancer_type = sys.argv[2]
    data_type = 'RPPA'
   
    run = pickle.load(open('/'.join([report_path, 'RunObject.p']), 'rb'))
    cancer = run.load_cancer(cancer_type) 
    initialize_rppa(run, cancer)
    sys.exit(0)