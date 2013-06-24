#! /usr/bin/env python

import os as os
import sys as sys
import pickle as pickle

from pandas import DataFrame
from matplotlib.pylab import subplots
from numpy import diag

from Data.Containers import Dataset
from Processing.Helpers import frame_svd, extract_pc


def extract_geneset_pcs(df, gene_sets):
    '''Extract PCs for all gene sets.'''
    pc = {p: extract_pc(df.ix[g]) for p,g, in gene_sets.iteritems()}
    pc = DataFrame({p: s for p,s in pc.iteritems() if s}).T
    pat_vec = DataFrame(pc.pat_vec.to_dict()).T
    return pc.gene_vec, pc.pct_var, pat_vec
    
def creat_pathway_figures(data):
    '''Create figures for all gene sets features in a Dataset object.'''
    pathway_plot_folder = data.path + '/Figures/PathwayPlots/'
    if not os.path.isdir(pathway_plot_folder):
        os.makedirs(pathway_plot_folder)
        
    for feature in data.features.index:
        fig, axs = subplots(1,2, figsize=(10,3))
        data.loadings.ix[feature].order().plot(kind='bar', ax=axs[0])
        axs[0].annotate('Explained Variation: %.2f' % data.pct_var[feature], 
                        (.03,.97), xycoords='axes fraction', ha='left', va='top')
        axs[0].set_ylabel('Eigen-Patient Loading')
        data.features.ix[feature].hist(ax=axs[1])
        axs[1].set_xlabel('Eigen-Gene Loading')
        fig.tight_layout()
        fig.savefig(pathway_plot_folder + feature)

def initialize_real(cancer_type, report_path, data_type, patient_filter=None, 
                    drop_pc1=False, create_meta_features=True, draw_figures=False):
    '''
    Initialize real-valued data for down-stream analysis.
    '''
    run = pickle.load(open(report_path + '/RunObject.p', 'rb'))
    cancer = run.load_cancer(cancer_type)
    data = Dataset(cancer, run, data_type)
    
    #Normalize data and calculate principal components
    df = data.df.xs('01', axis=1, level=1)
    if patient_filter is not None:
        df = df.ix[patient_filter]
    norm = ((df.T - df.mean(1)) / df.std(1)).T
    U,S,vH = frame_svd(norm)
    data.pc1, data.pc2 = vH[0], vH[1]
    data.loading1, data.loading2 = U[0], U[1]
    
    if drop_pc1 is True:
        S_n = S.copy()
        S_n[0] = 0
        norm = U.dot(DataFrame(diag(S_n)).dot(vH.T))
        
    if data_type is 'miRNASeq':
        create_meta_features = False
        draw_figures=False
        
    #return norm
        
    if create_meta_features is True:
        gs = extract_geneset_pcs(norm, run.gene_sets)
        data.loadings, data.pct_var, data.features = gs
    
    if draw_figures is True:
        creat_pathway_figures(data)
    
    data.save()
    return data