#! /usr/bin/env python

import os as os
import pickle as pickle

import pandas as pd
from matplotlib.pylab import subplots
from numpy import diag
import Stats.Scipy as Tests

from Data.Containers import Dataset
from Processing.Helpers import frame_svd, extract_pc
from Processing.Helpers import true_index

import Data.Intermediate as IM

def exp_change(s):
    '''
    Calculates an anova for the change in expression across a variable
    on the second level of a MultiIndex. (eg. tumor/normal).
    '''
    return Tests.anova(pd.Series(s.index.get_level_values(1), s.index), s)


def extract_pc_filtered(df, pc_threshold=.2):
    '''
    First pre-filters for patients with no tumor/normal change.
    Then normalizes by normals. 
    '''
    if '11' in df.columns.levels[1]:
        tt = df.xs('11', axis=1, level=1)
        rr = df.apply(exp_change, 1).sort('p')
        m, s = tt.mean(1), tt.std(1)
        df_n = df.xs('01', axis=1, level=1)
        df_n = ((df_n.T - m) / s).T
        df_n = df_n.ix[true_index(rr.p < .05)]
    else: #No matched normals
        df_n = df.xs('01', axis=1, level=1)
        df_n = ((df_n.T - df_n.mean(1)) / df_n.std(1)).T
    pc = extract_pc(df_n, pc_threshold, standardize=False)
    return pc

def extract_geneset_pcs(df, gene_sets, filter_down=True):
    '''Extract PCs for all gene sets.'''
    pc = {p: extract_pc_filtered(df.ix[g].dropna()) for p,g, 
          in gene_sets.iteritems()}
    pc = pd.DataFrame({p: s for p,s in pc.iteritems() if s}).T
    pat_vec = pd.DataFrame(pc.pat_vec.to_dict()).T
    return pc.gene_vec, pc.pct_var, pat_vec

        
class RealDataset(Dataset):
    '''
    Inherits from Dataset class.  Adds some added processing for real valued
    data.
    '''
    def __init__(self, run, cancer, data_type, patients=None, drop_pc1=False,
                 create_meta_features=True, draw_figures=False):
        '''
        '''
        Dataset.__init__(self, cancer.path, data_type, compressed=True)
        self.df = IM.read_data(run.data_path, cancer.name, data_type, 
                               tissue_code='All')
        
        self._calc_global_pcs(patients, drop_pc1)
        
        if create_meta_features is True:
            gs = extract_geneset_pcs(self.df, run.gene_sets)
            self.loadings, self.pct_var, self.features = gs
            
        if draw_figures is True:
            self._creat_pathway_figures()
    
    def _calc_global_pcs(self, patients=None, drop_pc1=False):
        '''
        Normalize data and calculate principal components. If drop_pc1 is
        set to True, also reconstructs the normalized data without the
        first PC. 
        '''
        df = self.df.xs('01', axis=1, level=1)
        if patients is not None:
            df = df.ix[patients]
        norm = ((df.T - df.mean(1)) / df.std(1)).T
        U,S,vH = frame_svd(norm)
        self.pc1, self.pc2 = vH[0], vH[1]
        self.loading1, self.loading2 = U[0], U[1]
        
        if drop_pc1 is True:
            S_n = S.copy()
            S_n[0] = 0
            norm = U.dot(pd.DataFrame(diag(S_n)).dot(vH.T))
            
        return norm
    
    def _create_pathway_figures(self):
        '''
        Create figures for all gene sets features in a Dataset object.
        '''
        pathway_plot_folder = self.path + '/Figures/PathwayPlots/'
        if not os.path.isdir(pathway_plot_folder):
            os.makedirs(pathway_plot_folder)
            
        for feature in self.features.index:
            fig, axs = subplots(1,2, figsize=(10,3))
            self.loadings.ix[feature].order().plot(kind='bar', ax=axs[0])
            axs[0].annotate('Explained Variation: %.2f' % self.pct_var[feature], 
                            (.03,.97), xycoords='axes fraction', ha='left', va='top')
            axs[0].set_ylabel('Eigen-Patient Loading')
            self.features.ix[feature].hist(ax=axs[1])
            axs[1].set_xlabel('Eigen-Gene Loading')
            fig.tight_layout()
            fig.savefig(pathway_plot_folder + feature)
            
        
def initialize_real(cancer_type, report_path, data_type, patients=None, 
                    drop_pc1=False, create_meta_features=True, 
                    draw_figures=False, save=True):
    '''
    Initialize real-valued data for down-stream analysis.
    '''
    run = pickle.load(open(report_path + '/RunObject.p', 'rb'))
    cancer = run.load_cancer(cancer_type)
    
    if data_type is 'miRNASeq':
        create_meta_features = False
        draw_figures=False
    
    data = RealDataset(run, cancer, data_type, patients, drop_pc1, 
                       create_meta_features, draw_figures)
        
    if save is True:
        data.save()
    return data