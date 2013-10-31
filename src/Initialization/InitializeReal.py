#! /usr/bin/env python

import os as os
import pickle as pickle

import pandas as pd
from matplotlib.pylab import subplots
from numpy import diag
import Stats.Scipy as Tests

from Data.Containers import Dataset
from Processing.Helpers import frame_svd, extract_pc
from Processing.Helpers import true_index, screen_feature
from Stats.Scipy import pearson_pandas

import Data.Intermediate as IM

def exp_change(s):
    '''
    Calculates an anova for the change in expression across a variable
    on the second level of a MultiIndex. (eg. tumor/normal).
    '''
    return Tests.anova(pd.Series(s.index.get_level_values(1), s.index), s)


def extract_pc_filtered(df, pc_threshold=.2, filter_down=True):
    '''
    First pre-filters for patients with no tumor/normal change.
    Then normalizes by normals. 
    '''
    if ('11' in df.columns.levels[1]) and filter_down:
        tt = df.xs('11', axis=1, level=1)
        rr = df.apply(exp_change, 1).sort('p')
        m, s = tt.mean(1), tt.std(1)
        df_n = df.xs('01', axis=1, level=1)
        df_n = ((df_n.T - m) / s).T
        df_n = df_n.ix[true_index(rr.p < .05)]
    else:  # No matched normals
        df_n = df.xs('01', axis=1, level=1)
        df_n = ((df_n.T - df_n.mean(1)) / df_n.std(1)).T
    pc = extract_pc(df_n, pc_threshold, standardize=False)
    return pc

def extract_geneset_pcs(df, gene_sets, filter_down=True):
    '''Extract PCs for all gene sets.'''
    pc = {p: extract_pc_filtered(df.ix[g].dropna(), .2, filter_down) for p, g,
          in gene_sets.iteritems()}
    pc = pd.DataFrame({p: s for p, s in pc.iteritems() if s}).T
    pat_vec = pd.DataFrame(pc.pat_vec.to_dict()).T
    return pc.gene_vec, pc.pct_var, pat_vec

def get_mirna_features(df):
    binary = (df[(df < -1).sum(1) > (df.shape[1] / 2)] >= -1) * 1.
    binary = binary[binary.sum(1).isin(range(20, df.shape[1] / 2))]
    
    real = df[((df.max(1) - df.min(1)) > 2)]
    real = real.ix[(real == -3).sum(1) < real.shape[1] / 2.]
    features = pd.concat([real, binary], keys=['real', 'binary'])
    return features

def extract_features(df):
    df_n = df.xs('01', level=1, axis=1)
    binary = df_n > -1
    binary = binary[binary.sum(1).isin(range(20, df.shape[1] / 2))]
    rr = df.ix[binary.index].apply(exp_change, 1)
    binary = binary.ix[true_index(rr.p < .05)]
    
    real = df_n.ix[df_n.index.diff(binary.index)]
    singles = real[((real.max(1) - real.min(1)) > 1)]
    singles = singles[(singles.std(1) > .25)]
    ch = df.ix[singles.index].apply(exp_change, 1)
    singles = df_n.ix[true_index(ch.p < .01)]
    return binary, singles, real
        
class RealDataset(Dataset):
    '''
    Inherits from Dataset class.  Adds some added processing for real valued
    data.
    '''
    def __init__(self, run, cancer, data_type, patients=None, drop_pc1=False,
                 create_real_features=True, create_meta_features=True,
                 filter_down=True, draw_figures=False):
        '''
        '''
        Dataset.__init__(self, cancer.path, data_type, compressed=True)
        self.df = IM.read_data(run.data_path, cancer.name, data_type,
                               tissue_code='All')
        if patients is not None:
            self.df = self.df.ix[:, patients].dropna(axis=1, how='all')
            self.patients = patients
        else:
            self.patients = self.df.xs('01', 1, 1).columns
        
        self.global_vars = pd.DataFrame(index=self.patients)
        self.features = {}
        self.global_loadings = pd.DataFrame(index=self.df.index)
        self._calc_global_pcs(drop_pc1)
        
        if create_real_features is True:
            self._get_real_features()
        
        if create_meta_features is True:
            self._get_meta_features(run.gene_sets, filter_down)
            
        self.features = pd.concat(self.features)
            
        if draw_figures is True:
            self._creat_pathway_figures()
            
    def _get_real_features(self):
        binary, singles, real = extract_features(self.df)
        background_df = real.ix[real.index.diff(singles.index)].dropna()
        background = extract_pc(background_df, 0)
        ss = screen_feature(background['pat_vec'], pearson_pandas, singles)
        singles = singles.ix[ss.p > 10e-5]
        
        singles = ((singles.T - singles.mean(1)) / singles.std(1)).T
        U, S, pc = frame_svd(singles)
        
        self.features['binary'] = binary
        self.features['real'] = singles
        self.global_vars['background'] = background['pat_vec']
        self.global_vars['filtered_pc1'] = pc[0]
        self.global_vars['filtered_pc2'] = pc[1]
        self.global_loadings['background'] = background['gene_vec']
        self.global_loadings['filtered_pc1'] = U[0]
        self.global_loadings['filtered_pc2'] = U[1]
        
    def _get_meta_features(self, gene_sets, filter_down):
        gs = extract_geneset_pcs(self.df, gene_sets, filter_down)
        self.loadings, self.pct_var, pathways = gs
        if hasattr(self.global_vars, 'background'):
            r = screen_feature(self.global_vars.background, pearson_pandas,
                               pathways)
            pathways = pathways.ix[r.p > 10e-5]
        pathways = ((pathways.T - pathways.mean(1)) / pathways.std(1)).T
        U, S, pc = frame_svd(pathways)
        
        self.pathways = pathways
        self.features['pathways'] = pathways
        self.global_vars['pathway_pc1'] = pc[0]
        self.global_vars['pathway_pc2'] = pc[1]
        self.global_loadings['pathway_pc1'] = U[0]
        self.global_loadings['pathway_pc2'] = U[1]
        
    
    def _calc_global_pcs(self, drop_pc1=False):
        '''
        Normalize data and calculate principal components. If drop_pc1 is
        set to True, also reconstructs the normalized data without the
        first PC. 
        '''
        df = self.df.xs('01', axis=1, level=1)
        norm = ((df.T - df.mean(1)) / df.std(1)).T
        U, S, vH = frame_svd(norm)
        self.global_vars['pc1'] = vH[0]
        self.global_vars['pc2'] = vH[1]
        self.global_loadings['pc1'] = U[0]
        self.global_loadings['pc2'] = U[1]        
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
            fig, axs = subplots(1, 2, figsize=(10, 3))
            self.loadings.ix[feature].order().plot(kind='bar', ax=axs[0])
            axs[0].annotate('Explained Variation: %.2f' % self.pct_var[feature],
                            (.03, .97), xycoords='axes fraction', ha='left', va='top')
            axs[0].set_ylabel('Eigen-Patient Loading')
            self.features.ix[feature].hist(ax=axs[1])
            axs[1].set_xlabel('Eigen-Gene Loading')
            fig.tight_layout()
            fig.savefig(pathway_plot_folder + feature)
            
        
def initialize_real(cancer_type, report_path, data_type, patients=None,
                    drop_pc1=False, create_real_features=True,
                    create_meta_features=True, filter_down=False,
                    draw_figures=False, save=True):
    '''
    Initialize real-valued data for down-stream analysis.
    '''
    run = pickle.load(open(report_path + '/RunObject.p', 'rb'))
    cancer = run.load_cancer(cancer_type)
    
    if data_type is 'miRNASeq':
        create_meta_features = False
        draw_figures = False
    
    data = RealDataset(run, cancer, data_type, patients, drop_pc1,
                       create_real_features, create_meta_features, filter_down,
                       draw_figures)

    if save is True:
        data.save()
    return data
