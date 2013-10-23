'''
Created on Aug 27, 2013

@author: agross
'''
import pandas as pd
import numpy as np

import Processing.Helpers as H
import Stats.Survival as Surv
from Stats.Scipy import rev_kruskal, fisher_exact_test
from Stats.Scipy import spearman_pandas
from Initialization.InitializeReal import exp_change

def fc(hit_vec, response_vec):
    f = response_vec.groupby(hit_vec).median()
    return f[0] > f[1]

def mut_filter(df, rate, binary_cutoff=12):
    '''
    Filter out mutation features, ensuring that a feature
    is not entirely an artifact of mutation rate.
    '''
    df = df[df.sum(1) >= binary_cutoff]
    cc = H.screen_feature(rate, rev_kruskal, df)
    
    fc_apply = lambda s: fc(s, rate)
    direction = df.apply(fc_apply, axis=1)
    direction.name = 'direction'
    
    cc = cc.join(direction)
    cc = cc[cc.direction==False]
    
    df = df.ix[H.true_index(cc.p > .01)]
    df = df.dropna(axis=1)
    return df

def cn_filter(df, binary_cutoff=12):
    '''
    Extract copy number features.
    '''
    del_df = (df.ix['Deletion'].dropna(1) < 0).astype(int)
    del_df = del_df[del_df.sum(1) >= binary_cutoff]
    del_df.index = del_df.index.droplevel(1)
    del_df = del_df.T
    amp_df = (df.ix['Amplification'].dropna(1) > 0).astype(int)
    amp_df = amp_df[amp_df.sum(1) >= binary_cutoff]
    amp_df.index = amp_df.index.droplevel(1)
    amp_df = amp_df.T
    return amp_df, del_df

def process_real(df):
    '''
    Process real valued feature into binary feature.
    '''
    df_c = df.copy()
    df_c = df_c.apply(lambda s: H.to_quants(s, std=1), axis=1)
    df_c = df_c > 0
    if type(df.index) == pd.MultiIndex:
        df_c.index = map(lambda s: '_'.join(s), df_c.index)
    return df_c.T


def binary_filter_fx(s):
    vc = s.value_counts()
    if len(vc) != 2:
        return -1
    else:
        return vc.min()

def filter_binary(df, cutoff):
    df = df.dropna(how='all')
    vc = df.apply(binary_filter_fx, axis=1)
    binary = df[vc > cutoff]
    return binary

def binarize_feature(f):
    '''
    Binarize a feature to minimize the difference in sum of squares between 
    the two resulting groups.  
    '''
    f = f - f.mean()
    f2 = (f.order()**2)
    split = f.ix[(f2.cumsum() - (f2.sum() / 2.)).abs().idxmin()]
    return f > split

def remove_redundant_pathways(pathways, background, cutoff=.7,
                              binarize=False):
    '''
    Screens out redundant pathways with high correlation above _cutoff_.
    Pathways are ranked based on lack of correlation to the background signal.
    Then if two pathways have high correlation the lower ranked pathway is 
    removed.  
    '''
    bg = H.screen_feature(background, spearman_pandas, pathways)
    dd = pathways.ix[bg.index[::-1]].T.corr()
    dd = pd.DataFrame(np.triu(dd, 1), dd.index, dd.index)
    dd = dd.replace(0, np.nan).stack()
    drop = dd[dd.abs() > cutoff].index.get_level_values(1)
    pathways_to_keep = pathways.index.diff(drop.unique())
    pathways =  pathways.ix[pathways_to_keep]
    if binarize is False:
        return pathways
    else:
        binary_pathways = pathways.apply(binarize_feature, 1)
        return binary_pathways
    

def extract_diff_exp_rna(rna, n=300, binarize=False):
    '''
    Pull the most differentially expressed genes from the rna expression 
    object. 
    '''
    genes = rna.features.ix[['real','binary']].index.get_level_values(1)
    dd = rna.df.ix[genes].dropna()
    rr = dd.apply(exp_change, 1)
    d2 = dd.ix[rr.sort('F').index[-n:]].xs('01',1,1)
    if binarize is False:
        return d2
    else:
        real_genes = rna.features.ix['real'].index
        tf = lambda s: binarize_feature(s) if s.name in real_genes else s < -1
        d3 = d2.apply(tf, 1)
        return d3
    
def corrections(vec):
    '''
    Correct p-values multiple ways along multi-index.
    '''
    bonf_all = vec * len(vec)
    bonf_within = vec.groupby(level=0).apply(lambda s: s*len(s))
    
    bh_all = H.bhCorrection(vec)
    bh_within = vec.groupby(level=0).apply(H.bhCorrection).order()
    
    two_step = bh_within * len(vec.groupby(level=0).size())
    q = pd.concat([vec, bh_within, bh_all, bonf_all, bonf_within, two_step],
                  keys=['uncorrected', 'bh_within', 'bh_all', 'bonf_all', 'bonf_within',
                        'two_step'], axis=1)
    return q

class Screen(object):
    '''
    Object to hold data and results for survival screen.
    '''
    
    def __init__(self, mut, cn, rna, mirna, clinical_df): 
        
        '''Process Gene / Pathway Expression'''
        pathways = remove_redundant_pathways(rna.pathways, rna.global_vars.background, 
                                             binarize=True)
        rna_gene_df = extract_diff_exp_rna(rna, n=300, binarize=True)
        self.rna_df = pd.concat([rna_gene_df, pathways]).T
        
        '''Process miRNA Expression'''
        mirna_binarized = mirna.features.ix['real'].apply(binarize_feature, 1)
        self.mirna_df = pd.concat([mirna.features.ix['binary'], mirna_binarized]).T
        
        '''Process mutation data'''
        self.rate = mut.df.sum()
        self.mut_df = mut.features
        
        '''Process CNA data'''
        self.cna_df = cn.features
        
        '''Process Clinical Data'''
        self.clinical_df  = clinical_df
        
    def get_patient_set(self, filters):
        f1 = list(filters)
        filter_df = pd.concat(f1, axis=1)
        clinical_filter = filter_df.dropna().sum(1) == 0
        keepers_o = H.true_index(clinical_filter)
        keepers_o = keepers_o.intersection(self.mut_df.columns)
        keepers_o = keepers_o.intersection(self.cna_df.columns)
        return keepers_o
        
    def get_data(self, keepers_o, cutoff=12):
        mut_df = mut_filter(self.mut_df.ix[:,keepers_o], self.rate, 
                            cutoff).T
        amp_df, del_df = cn_filter(self.cna_df.ix[:, keepers_o], cutoff)
        cna_df = pd.concat([del_df, amp_df], keys=['del','amp'], axis=1)
        cna_df.columns = map(lambda i: '_'.join(i), cna_df.columns)
        
        df = pd.concat({'clinical': self.clinical_df,
                        'mutation': mut_df,
                        'cna': cna_df,
                        'rna': self.rna_df,
                        'mirna': self.mirna_df}, axis=1).T
        df = df.ix[:, keepers_o]
        df = filter_binary(df, cutoff)
        return df
    
class ScreenResult(object):
    '''
    Object to hold results for survival screen iteration.
    '''
    def __init__(self, results, pairs, full, univariate, patients, df):
        self.patients = patients
        self.features = univariate.index
        
        self.univariate = univariate
        self.full = full
        self.pairs = pairs
        self.results = results 
        self.df = df
        
    def __repr__(self):
        hits = sum(self.full.p.bh_all < .1)
        best = self.results.p.uncorrected.order().index[0]
        best_q = self.results.p.bh_all.ix[best]
        s = 'Screen Result: \n'
        s += '    {} events tested across {} patients\n'.format(len(self.features), len(self.patients))
        s += '    {} events were significant above .1 FDR.\n'.format(hits)
        s += '    {} pairs of events were significantly overlapping.\n'.format(len(self.pairs))
        s += '    {} was the top association with a q value of {}.\n'.format(best, best_q)
        return s

