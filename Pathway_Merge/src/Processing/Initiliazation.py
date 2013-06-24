'''
Created on May 8, 2013

@author: agross
'''

from Data.Containers import Dataset
from Processing.Helpers import true_index, frame_svd, extract_pc

import pandas as pd
import numpy as np

def initialize_real(cancer, run, data_type, patient_filter=None, drop_pc1=False):
    data = Dataset(cancer, run, data_type)
    if data_type == 'mRNASeq':
        data.df = data.df.groupby(by=lambda n: n.split('|')[0]).mean()
    
    '''Normalize data and calculate principal components'''
    if type(patient_filter) == pd.Series:
        data.df  = data.df.ix[:,true_index(patient_filter)].dropna(axis=1, thresh=100)
        
    df = data.df
    norm = ((df.T - df.mean(1)) / df.std(1)).T
    U,S,vH = frame_svd(norm)
    data.pc1, data.pc2 = vH[0], vH[1]
    data.loading1, data.loading2 = U[0], U[1]
    
    '''Reconstruct data without first pc'''
    S_n = S.copy()
    if drop_pc1:
        S_n[0] = 0
    rest = U.dot(pd.DataFrame(np.diag(S_n)).dot(vH.T))
    
    '''Extract PCs for all gene sets'''
    pc = {p: extract_pc(rest.ix[g]) for p,g, in run.gene_sets.iteritems()}
    pc = pd.DataFrame({p: s for p,s in pc.iteritems() if s}).T
    
    data.pc_genes = pc.gene_vec
    data.pct_var = pc.pct_var
    data.features = pd.DataFrame(pc.pat_vec.to_dict()).T
    return data