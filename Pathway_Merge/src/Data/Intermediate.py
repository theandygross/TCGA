'''
Created on Jun 23, 2013

Module for loading data from processed files.
This is separate from Data.Firehose due to its dependence on upstream
processing pipelines or other functions in the package.

@author: agross
'''
import os as os
import numpy as np
import pandas as pd

import Data.Firehose as FH

from Processing.Tests import kruskal_pandas
from Processing.Helpers import true_index, bhCorrection


def get_beta_values(data_path, cancer, patients=None, tissue_code='01'):
    '''
    Retrieve methylation beta-values from my pre-processed file.  
    TCGA has a lot more columns that eat up memory, so I parse out the 
    beta-values in preprocess_methylation.py.  
    This file still has all of the probes by pateints so it still eats my
    kill a computer without a lot of memory (takes ~2GB for HNSC). 
    '''
    path = '{}/ucsd_processing/{}/methylation450/'.format(data_path, cancer)
    t = pd.read_table(path + 'beta_values.txt', skiprows=[1], index_col=[0])
    t = t.rename(columns=lambda s: s if s!=t.columns[0] else 'symbol')
    t = t.set_index('symbol', append=True)
    t = t.swaplevel(0,1)
    t = FH.fix_barcode_columns(t, patients, tissue_code)
    return t

def read_methylation(data_path, cancer, patients=None, tissue_code='01'):
    '''
    Reads in gene by patient methylation matrix.  
    Here genes levels are pre-computed by taking the principal component of 
    the probes in the gene (see preprocess_methylation.py for details).  

    tissue_code: ['01','11','All']  #if all returns MultiIndex
    '''
    path = '{}/ucsd_processing/{}/'.format(data_path, cancer)
    data_types = filter(lambda f: f[:11] == 'methylation', os.listdir(path))
    data_type = sorted([d for d in data_types 
                        if os.path.isfile(path + d + '/meta_probes.csv')])[-1]
    meth = pd.read_csv(path +  data_type + '/meta_probes.csv', index_col=0)
    ft = pd.MultiIndex.from_tuples
    meth.columns = ft(map(lambda s: eval(s,{},{}), meth.columns))
    
    if patients is not None:
        meth  = meth.ix[:, patients]
    meth = meth.astype(np.float)
    meth = meth.dropna(thresh=100)
    if tissue_code != 'All':
        meth = meth.T.xs(tissue_code, level=1).T #pandas bug
    return meth

def rna_filter(cn, val, rna):
    '''
    Filter copy number events with rna expression data.
    Here we test whether the event is associated with a subsequent
    change in expression in those patients. 
    
    cn: copy number matrix, should have a MultiIndex, with the gene name
        in the last level
    val: value of the copy number to test in [-2, -1, 1, 2] 
    '''
    assert val in [-2, -1, 1, 2]
    change = pd.DataFrame({g: kruskal_pandas(vec == val, rna.ix[g[-1]])
                           for g,vec in cn.iterrows() 
                           if g[-1] in rna.index}).T
    q_vals = bhCorrection(change.p)
    filtered = cn.ix[true_index(q_vals < .1)]
    return filtered

def get_gistic_genes(data_path, cancer, filter_with_rna=True, 
                     collapse_on_bands=True, min_patients=5):
    '''
    Gets a matrix of events for high grade amplifications and homozygous 
    deletions. 
    We filter down this list by asserting that a copy number event corresponds
    with a resultant expression change. 
    The final matrix merges gene-level events on the same band to combine
    redundant events and reduce the test space.     
    '''
    gistic = FH.get_gistic_gene_matrix(data_path, cancer, '01')
    deletion = gistic[(gistic == -2).sum(1) > min_patients]
    amp = gistic[(gistic == 2).sum(1) > min_patients]
    ft = pd.MultiIndex.from_tuples #rediculously long pandas names
    deletion.index = ft([('Deletion', s[0], s[2]) for s in deletion.index])
    amp.index = ft([('Amplification', s[0], s[2]) for s in amp.index])
    
    if filter_with_rna:
        rna = FH.read_rnaSeq(data_path, cancer)
        deletion = rna_filter(deletion, -2, rna)
        amp = rna_filter(amp, 2, rna)
   
    cna_genes = amp.append(deletion)
    if collapse_on_bands == False:
        return cna_genes
    
    cna_genes = pd.DataFrame({(a[0], a[1], tuple(b.index.get_level_values(2))): 
                           b.mean().round() for a,b in 
                           cna_genes.groupby(level=[0,1])}).T
    cna_genes.index = pd.MultiIndex.from_tuples(cna_genes.index)
    return cna_genes

def get_gistic(data_path, cancer, filter_with_rna=True, 
               collapse_on_bands=True, min_patients=5):
    '''
    Get the combined GISTIC feature matrix for testing. 
    '''
    lesions = FH.get_gistic_lesions(cancer, data_path)
    cna_genes = get_gistic_genes(data_path, cancer, filter_with_rna, 
                                 collapse_on_bands, min_patients)
    cna = cna_genes.append(lesions)
    return cna