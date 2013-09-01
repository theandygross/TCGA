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

from Stats.Scipy import kruskal_pandas
from Processing.Helpers import true_index, bhCorrection
from Processing.Helpers import frame_svd

def get_beta_values(data_path, cancer, patients=None, tissue_code='All'):
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
    t = t.sort_index() #think this is a bug that it needs to be sorted
    t = FH.fix_barcode_columns(t, patients, tissue_code)
    return t

def read_methylation(data_path, cancer, patients=None, tissue_code='01'):
    '''
    Reads in gene by patient methylation matrix.  
    Here genes levels are pre-computed by taking the principal component of 
    the probes in the gene (see preprocess_methylation.py for details).  

    tissue_code: ['01','11','All']  #if all returns MultiIndex
    '''
    path = '{}ucsd_processing/{}/'.format(data_path, cancer)
    data_types = filter(lambda f: f[:11] == 'methylation', os.listdir(path))
    data_type = sorted([d for d in data_types 
                        if os.path.isfile(path + d + '/meta_probes.csv')])[-1]
    meth = pd.read_csv(path +  data_type + '/meta_probes.csv', index_col=0)
    ft = pd.MultiIndex.from_tuples
    meth.columns = ft(map(lambda s: eval(s,{},{}), meth.columns))
    
    if patients is not None:
        meth  = meth.ix[:, patients]
    meth = meth.astype(np.float)
    meth = meth.dropna(how='all')
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

def build_meta_matrix(gene_sets, gene_matrix, min_size=4, set_filter=None):
    '''
    Builds meta-gene matrix from gene-set definitions and single gene matrix.
    
    Input
        gene_sets: dict {name: genes} or list of gene-set definitions
        gene_matrix: DataFrame of single gene profiles.
        min_size: minumum number of patient alterations for meta-feature 
                  (default 4)
        set_filter: filter to run on meta-features after they are constructed
        
    Returns
        meta_matrix
    '''
    if not isinstance(gene_sets, dict):
        gene_sets = dict(list(enumerate(gene_sets)))
    meta_matrix = pd.DataFrame({group: gene_matrix.ix[genes].sum(0) 
                                for group,genes in gene_sets.items()})
    if set_filter is not None: 
        meta_matrix = meta_matrix.apply(set_filter, axis=0)
    keepers = meta_matrix.sum().isin(range(min_size, 
                                           len(meta_matrix) - min_size))
    meta_matrix = meta_matrix.ix[:,keepers].T
    return meta_matrix

def get_mutation_rates(data_path, cancer, patients=None):
    '''
    Get mutation rate data from MutSig processing pipeline. This function
    depends on the current Firehose output of this program as of July 2013.
    '''
    path = '{}/analyses/{}/MutSigNozzleReport2/'.format(data_path, cancer)
    path = path + cancer + '-TP.'
    try:
        rates = pd.read_table(path + 'patients.counts_and_rates.txt')
    except:
        return pd.DataFrame([])
    rates = rates.set_index('name')
    fix_barcode = lambda s: '-'.join(['TCGA'] + s.split('-')[1:3])
    rates.index = rates.index.map(fix_barcode)
    rates = rates[['rate_dbsnp', 'rate_sil', 'rate_non']]
    
    maf = pd.read_table(path + 'final_analysis_set.maf')
    maf = maf.dropna(how='all', axis=[0,1])
    barcode = 'Tumor_Sample_Barcode'  #long Firehose variable names
    maf[barcode] = maf[barcode].map(lambda s: s[:12])
    maf = maf.set_index(['Hugo_Symbol','Tumor_Sample_Barcode'])
    non_silent = maf[maf.is_silent == 0]
    df = pd.DataFrame({pat: df.categ.value_counts() for pat, df in 
                    non_silent.groupby(level=1)}).T
    labels = pd.read_table(path + 'mutation_rates.txt')
    df = df.fillna(0)
    df = df.rename(columns=lambda s: labels.category[s-1])
    pct = (df.T / df.sum(1)).T 
    rates = rates.join(pct)
    rates = rates.rename(columns=lambda s: s.replace('/','_'))
    if patients is not None:
        rates = rates.ix[patients].dropna()
    return rates

def get_cna_rates(data_path, cancer, patients=None):
    '''
    Get copy-number aberration rates from GISTIC processing pipeline.  
    This function depends on the current Firehose output of this program 
    as of July 2013.
    '''
    gistic = FH.get_gistic_gene_matrix(data_path, cancer)
    amp_gene_all = (gistic >= 1).astype(int).sum()
    amp_gene_high = (gistic == 2).astype(int).sum()
    del_gene_all = (gistic <= -1).astype(int).sum()
    del_gene_homo = (gistic <= -2).astype(int).sum()
    
    lesions = FH.get_gistic_lesions(data_path, cancer)
    amp_lesion_all = (lesions.ix['Amplification'] >= 1).sum()
    amp_lesion_high = (lesions.ix['Amplification'] == 2).sum()
    del_lesion_all = (lesions.ix['Deletion'] <= -1).sum()
    del_lesion_homo = (lesions.ix['Deletion'] == -2).sum()
    
    arm_cn = FH.get_gistic_arm_values(data_path, cancer)
    chromosomal_instability = arm_cn.abs().mean()
    
    cna_df = {'gene_amp': amp_gene_all, 'gene_amp_high': amp_gene_high, 
              'gene_del': del_gene_all, 'gene_del_homo': del_gene_homo, 
              'lesion_amp': amp_lesion_all, 'lesion_amp_high': amp_lesion_high, 
              'lesion_del': del_lesion_all, 'lesion_del_homo': del_lesion_homo,
              'chrom_instability': chromosomal_instability}
    cna_df = pd.DataFrame(cna_df)
    if patients is not None:
        cna_df = cna_df.ix[patients].dropna()
    return cna_df

def get_global_vars(data_path, cancer, patients=None):
    '''
    Get compiled DataFrame of global molecular variables from Firehose
    data.  Returns a feature by patient DataFrame with (data-type, variable)
    on the columns and patient barcodes on the index.
    '''
    try:
        data_matrix = FH.read_rnaSeq(data_path, cancer, patients)
        U, S, vH = frame_svd(data_matrix)
        exp_pc = pd.DataFrame({'pc1': vH[0], 'pc2': vH[1]})
    except:
        exp_pc = pd.DataFrame()
        
    try:
        data_matrix = read_methylation(data_path, cancer, patients)
        U, S, vH = frame_svd(data_matrix)
        meth_pc = pd.DataFrame({'pc1': vH[0], 'pc2': vH[1]})
    except:
        meth_pc = pd.DataFrame()
        
    try:
        meth_age, amar = 'FAIL','FAIL'
        #meth_age, amar = get_age_signal(data_path, cancer) 
        meth_pc = meth_pc.join(meth_age).join(amar)
        print 'Should probably check this out'
    except:
        pass
    
    cna_rates = get_cna_rates(data_path, cancer, patients)
    mutation_rates = get_mutation_rates(data_path, cancer, patients)
    
    gv = pd.concat([exp_pc, meth_pc, cna_rates, mutation_rates], 
                    keys=['mRNASeq','methylation', 'cna', 'mutation'], axis=1)
    gv = gv.dropna(how='all', axis=1)
    return gv

def read_data(data_path, cancer, data_type, patients=None, tissue_code='01'):
    fx = {'miRNASeq': FH.read_miRNASeq,
          'mRNASeq': FH.read_rnaSeq,
          'RPPA': FH.read_rppa,
          'Methylation': read_methylation}
    f = fx[data_type]
    return f(data_path, cancer, patients=patients, tissue_code=tissue_code)