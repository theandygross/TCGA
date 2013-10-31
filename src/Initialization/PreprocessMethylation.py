'''
Created on Jun 22, 2013

@author: agross
'''
import os as os
import pandas as pd
import numpy as np

from subprocess import call, check_output
from Processing.Helpers import extract_pc, true_index
from Data.Firehose import fix_barcode_columns
from Data.Intermediate import get_beta_values

def extract_betas(folder, meth_file, probes='All', outfile='beta_values.txt'):
    '''
    TCGA has a lot of unnecessary columns that eat up memory, so I parse out 
    the beta-values and save them to an intermediate file. 
    I do this with system calls because of the file size.
    '''
    curdir = os.getcwd()
    os.chdir(folder)
    
    if probes != 'All':
        target_file = 'picked.txt'
        query = '\\|'.join(['Hybridization REF'] + list(probes))
        call(['grep', query, meth_file],
             stdout=open(target_file, 'wb'))
    else:
        target_file = meth_file
    cols = check_output(['awk', '-F\t', '{print NF; exit}', target_file])
    cols = int(cols.strip())
    cols = ','.join(map(str, range(2, int(cols) + 1, 4)))
    
    call(['cut', '-f' + cols, target_file], stdout=open('tmp', 'wb'))
    call(['cut', '-f1,3', target_file], stdout=open('idx', 'wb'))
    call(['paste', 'idx', 'tmp', '-d\t'], stdout=open(outfile, 'wb'))
    call(['rm', 'tmp', 'idx'])
    
    os.chdir(curdir)

def peel_pc(df):
    '''
    Wrapper around extract_pc.
    Flips the PC slightly differently based on correlation with the mean. 
    Does not standardize data for PCA due to underlying distribution of
    beta values.
    '''
    try:
        r = extract_pc(df - .5)
        l, r, p = r['gene_vec'], r['pat_vec'], r['pct_var']
        mean = df.mean(1)
        if l.corr(mean) < 0:
            l = l * -1
            r = r * -1
        return l, r, p
    except:
        r = df.mean()
        return np.nan, r, np.nan

def process_meth(data_path, cancer, probeset='All'):
    '''
    Gets data into a more usable format for our analysis. 
    Processes large probe by patient matrix and extracts principal components
    for each gene's constituent probes. 
    Saves the data in the ucsd_processing folder.
    '''
    path = '{}stddata/{}/methylation'.format(data_path, cancer)
    meth_path = [f[0] for f in os.walk(path)
                 if 'humanmethylation450' in f[0]][-1]
    meth_file = [f for f in os.listdir(meth_path) if 'data' in f][0]
    
    outpath = '{}/ucsd_processing/{}/methylation450/'.format(data_path, cancer)
    if not os.path.isdir(outpath):
        os.makedirs(outpath)
    if os.path.isfile(outpath + 'meta_probes.csv'):  # file already exists
        return
    outfile = outpath + ('beta_values.txt' if probeset == 'All' else 
                         'beta_values_picked.txt')
    
    if not os.path.isfile(outfile):
        extract_betas(meth_path, meth_file, probeset, outfile) 
    
    table = get_beta_values(data_path, cancer)
    
    loadings = {}
    meta_probes = {}
    pct_var = {}
    for g in table.index.levels[0]:
        loadings[g], meta_probes[g], pct_var[g] = peel_pc(table.ix[g])
    
    loadings = {k:v for k, v in loadings.iteritems() if type(v) == pd.Series}
    loadings = pd.concat(loadings.values(), keys=loadings.keys())
    meta_probes = pd.DataFrame(meta_probes).T
    pct_var = pd.Series(pct_var)
    
    loadings.to_csv(outpath + 'loadings.csv')
    meta_probes.to_csv(outpath + 'meta_probes.csv')
    pct_var.to_csv(outpath + 'pct_var.csv')
