
import os as os
import sys as sys
import pickle as pickle
from pandas import Series, DataFrame
import pandas as pd
from Data.Methylation import *

METH_FOLDER = ('methylation/humanmethylation450/jhu_usc_edu/Level_3/' + 
               'within_bioassay_data_set_function/data/')
METH_FILE = 'data.txt'

result_path = '/scratch/TCGA/Firehose__2012_01_16/ucsd_analyses'
run = sorted(os.listdir(result_path))[0]
run = pickle.load(open('/'.join([result_path, run, 'RunObject.p']), 'rb'))


def pull_out_beta_values(folder, probes='All', outfile='beta_values.txt'):
    
    curdir = os.getcwd()
    os.chdir(folder)
    
    if probes != 'All':
        target_file = 'picked.txt'
        query = '\\|'.join(['Hybridization REF'] + list(probes))
        call(['grep', query, METH_FILE], 
             stdout=open(target_file,'wb'))
    else:
        target_file = METH_FILE
    cols = check_output(['awk', '-F\t', '{print NF; exit}', target_file])
    cols = int(cols.strip())
    cols = ','.join(map(str, range(2, int(cols)+1,4)))
    
    call(['cut', '-f' + cols, target_file], stdout=open('tmp','wb'))
    call(['cut', '-f1,3', target_file], stdout=open('idx','wb'))
    call(['paste', 'idx', 'tmp', '-d\t'], stdout=open(outfile,'wb'))
    call(['rm', 'tmp', 'idx'])
    
    os.chdir(curdir)
    
from Processing.Helpers import frame_svd
def extract_pc(df, pc_threshold=.2, standardize=True):
    try:
        if standardize:
            df = ((df.T - df.mean(1)) / df.std(1)).T
        U,S,vH = frame_svd(df)
    except LinAlgError:
        return None
    p = S**2/sum(S**2)
    pat_vec = vH[0]
    gene_vec = U[0]
    pct_var = p[0]
    if sum(gene_vec) < 0:
        gene_vec = -1*gene_vec
        pat_vec = -1*pat_vec
    #pat_vec = (pat_vec - pat_vec.mean()) / pat_vec.std()
    ret = {'pat_vec': pat_vec, 'gene_vec': gene_vec, 'pct_var': pct_var}
    return  ret if pct_var > pc_threshold else None

def peel_pc(vec):
    r = extract_pc(vec-.5, standardize=False)
    l,r,p = r['gene_vec'], r['pat_vec'], r['pct_var']
    mean = vec.mean(1)
    if l.corr(mean) < 0:
        l = l*-1
        r = r*-1
    return l,r,p

def process_meth(cancer, data_path, probeset = 'All'):
    folder = data_path + '/'.join(['stddata', cancer, METH_FOLDER])
    
    outpath = data_path + '/'.join(['ucsd_processing', cancer, 
                                    'methylation450', ''])
    if not os.path.isdir(outpath):
        os.makedirs(outpath)
    outfile = outpath + ('beta_values.txt' if probeset == 'All' else 
                         'beta_values_picked.txt')
    
    pull_out_beta_values(folder, probeset, outfile)
    
    table = read_table(outpath + 'beta_values.txt', skiprows=[1], index_col=[0])
    table = table.rename(columns=lambda s: s if s!=table.columns[0] else 'symbol')
    t = table.set_index('symbol', append=True)
    t = t.swaplevel(0,1)
    t = t.select(lambda s: s.split('-')[3].startswith('01'), axis=1) 
    
    loadings = {}
    meta_probes = {}
    pct_var = {}
    for gene in t.index.levels[0]:
        try:
            loadings[gene], meta_probes[gene], pct_var[gene] = peel_pc(t.ix[gene])
        except:
            meta_probes[gene] = t.ix[gene].mean()
    loadings = pd.concat(loadings.values(), keys=loadings.keys())
    meta_probes = DataFrame(meta_probes).T
    pct_var = Series(pct_var)
    
    loadings.to_csv(outpath + 'loadings.csv')
    meta_probes.to_csv(outpath + 'meta_probes.csv')
    pct_var.to_csv(outpath + 'pct_var.csv')
    
    process_meth(sys.argv[1], run.data_path)