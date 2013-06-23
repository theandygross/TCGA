
import os as os
import sys as sys
from pandas import Series, DataFrame
import pandas as pd
import numpy as np
from subprocess import call, check_output

def pull_out_beta_values(folder, meth_file, probes='All', outfile='beta_values.txt'):
    
    curdir = os.getcwd()
    os.chdir(folder)
    
    if probes != 'All':
        target_file = 'picked.txt'
        query = '\\|'.join(['Hybridization REF'] + list(probes))
        call(['grep', query, meth_file], 
             stdout=open(target_file,'wb'))
    else:
        target_file = meth_file
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
    except np.linalg.LinAlgError:
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
    meth_path = [f[0] for f 
                 in os.walk('{}stddata/{}/methylation'.format(data_path, cancer))
                 if 'humanmethylation450' in f[0]][-1]
    meth_file = [f for f in os.listdir(meth_path) if 'data' in f][0]
    
    outpath = data_path + '/'.join(['ucsd_processing', cancer, 
                                    'methylation450', ''])
    if not os.path.isdir(outpath):
        os.makedirs(outpath)
    outfile = outpath + ('beta_values.txt' if probeset == 'All' else 
                         'beta_values_picked.txt')
    
    if os.path.isfile(outpath + 'meta_probes.csv'):
        return
    
    if not os.path.isfile(outpath + outfile):
        pull_out_beta_values(meth_path, meth_file, probeset, outfile) 
    
    table = pd.read_table(outpath + 'beta_values.txt', skiprows=[1], index_col=[0])
    table = table.rename(columns=lambda s: s if s!=table.columns[0] else 'symbol')
    t = table.set_index('symbol', append=True)
    t = t.swaplevel(0,1)
    t.columns = pd.MultiIndex.from_tuples([(i[:12], i[13:15]) for i in t.columns])
    
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

if __name__ == "__main__":
    cancer = sys.argv[0]
    data_path = sys.argv[1]
    process_meth(cancer, data_path)
    sys.exit(0)