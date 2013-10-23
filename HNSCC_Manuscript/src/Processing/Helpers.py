'''
Created on Jul 12, 2012

@author: agross
'''
import os as os
import pickle as pickle

import numpy as np
import pandas as pd

import itertools as it

import matplotlib.pyplot as plt
from numpy.linalg import LinAlgError, svd
from numpy import array, diag, sort
from numpy.random import random_integers

from scipy.cluster import hierarchy
from scipy.spatial import distance

from statsmodels.sandbox.stats import multicomp


def transferIndex(source,target):
    return pd.Series(list(target), index=source.index)

def bhCorrection(s, n=None):
    s = s.fillna(1.)
    if n > len(s):
        p_vals = list(s) + [1]*(n-len(s))
    else:
        p_vals = list(s)
    q = multicomp.multipletests(p_vals, method='fdr_bh')[1][:len(s)]
    q = pd.Series(q[:len(s)], s.index, name='p_adj')
    return q
    
def match_series(a,b):
    a, b = a.align(b, join='inner', copy=False)
    valid = pd.notnull(a) & pd.notnull(b)
    a = a[valid].groupby(lambda s: s).first() #some sort of duplicate index bug
    b = b[valid].groupby(lambda s: s).first()
    return a,b

def split_a_by_b(a,b):
    a, b = match_series(a, b)
    groups = [a[b==num] for num in set(b)]
    return groups

def screen_feature(vec, test, df):
    s = pd.DataFrame({f: test(vec, feature) for f,feature in df.iterrows()}).T
    s['q'] = bhCorrection(s.p)
    s = s.sort(columns='p')
    return s

def frame_svd(data_frame, impute='mean'):
    '''
    Wrapper for taking in a pandas DataFrame, preforming SVD
    and outputting the U, S, and vH matricies in DataFrame form.
    '''
    if impute == 'mean':
        data_frame = data_frame.dropna(thresh=int(data_frame.shape[1] * .75))
        data_frame = data_frame.fillna(data_frame.mean())
    
    U,S,vH = svd(data_frame.as_matrix(), full_matrices=False)
    U = pd.DataFrame(U, index=data_frame.index)
    vH = pd.DataFrame(vH, columns=data_frame.columns).T
    return U,S,vH

def extract_pc_old(data_frame, pc_threshold=.2):
    try:
        U,S,vH = frame_svd(((data_frame.T - data_frame.mean(1)) / data_frame.std(1)).T)
    except LinAlgError:
        return None
    p = S**2/sum(S**2)
    return vH[0] if p[0] > pc_threshold else None

def extract_pc(df, pc_threshold=.2, standardize=True):
    if standardize:
        df = ((df.T - df.mean(1)) / df.std(1)).T
    try:
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
    ret = {'pat_vec': pat_vec, 'gene_vec': gene_vec, 'pct_var': pct_var}
    return  ret if pct_var > pc_threshold else None

def df_to_binary_vec(df):
    cutoff = sort(df.sum())[-int(df.sum(1).mean())]
    if (len(df) > 2) and (cutoff == 1.):
        cutoff = 2
    vec = (df.sum() >= cutoff).astype(int)
    return vec

def drop_first_norm_pc(data_frame):
    '''
    Normalize the data_frame by rows and then reconstruct it without the first 
    principal component.  (Idea is to drop the biggest global pattern.)
    '''
    norm = ((data_frame.T - data_frame.mean(1)) / data_frame.std(1)).T
    U,S,vH = frame_svd(norm)
    S[0] = 0   #zero out first pc
    rest = U.dot(pd.DataFrame(diag(S)).dot(vH.T))
    return rest

def cluster_down(df, agg_function, dist_metric='euclidean', num_clusters=50,
                 draw_dendrogram=False):
    '''
    Takes a DataFrame and uses hierarchical clustering to group along the index.
    Then aggregates the data in each group using agg_function to produce a matrix
    of prototypes representing each cluster.          
    '''
    d = distance.pdist(df.as_matrix(), metric=dist_metric)
    D = distance.squareform(d)
    Y = hierarchy.linkage(D, method='complete') 
    c = hierarchy.fcluster(Y, num_clusters, criterion='maxclust')
    c = pd.Series(c, index=df.index, name='cluster')
    clustered = df.join(c).groupby('cluster').aggregate(agg_function)
    if draw_dendrogram:
        fig, ax = plt.subplots(1,1, figsize=(14,2))
        hierarchy.dendrogram(Y, color_threshold=sort(Y[:,2])[-50], no_labels=True, 
                           count_sort='descendent')
        ax.set_frame_on(True)
        ax.set_yticks([])
        return clustered, c, fig
    return clustered, c

def get_random_genes(bp, lengths):
    s = 0
    genes = []
    new_gene = 0
    while s < (bp + new_gene/2.):
        i = random_integers(0, len(lengths)-1)
        genes.append(i)
        new_gene = lengths.ix[i]
        s += new_gene
    genes = lengths.index[genes]
    return genes

def do_perm(f, vec, hit_mat, bp, lengths, iterations):
    real_val = f(vec > 0)
    results = []
    for i in range(iterations):
        perm = hit_mat.ix[get_random_genes(bp, lengths)].sum() > 0
        results.append(f(perm))
    return sum(array(results) <  real_val) / float(len(results))

def run_rate_permutation(df, hit_mat, gene_sets, lengths, f):
    res = {}
    for p,genes in gene_sets.iteritems():
        if p not in df.index:
            continue
        bp = lengths[lengths.index.isin(genes)].sum()
        iterations = 10
        res[p] = do_perm(f, df.ix[p], hit_mat, bp, lengths, iterations)
        while (res[p] <= (10./iterations)) and (iterations <= 2500):
            res[p] = do_perm(f, df.ix[p], hit_mat, bp, lengths, iterations)
            iterations = iterations * 5
    res = sort(pd.Series(res))
    return res

def get_vec_type(vec):
    if vec.count() < 10:
        return
    elif vec.dtype in [float, int]:
        return 'real'
    vc = vec.value_counts()
    if len(vc) == 1 or vc.order()[-2] <= 5:
        return 
    elif len(vc) == 2:
        return 'boolean'
    elif vec.dtype == 'object':
        return 'categorical'   
    
def make_path_dump(obj, file_path):
    dir_path = '/'.join(file_path.split('/')[:-1])
    if not os.path.isdir(dir_path):
        os.makedirs(dir_path)
    pickle.dump(obj, open(file_path, 'wb'))
    
def merge_redundant(df):
    d = df.sort(axis=1).duplicated()
    features = {n: [n] for n,b in d.iteritems() if b == False}
    place=d.index[0]
    for idx,b in d.iteritems():
        if b == True:
            features[place] = features[place] + [idx]
        else:
            place = idx
    features = pd.Series(features)
    
    df = df.ix[d==False]
    df = df.rename(index=features.map(lambda s: '/'.join(s)))
    return df

def add_column_level(tab, arr, name):
    tab = tab.T
    tab[name] = arr
    tab = tab.set_index(name, append=True)
    tab.index = tab.index.swaplevel(0,1)
    return tab.T

def to_quants(vec, q=.25, std=None, labels=False):
    vec = (vec - vec.mean()) / vec.std()
    if q == .5: 
        vec = (vec > 0).astype(int)
        if labels:
            vec = vec.map({0:'Bottom 50%', 1:'Top 50%'})
    elif std is None:
        vec = ((vec > vec.quantile(1-q)).astype(int) - 
               (vec <= vec.quantile(q)).astype(int)).astype(float)
        if labels:
            vec = vec.map({-1:'Bottom {}%'.format(int(q*100)), 0:'Normal', 
                           1:'Top {}%'.format(int(q*100))})
    else:
        vec = (vec - vec.mean()) / vec.std()
        vec = (1.*(vec > std) - 1.*(vec <= (-1*std)))
        if labels:
            vec =  vec.map({-1: 'low', 0: 'normal', 1:'high'})
    return vec

def combine(a,b):
    '''
    Combine two categorical features.
    '''
    combo = (a*1.).add(b*2.)
    combo = combo.dropna()
    if not a.name:
        a.name = 'first'
    if not b.name:
        b.name = 'second'
    if a.name != b.name:
        combo = combo.map({0: 'neither', 1: a.name, 2: b.name, 3:'both'})
    else:
        combo = combo.map({0: 'neither', 1: 'first', 2: 'second', 3:'both'})
    return combo


def true_index(s):
    '''Return indicies for which the variable is true'''
    return s[s].index

def powerset(iterable):
    "http://docs.python.org/2/library/itertools.html#recipes"
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return it.chain.from_iterable(it.combinations(s, r) for r in 
                                  range(len(s)+1))

ti = true_index