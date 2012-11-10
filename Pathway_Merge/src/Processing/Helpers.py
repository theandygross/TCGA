'''
Created on Jul 12, 2012

@author: agross
'''
import matplotlib.pyplot as plt
from pandas import Series, DataFrame, notnull
from numpy.linalg import LinAlgError
from numpy import array, diag, sort
from numpy.random import random_integers

from scipy.cluster import hierarchy
from scipy.spatial import distance
from rpy2.robjects import r, FloatVector

transferIndex = lambda source,target: Series(list(target), index=source.index)

def bhCorrection(s): 
    return Series(r('p.adjust')(FloatVector(s), method='BH'),index=s.index, 
                  name=s.name)
    
def match_series(a,b):
    a, b = a.align(b, join='inner', copy=False)
    valid = notnull(a) & notnull(b)
    a = a[valid]
    b = b[valid]
    return a,b

def split_a_by_b(a,b):
    a, b = match_series(a, b)
    groups = [a[b==num] for num in set(b)]
    return groups

def frame_svd(data_frame):
    '''
    Wrapper for taking in a pandas DataFrame, preforming SVD
    and outputting the U, S, and vH matricies in DataFrame form.
    '''
    from numpy.linalg import svd #@UnresolvedImport
    from pandas import DataFrame
    
    U,S,vH = svd(data_frame.as_matrix(), full_matrices=False)
    U = DataFrame(U, index=data_frame.index)
    vH = DataFrame(vH, columns=data_frame.columns).T
    return U,S,vH

def extract_pc(data_frame, pc_threshold=.2):
    try:
        U,S,vH = frame_svd(((data_frame.T - data_frame.mean(1)) / data_frame.std(1)).T)
    except LinAlgError:
        return None
    p = S**2/sum(S**2)
    return vH[0] if p[0] > pc_threshold else None

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
    U,S,vH = frame_svd(norm.dropna())
    S[0] = 0   #zero out first pc
    rest = U.dot(DataFrame(diag(S)).dot(vH.T))
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
    c = Series(c, index=df.index, name='cluster')
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
    res = sort(Series(res))
    return res