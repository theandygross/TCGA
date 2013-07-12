#! /usr/bin/env python

import os as os
import sys as sys
import pickle as pickle
import gc as gc

from pandas import Series, DataFrame
from matplotlib.pylab import savefig
import matplotlib.pyplot as plt

from Data.Containers import Dataset
from Reports.Figures import pathway_plot
from Processing.Helpers import merge_redundant



report_path = sys.argv[1]
cancer_type = sys.argv[2]
data_type = sys.argv[3]
data_type = data_type[3:]

'''Load in run and CN data'''
run = pickle.load(open(report_path + '/RunObject.p', 'rb'))
cancer = run.load_cancer(cancer_type)

if data_type == 'broad':
    data = Dataset(cancer, run, 'CN_broad')
    data.features = data.df
    data.save()
    sys.exit(0)
    
data = Dataset(cancer, run, 'CN')
data.path = '_'.join([data.path, data_type])

if data_type == 'deletion':
    data.hit_val = -2
elif data_type == 'amplification':
    data.hit_val = 2
elif data_type == 'amplification_low':
    data.df = data.df.replace(1,2)
    data.hit_val = 2

hit_matrix = (data.df==data.hit_val).astype(float)

genes_in_bands = hit_matrix.groupby(level=0).size()
frac = (hit_matrix.groupby(level=0).sum().T / (genes_in_bands + 1)).T
hit_matrix = hit_matrix.mul((frac <= .25).astype(float), level=0, 
                            fill_value=0)
hit_matrix = hit_matrix[hit_matrix.sum(1) > 0]
hit_matrix.index = hit_matrix.index.reorder_levels([2,1,0])

def size_filter(s):
    '''Make sure features covers a minimum number of patients'''
    min_p = run.parameters['min_patients']
    return s.sum(1).isin(range(min_p, meta_matrix.shape[1] - min_p))

def is_one_gene(p):
    '''Test to see if most mutations are due to single gene'''
    counts = hit_matrix.ix[run.gene_sets[p]].sum(1).dropna().order()
    with_top = hit_matrix.ix[run.gene_sets[p]].sum().clip_upper(1).sum()
    without = hit_matrix.ix[run.gene_sets[p] - 
                            {counts.idxmax()}].sum().clip_upper(1).sum()
    return ((with_top - without) / without) > .5

meta_matrix = DataFrame({p: hit_matrix.ix[g].sum() for p,g in 
                         run.gene_sets.iteritems()}).T
meta_matrix = meta_matrix.fillna(0).clip_upper(1.)
min_p = run.parameters['min_patients']

meta_matrix = meta_matrix[size_filter(meta_matrix)]
s = Series({p: is_one_gene(p) for p in meta_matrix.index})
meta_matrix = meta_matrix.ix[s==False]
hit_matrix = hit_matrix[size_filter(hit_matrix)]

hit_genes = hit_matrix.copy()
hit_genes.index = hit_genes.index.get_level_values(0)
non_redundant = merge_redundant(hit_genes.append(meta_matrix))

'''Add passing features to the Data Object'''
data.features = non_redundant

'''Save updated Data Object (with additional features field)'''
data.df = data.df.replace([1,-1], 0)
data.save()
data.uncompress()


'''Draw pathway_plots for pathway level features'''
meta_features = [f for f in data.features.index if f in run.gene_sets]
pathway_plot_folder = data.path + '/Figures/PathwayPlots/'
if not os.path.isdir(pathway_plot_folder):
    os.makedirs(pathway_plot_folder)
    
hit_mat = data.df.copy()
hit_mat.index = hit_mat.index.get_level_values(2)
hit_mat = (hit_mat == data.hit_val).astype(float)
        
for i,p in enumerate(meta_features):
    if os.path.isfile(pathway_plot_folder + p):
        continue
    df = hit_mat.ix[run.gene_sets[p]]
    pathway_plot(df) 
    savefig(pathway_plot_folder + p)
    
    plt.close()
    plt.gcf().clf()
    if (i%20) == 0:
        gc.collect()
