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

report_path = sys.argv[1]
cancer_type = sys.argv[2]
data_type = 'MAF'

'''Load in run and mutation data'''
run = pickle.load(open(report_path + '/RunObject.p', 'rb'))
cancer = run.load_cancer(cancer_type)
mut = Dataset(cancer, run, data_type)

'''Create hit_matrix and meta_matrix, filter out genes for features'''
hit_matrix = mut.df.fillna(0).clip_upper(1.)
meta_matrix = DataFrame({p: mut.df.ix[g].sum() for p,g in 
                         run.gene_sets.iteritems()}).T
meta_matrix = meta_matrix.fillna(0).clip_upper(1.)

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

meta_matrix = meta_matrix[size_filter(meta_matrix)] 
s = Series({p: is_one_gene(p) for p in meta_matrix.index})
meta_matrix = meta_matrix.ix[s==False]
hit_matrix = hit_matrix[size_filter(hit_matrix)] 

'''Add passing features to the Data Object''' 
mut.features = meta_matrix.append(hit_matrix)
mut.compress()
mut.uncompress()

'''Save updated Data Object (with additional features field'''
mut.save()
mut.uncompress()

'''Draw pathway_plots for pathway level features'''
meta_features = [f for f in mut.features.index if f in run.gene_sets]
pathway_plot_folder = mut.path + '/Figures/PathwayPlots/'
if not os.path.isdir(pathway_plot_folder):
    os.makedirs(pathway_plot_folder)
        
for i,p in enumerate(meta_features):
    df = mut.df.ix[run.gene_sets[p]]
    pathway_plot(df) 
    savefig(pathway_plot_folder + p)
    
    plt.close()
    plt.gcf().clf()
    if (i%20) == 0:
        gc.collect()

sys.exit()