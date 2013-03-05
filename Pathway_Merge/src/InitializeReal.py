#! /usr/bin/env python

import os as os
import sys as sys
import pickle as pickle

from pandas import DataFrame
from matplotlib.pylab import subplots
from numpy import diag

from Data.Containers import Dataset
from Processing.Helpers import frame_svd, extract_pc

report_path = sys.argv[1]
cancer_type = sys.argv[2]
data_type = sys.argv[3]

'''Load in run and mutation data'''
run = pickle.load(open(report_path + '/RunObject.p', 'rb'))
cancer = run.load_cancer(cancer_type)
data = Dataset(cancer, run, data_type)
if data_type == 'mRNASeq':
    data.df = data.df.groupby(by=lambda n: n.split('|')[0]).mean()

'''Normalize data and calculate principal components'''
df = data.df
norm = ((df.T - df.mean(1)) / df.std(1)).T
U,S,vH = frame_svd(norm)
data.pc1, data.pc2 = vH[0], vH[1]
data.loading1, data.loading2 = U[0], U[1]

'''Reconstruct data without first pc'''
S_n = S.copy()
S_n[0] = 0
rest = U.dot(DataFrame(diag(S_n)).dot(vH.T))

'''Extract PCs for all gene sets'''
pc = {p: extract_pc(rest.ix[g]) for p,g, in run.gene_sets.iteritems()}
pc = DataFrame({p: s for p,s in pc.iteritems() if s}).T

data.pc_genes = pc.gene_vec
data.pct_var = pc.pct_var
data.features = DataFrame(pc.pat_vec.to_dict()).T

'''Create figures for all gene sets'''
pathway_plot_folder = data.path + '/Figures/PathwayPlots/'
if not os.path.isdir(pathway_plot_folder):
    os.makedirs(pathway_plot_folder)
    
for feature in data.features.index:
    fig, axs = subplots(1,2, figsize=(10,3))
    data.pc_genes.ix[feature].order().plot(kind='bar', ax=axs[0])
    axs[0].annotate('Explained Variation: %.2f' % data.pct_var[feature], 
                    (.03,.97), xycoords='axes fraction', ha='left', va='top')
    axs[0].set_ylabel('Eigen-Patient Loading')
    data.features.ix[feature].hist(ax=axs[1])
    axs[1].set_xlabel('Eigen-Gene Loading')
    fig.tight_layout()
    fig.savefig(pathway_plot_folder + feature)
    
data.save()