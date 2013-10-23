'''
Created on Jul 2, 2013

@author: agross
'''

import os as os
import pickle as pickle
import gc as gc
from Data.Containers import Dataset
import Data.Firehose as FH

from pandas import Series, DataFrame
from matplotlib.pylab import savefig
import matplotlib.pyplot as plt
from Figures.MemoPlot import pathway_plot

def is_one_gene(genes, df):
    '''Test to see if most mutations are due to single gene'''
    top_hit = df.ix[genes].sum(1).idxmax()
    with_top = df.ix[genes].sum().clip_upper(1).sum()
    without = df.ix[genes - {top_hit}].sum().clip_upper(1).sum()
    return ((with_top - without) / (without + .1)) > .5

def size_filter(s, min_pat=10):
    '''Test if all sufficient feature diversity'''
    vc = s.clip_upper(1.).value_counts()
    return (len(vc)== 2) and (vc.min() >= min_pat)


class MutDataset(Dataset):
    '''
    Inherits from Dataset class.  Adds some added processing for mutation
    data.
    '''
    def __init__(self, run, cancer, patients=None, 
                 create_features=True, draw_figures=False):
        '''
        '''
        Dataset.__init__(self, cancer.path, 'Mutation', compressed=False)
        self.df = FH.get_mutation_matrix(run.data_path, cancer.name)
        if patients is not None:
            self.df = self.df.ix[:, patients].dropna(1, how='all')
        
        if create_features is True:
            min_pat = run.parameters['min_patients']
            self._create_feature_matrix(run.gene_sets, min_pat)
            
        if draw_figures is True:
            self._create_pathway_figures(run.gene_sets)
            
    def _create_feature_matrix(self, gene_sets, min_size=10):
        '''
        Create hit_matrix and meta_matrix, filter out genes for features.
        '''
        hit_matrix = self.df.fillna(0).clip_upper(1.)
        meta_matrix = DataFrame({p: self.df.ix[g].sum() for p,g in 
                                 gene_sets.iteritems()}).T
        meta_matrix = meta_matrix.fillna(0).clip_upper(1.)
        meta_matrix = meta_matrix.dropna()
        s = meta_matrix.apply(size_filter, args=(min_size,), axis=1)
        meta_matrix = meta_matrix.ix[s]
        #s = Series({p: is_one_gene(gene_sets[p], hit_matrix) for p in 
        #            meta_matrix.index})
        s = [p for p in meta_matrix.index 
             if is_one_gene(gene_sets[p], hit_matrix) == False]
        meta_matrix = meta_matrix.ix[s]
        s = hit_matrix.apply(size_filter, args=(min_size,), axis=1)
        hit_matrix = hit_matrix.ix[s]
        self.features = meta_matrix.append(hit_matrix)
        self.compress()
        self.uncompress()
        
    def _create_pathway_figures(self, gene_sets):       
        '''
        Draw pathway_plots for pathway level features.
        '''
        meta_features = [f for f in self.features.index if f in gene_sets]
        pathway_plot_folder = self.path + '/Figures/PathwayPlots/'
        if not os.path.isdir(pathway_plot_folder):
            os.makedirs(pathway_plot_folder)
                
        for i,p in enumerate(meta_features):
            df = self.df.ix[gene_sets[p]]
            pathway_plot(df) 
            savefig(pathway_plot_folder + p)
            
            plt.close()
            plt.gcf().clf()
            if (i%20) == 0:
                gc.collect()
                
def initialize_mut(cancer_type, report_path, patients=None, 
                   create_meta_features=True, draw_figures=False, save=True):
    '''
    Initialize mutation data for down-stream analysis.
    '''
    run = pickle.load(open(report_path + '/RunObject.p', 'rb'))
    cancer = run.load_cancer(cancer_type)
    
    data = MutDataset(run, cancer, patients, create_meta_features, 
                      draw_figures)
        
    if save is True:
        data.save()
    
    data.uncompress()
    return data