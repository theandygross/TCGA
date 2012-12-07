'''
Created on Oct 28, 2011

Module to run various enrichments for gene sets obtained by
sub-network find algorithms

@author: agross
'''
from pandas import DataFrame

def read_in_pathways(mSigDbFile):
    '''
    Reads in mSigDb pathways.
    File format should be that downloaded from website, tsv with format:
    pathway \t gene1 \t gene2 ... geneX \n
    input:
      mSigDbFile:        pathway file, can be downloaded from the web-site
    output:
      geneSets:          dict mapping pathway name to set of genes in the pathway
      geneLookup:        dict mapping genes to the pathways they occur in
    '''    
    file = open(mSigDbFile, 'r')

    geneSets = {}
    genes = []
    for line in file:
        line = line.replace('\"','')
        tmp = line.strip().split("\t")
        setName  = tmp[0]
        geneSets[setName]= set(tmp[2:])
        genes.extend(tmp[2:])
    genes = set(genes)

    geneLookup = dict([(gene, set()) for gene in genes])
    for pathway in geneSets: 
        for gene in geneSets[pathway]: 
            geneLookup[gene].add(pathway)
    return geneSets, geneLookup


def build_meta_matrix(gene_sets, gene_matrix, min_size=4, set_filter=None):
    if not isinstance(gene_sets, dict):
        gene_sets = dict(list(enumerate(gene_sets)))
    meta_matrix = DataFrame({group: gene_matrix.ix[genes].sum(0) 
                             for group,genes in gene_sets.items()})
    if set_filter is not None: 
        meta_matrix = meta_matrix.apply(set_filter, axis=0)
    keepers = meta_matrix.sum().isin(range(min_size, 
                                           len(meta_matrix) - min_size))
    meta_matrix = meta_matrix.ix[:,keepers].T
    return meta_matrix