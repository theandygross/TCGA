'''
Created on Oct 28, 2011

Module to run various enrichments for gene sets obtained by
sub-network find algorithms

@author: agross
'''
import pandas as pandas
from scipy.stats import hypergeom #@UnresolvedImport

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


def build_meta_matrix(geneSets, geneMatrix, minGroupSize=4, setFilter=None):
    if not isinstance(geneSets, dict):
        geneSets = dict(enumerate(geneSets))
    metaMatrix = dict((group, geneMatrix.ix[genes].sum(0)) for group,genes in
                      geneSets.items())
    metaMatrix = pandas.DataFrame(metaMatrix)
    if setFilter is not None: 
        metaMatrix = metaMatrix.apply(setFilter, axis=0)
    sizeFilter = lambda s: s.sum() in range(minGroupSize+1,len(metaMatrix)-minGroupSize)
    metaMatrix = metaMatrix.T[metaMatrix.apply(sizeFilter)]
    return metaMatrix