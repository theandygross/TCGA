'''
Created on Oct 28, 2011

Module to run various enrichments for gene sets obtained by
sub-network find algorithms

@author: agross
'''
import pandas as pandas
from scipy.stats import hypergeom #@UnresolvedImport

def readInPathways(mSigDbFile):
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

def findEnrichedPathways(moduleName, geneSetAll, geneSets, geneLookup, alteredGenes):
    background = list(set(geneLookup.keys()).intersection(alteredGenes))
    geneSet = list(set(geneSetAll).intersection(background))
    allPathways = set()
    pathwayCounts = []
    for gene in geneSet: 
        allPathways = allPathways.union(geneLookup[gene])
                                           
    for pathway in allPathways: 
        pathwayHits = len(geneSets[pathway].intersection(geneSet))
        moduleSize = len(geneSet)
        pathwaySize = len(set(geneSets[pathway]).intersection(alteredGenes))
        pValue = hypergeom.sf(pathwayHits, len(background), pathwaySize, 
                                    moduleSize)
        if(pValue < .00000001 and pathwayHits > 2):
            pathwayCounts.append((moduleName, pathway, pathwayHits, moduleSize, 
                                  len(background), pathwaySize, pValue))
    return sorted(pathwayCounts, key=lambda tup: tup[5], reverse=False)

def getEnrichments(modules, geneSets, geneLookup): 
    '''
    Function to get enrichments for a given network-module
    If output file is not given, just returns a string of enrichments
    
    inputs:
      modules:         list of ModuleNode objects
      geneSets:        dictionary mapping gene group names to genes in group
      geneLookup:      dictionary mapping gene to gene sets annotated for that gene
      output:          (optional) file to output enrichments
      
    output:
      string of gene set enrichments for each module
    ''' 
    pathwayCounts = []   
    for i,mod in enumerate(modules):
        background = list(set(geneLookup.keys()).intersection(mod.geneList))
        geneSet = list(set(mod.genes).intersection(background))
        try:
            allPathways = reduce(set.union, [geneLookup[gene] for gene in mod.genes 
                                             if gene in geneLookup])
        except:
            allPathways = []
        for pathway in allPathways: 
            pathwayHits = len(geneSets[pathway].intersection(geneSet))
            moduleSize = len(geneSet)
            pathwaySize = len(set(geneSets[pathway]).intersection(mod.geneList))
            pValue = hypergeom.sf(pathwayHits, len(background), pathwaySize, 
                                    moduleSize)
            if(pValue < .00000001 and pathwayHits > 2):
                pathwayCounts.append((mod.label, i, pathway, pathwayHits, moduleSize, 
                                      len(background), pathwaySize, pValue))
    dataFrame = pandas.DataFrame(pathwayCounts, 
                                  columns=['data type', 'module #', 'pathway','pathwayHits','moduleSize','backgroundSize', 'pathwaySize','pValue'])
        
    return dataFrame

def geneSetsToMetaMatrix(geneSets, geneMatrix, minGroupSize=4, setFilter=None):
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