'''
Created on Jun 23, 2013

@author: agross
'''
import os as os
import numpy as np
import pandas as pd

GENE_POS = pd.read_csv('/cellar/users/agross/Data/GeneSets/HGNC_chr_pos.txt')
GENE_LIST_FILE = '/cellar/users/agross/Data/GeneSets/HGNC_Genes'
GENES = open(GENE_LIST_FILE, 'rb').read().split('\n')
COM = {'T':'A', 'A':'T', 'G':'C', 'C':'G'}
BASE_CHANGES = ['C>T', 'C>A', 'C>G', 'A>G', 'A>T', 'A>C']
COLUMNS = ['reference_allele', 'tumor_seq_allele1', 'tumor_seq_allele2',
           'entrez_gene_id', 'chromosome', 'start_position', 'variant_type',
           'variant_classification']

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
    f = open(mSigDbFile, 'r')
    geneSets = {}
    genes = []
    for line in f:
        line = line.replace('\"', '')
        tmp = line.strip().split("\t")
        setName = tmp[0]
        geneSets[setName] = set(tmp[2:])
        genes.extend(tmp[2:])
    f.close()
    genes = set(genes)

    geneLookup = dict([(gene, set()) for gene in genes])
    for pathway in geneSets: 
        for gene in geneSets[pathway]: 
            geneLookup[gene].add(pathway)
    return geneSets, geneLookup

def map_splice_junction_to_gene(j):
    '''
    Takes a label from the Firehose junction_quantification file and maps to 
    a gene. 
    Needs to be refactored for speed. 
    '''
    parsed = j.split(':')
    chromosome = parsed[0][3:]
    start = int(parsed[1])
    end = int(parsed[3])
    strand = parsed[4]
    
    lu = GENE_POS[((GENE_POS['Chromosome Name'] == chromosome) * 
                   (GENE_POS['Gene Start (bp)'] < start) * 
                   (GENE_POS['Gene End (bp)'] > end))
                  ].sort('Gene Start (bp)')
    if len(lu) == 0:
        return '?'
    return lu.iloc[0].ix['HGNC symbol']

"""
MutationBaseChange
-------------------
"""
def good_mut(s):
    return ((s['variant_type'] in ['SNP']) & 
            (s['variant_classification'] not in ['Silent', 'RNA']) & 
            (s.name in GENES))
    
def get_base_change(s):
    s = s[['reference_allele', 'tumor_seq_allele1', 'tumor_seq_allele2']]
    if s[0] not in list('ATCG'):
        return np.nan
    cols = s[[0, 1]] if s[0] != s[1] else s[[0, 2]]
    return '>'.join(cols)

def get_base_change_matrix(mutation_dir):
    mutation_matrix = pd.DataFrame()
    for f in os.listdir(mutation_dir):
        if f[:4] != 'TCGA':
            continue
        tab = pd.read_table(mutation_dir + f, index_col=0)
        tab = tab.select(lambda c: c.lower() in COLUMNS, axis=1)
        tab = tab.rename(columns=str.lower)
        id_cols = ['entrez_gene_id', 'chromosome', 'start_position']
        tab = tab.drop_duplicates(cols=id_cols)
        tab = tab[tab.apply(good_mut, axis=1)]
        tab = tab.apply(get_base_change, axis=1)
        tab = tab.dropna().groupby(level=0).aggregate(lambda s: ','.join(s))
        tab.name = f[:12]
        mutation_matrix = mutation_matrix.join(tab, how='outer')
    return mutation_matrix
    
rev = lambda s: ''.join([COM[s[0]], s[1], COM[s[2]]])
def count(vec, change):
    vec = vec.dropna().apply(lambda s: (change in s) or (rev(change) in s))
    return vec.sum()
count_nucs = lambda r: pd.Series(dict((change, count(r, change)) for change in 
                                   BASE_CHANGES))
def get_base_change_counts(mutation_dir=None, base_change_matrix=None):
    '''
    Takes in a base change matrix and counts the changes (including RC).
    '''
    if mutation_dir is not None:
        base_change_matrix = get_base_change_matrix(mutation_dir)
    else:
        assert base_change_matrix is not None
    return base_change_matrix.apply(count_nucs)
