'''
Created on Sep 25, 2012

@author: agross
'''
from pandas import Series, DataFrame, read_table
from numpy import nan
import os as os

GENE_LIST_FILE = '/cellar/users/agross/Data/GeneSets/HGNC_Genes'
GENES = open(GENE_LIST_FILE, 'rb').read().split('\n')
COM = {'T':'A','A':'T','G':'C','C':'G'}
BASE_CHANGES = ['C>T','C>A','C>G','A>G','A>T','A>C']
COLUMNS = ['reference_allele','tumor_seq_allele1', 'tumor_seq_allele2',
           'entrez_gene_id', 'chromosome','start_position', 'variant_type',
           'variant_classification']

good_mut = lambda s: ((s['variant_type'] in ['SNP']) & 
                     (s['variant_classification'] not in ['Silent', 'RNA']) & 
                     (s.name in GENES))
    
    
def get_base_change(s):
    s = s[['reference_allele','tumor_seq_allele1', 'tumor_seq_allele2']]
    if s[0] not in list('ATCG'):
        return nan
    cols = s[[0,1]] if s[0] != s[1] else s[[0,2]]
    return '>'.join(cols)

def get_base_change_matrix(mutation_dir):
    mutation_matrix = DataFrame()
    for f in os.listdir(mutation_dir):
        if f[:4] != 'TCGA':
            continue
        tab = read_table(mutation_dir + f, index_col=0)
        tab = tab.select(lambda c: c.lower() in COLUMNS, axis=1)
        tab = tab.rename(columns=str.lower)
        tab = tab.drop_duplicates(cols=['entrez_gene_id', 'chromosome','start_position'])
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
count_nucs = lambda r: Series(dict((change, count(r, change)) for change in 
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