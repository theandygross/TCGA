'''
Created on Nov 9, 2012

@author: agross
'''
from pandas import read_csv
from numpy import array

HUGO_INFO_PATH = '/cellar/users/agross/Data/GeneSets/'

def get_length(s):
    starts = array(map(int, s['exon_starts'].split(',')[:-1]))
    ends = array(map(int, s['exon_ends'].split(',')[:-1]))
    gene_len = sum(ends - starts)
    return gene_len

def process_gene_lengths(path):
    hugo_info = read_csv(path + 'hugo_info.csv', index_col=0)    
    lengths = hugo_info.apply(get_length, 1)
    lengths = lengths.groupby(level=0).max()
    lengths.to_csv(path + 'coding_lengths.csv')