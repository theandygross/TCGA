#! /cellar/users/agross/semillon/epd/bin/python

import sys
import pickle as pickle
from collections import defaultdict

from pandas import DataFrame

from Reports.Reports import create_clinical_report

clinical_obj = sys.argv[1]
next_cancer = sys.argv[2]
prev_cancer = sys.argv[3]

def get_gene_lookup(gene_sets):
    gene_lookup = defaultdict(set)
    for pathway, gene_set in gene_sets.iteritems():
        for gene in gene_set:
            gene_lookup[gene].add(pathway)
    return gene_lookup

cancer = pickle.load(open(clinical_obj, 'rb'))

if hasattr(cancer, 'p_genes'):
    if hasattr(cancer, 'filter_bad_pathways'):
        gene_lookup = get_gene_lookup(cancer.gene_sets)
        cancer.filter_bad_pathways(gene_lookup)
        cancer.q_pathways = cancer.q_pathways.fillna(1.)
else:
    cancer.q_genes = DataFrame(columns=cancer.q_pathways.columns)
    cancer.p_genes = DataFrame(columns=cancer.q_pathways.columns)
    
if hasattr(cancer, '_res'):
    cancer.q_pathways.rate = cancer._res    
'''
try:
    create_clinical_report(cancer, next_cancer, prev_cancer)
    print cancer.cancer + ' sucess!'    
except:
    print cancer.cancer + ' fail!'
    sys.exit(-1)
'''   
create_clinical_report(cancer, next_cancer, prev_cancer)
sys.exit(0)