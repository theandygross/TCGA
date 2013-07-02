#! /cellar/users/agross/semillon/epd/bin/python

import sys
import pickle as pickle

from Data.Pathways import read_in_pathways
from Processing.GlobalFeatures import RateObject

PATHWAY_FILE = '/cellar/users/agross/Data/GeneSets/c2.cp.v3.0.symbols_edit.csv'
gene_sets, gene_lookup = read_in_pathways(PATHWAY_FILE)

def pickle_cancer_obj(c):
    '''
    Make the cancer object picklable, and drop the file size
    by getting rid of genes with no mutations (note these are
    still stored in the gene_counts attribute). 
    '''
    if hasattr(c, 'hit_matrix'):
        altered = c.hit_matrix.index[c.hit_matrix.sum(1)>0]
        c.hit_matrix = c.hit_matrix.ix[altered]
    if hasattr(c, 'single_matrix'):
        altered = c.single_matrix.index[c.single_matrix.sum(1)>0]
        c.single_matrix = c.single_matrix.ix[altered]

    del c.tests
    pickle.dump(c, open(c.data_path + c.report_folder + 'RateObject.p', 
                        'wb'))
    return c

cancer = sys.argv[1]
data_type = sys.argv[2]
data_path = sys.argv[3]
report_folder = sys.argv[4]

try:
    cancer_obj = RateObject(cancer, data_path, gene_sets, data_type)
    cancer_obj.report_folder = report_folder
    pickle_cancer_obj(cancer_obj)
    print cancer + ' sucess!'

except:
    print cancer + ' fail!'
    sys.exit(-1)

sys.exit(0)
