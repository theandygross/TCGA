#! /cellar/users/agross/semillon/epd/bin/python

import sys
import pickle as pickle

from Data.Pathways import read_in_pathways
from Processing.Clinical import ClinicalObject

PATHWAY_FILE = '/cellar/users/agross/Data/GeneSets/c2.cp.v3.0.symbols_edit.csv'
gene_sets, gene_lookup = read_in_pathways(PATHWAY_FILE)

REAL_VARIABLES = ['AMAR','rate','age','karnofsky_performance', 'pack_years',
                  'pct_tumor_invasion', 'height', 'weight', 'neo_area']
BINARY_VARIABLES = ['gender', 'therapy', 'radiation', 'triple_neg', 'triple_pos',
                    'ER_pos','PR_pos','her2_pos', 'lymphnode_n0n1', 'tumor_t1t2',
                    'post_menopause', 'histo_g1g2', 'neo_status', 'chemo', 
                    'hormones''complete_response', 'metastatic_recurrence',
                    'new_tumor','smoker', 'drinker', 'aml_cyto_risk_favorable',
                    'morphology_m1m2', 'normal_cyto', 'abnormallymphocyte', 
                    'bonemarrowbandcell', 'bonemarrowbasophil', 'bonemarrowblastcell',
                    'bonemarrowcellularity', 'bonemarrowlabeosinophil', 
                    'bonemarrowlymphocyte', 'bonemarrowmyelocyte', 'bonemarrowneutrophil', 
                    'bonemarrowprolymphocyte', 'bonemarrowpromonocytecount', 
                    'bonemarrowpromyelocyte', 'monocyte', 'venous_invasion',
                    'lymphatic_invasion', 'calcium_level', 'white_cell_count',
                    'tumor_focality']
SURVIVAL_TESTS = {'survival_3y' : {'event_var' : 'deceased_3y', 'time_var' : 'days_3y', 
                                   'covariates' : ['age', 'rate']},
                  'survival_5y' : {'event_var' : 'deceased_5y', 'time_var' : 'days_5y', 
                                   'covariates' : ['age', 'rate']},
                  'event_free_survival_3y' : {'event_var' : 'event_3y', 
                                           'time_var' : 'event_free_survival_3y', 
                                           'covariates' : ['age', 'rate']},
                  'event_free_survival_5y' : {'event_var' : 'event_5y', 
                                           'time_var' : 'event_free_survival_5y', 
                                           'covariates' : ['age', 'rate']}
                  }

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
    pickle.dump(c, open(c.data_path + c.report_folder + 'ClinicalObject.p', 
                        'wb'))
    return c

cancer = sys.argv[1]
data_type = sys.argv[2]
data_path = sys.argv[3]
report_folder = sys.argv[4]
drop_pc = True

try:
    cancer_obj = ClinicalObject(cancer, data_path, gene_sets, data_type,
                                survival_tests=SURVIVAL_TESTS, 
                                real_variables=REAL_VARIABLES,
                                binary_variables=BINARY_VARIABLES,
                                drop_pc=drop_pc)
    cancer_obj.report_folder = report_folder
    pickle_cancer_obj(cancer_obj)
    print cancer + ' sucess!'

except:
    print cancer + ' fail!'
    sys.exit(-1)

sys.exit(0)
