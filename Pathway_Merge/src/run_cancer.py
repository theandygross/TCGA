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
SURVIVAL_TESTS = {'survival' : {'event_var' : 'deceased', 'time_var' : 'days', 
                                'covariates' : ['age', 'rate']},
                  'event_free_survival' : {'event_var' : 'event', 
                                           'time_var' : 'event_free_survival', 
                                           'covariates' : ['age', 'rate']}
                  }

cancer = sys.argv[1]
data_type = sys.argv[2]
data_path = sys.argv[3]
report_ext = sys.argv[4]
drop_pc = True

try:
	cancer_obj = ClinicalObject(cancer, data_path, gene_sets, data_type,
								survival_tests=SURVIVAL_TESTS, 
								real_variables=REAL_VARIABLES,
                        		binary_variables=BINARY_VARIABLES,
                        		drop_pc=drop_pc)
	if cancer_obj.data_type in ['mutation','amplification','deletion']:
		cancer_obj.filter_bad_pathways(gene_lookup)
	cancer_obj.report_folder = (cancer_obj.data_path + cancer_obj.data_type + 
								'_' + report_ext)
	del cancer_obj.tests #Can't pickle functions at module level, regenerate with Clinical.get_tests_x
	pickle.dump(cancer_obj, open(cancer_obj.report_folder + '/ClinicalObject.p', 
								'wb'))
	print cancer + ' sucess!'
except:
	print cancer + ' fail!'
	sys.exit(-1)

sys.exit(0)
