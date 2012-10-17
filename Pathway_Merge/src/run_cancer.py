#! /cellar/users/agross/semillon/epd/bin/python

import sys
import pickle as pickle

from Data.Pathways import read_in_pathways
from Processing.Clinical import ClinicalObject

PATHWAY_FILE = '/cellar/users/agross/Data/GeneSets/c2.cp.v3.0.symbols_edit.csv'
gene_sets, gene_lookup = read_in_pathways(PATHWAY_FILE)

cancer = sys.argv[1]
data_type = sys.argv[2]
data_path = sys.argv[3]
report_ext = sys.argv[4]
drop_pc = True

try:
	cancer_obj = ClinicalObject(cancer, data_path, gene_sets, data_type,
								drop_pc)
	if False and cancer_obj.data_type in ['mutation','amplification','deletion']:
		cancer_obj.filter_bad_pathways(gene_lookup)
	cancer_obj.report_folder = cancer_obj.data_path + cancer_obj.data_type + '_' + report_ext
	del cancer_obj.tests #Can't pickle functions at module level, regenerate with Clinical.get_tests_x
	pickle.dump(cancer_obj, open(cancer_obj.report_folder + '/ClinicalObject.p', 'wb'))
	print cancer + ' sucess!'
except:
	print cancer + ' fail!'
	sys.exit(-1)

sys.exit(0)
