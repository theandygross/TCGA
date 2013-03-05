#! python 

import sys

import pickle as pickle
from Processing.Tests import run_feature_matrix

test = pickle.load(open('test.p','rb'))
test.first_pass = test.get_first_pass()
test.full_test = test.get_full_test()

perm_file = sys.argv[1]
mat = pickle.load(open(perm_file, 'rb'))

vals = run_feature_matrix(mat, test)
pickle.dump(vals, open('res_' + perm_file, 'wb'))

sys.exit(0)
