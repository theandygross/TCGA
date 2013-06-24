'''
Created on Jun 18, 2013

@author: agross
'''
import os as os
import pickle as pickle
import pandas as pd

from Reports.NotebookTools import *
from Processing.Tests import *
from Reports.Figures import *
from Processing.Helpers import *
from Figures.Survival import draw_survival_curve

pd.set_option('precision',3)
pd.set_option('display.line_width', 100)
pd.set_option('display.width', 300)

def get_run(date):
    result_path = '/cellar/data/TCGA/Firehose__{}/ucsd_analyses'.format(date)
    run = sorted(os.listdir(result_path))[-1]
    run = pickle.load(open('/'.join([result_path, run, 'RunObject.p']), 'rb'))
    return run