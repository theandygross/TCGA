'''
Created on Jun 18, 2013

@author: agross
'''
import os as os
import pickle as pickle
import pandas as pd

from Reports.NotebookTools import *
from Stats.Scipy import *
from Stats.Survival import *
#from Reports.Figures import *
from Processing.Helpers import *
from Figures.Pandas import *
from Figures.Boxplots import *
from Figures.R_Wrappers import *
from Figures.Survival import draw_survival_curve

pd.set_option('precision',3)
pd.set_option('display.line_width', 100)
pd.set_option('display.width', 300)

def get_run(firehose_dir, version='Latest'):
    '''
    Helper to get a run from the file-system. 
    '''
    path = '{}/ucsd_analyses'.format(firehose_dir)
    if version is 'Latest':
        version = sorted(os.listdir(path))[-1]
    run = pickle.load(open('{}/{}/RunObject.p'.format(path, version), 'rb'))
    return run