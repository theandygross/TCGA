'''
Created on Jun 18, 2013

@author: agross
'''
import os as os
import pickle as pickle
import pandas as pd

# from Reports.NotebookTools import *
from Stats.Scipy import *
from Stats.Survival import *
# from Reports.Figures import *
from Processing.Helpers import *
from Figures.Helpers import *
from Figures.Pandas import *
from Figures.Boxplots import *
# from Figures.R_Wrappers import *
from Figures.Survival import draw_survival_curve, survival_and_stats
from Figures.Survival import draw_survival_curves

pd.set_option('precision', 3)
pd.set_option('display.line_width', 100)
pd.set_option('display.width', 300)
plt.rcParams['font.size'] = 12

'''Color schemes for paper taken from http://colorbrewer2.org/'''
colors = plt.rcParams['axes.color_cycle']
colors_st = ['#CA0020', '#F4A582', '#92C5DE', '#0571B0']
colors_th = ['#E66101', '#FDB863', '#B2ABD2', '#5E3C99']

def get_run(firehose_dir, version='Latest'):
    '''
    Helper to get a run from the file-system. 
    '''
    path = '{}/ucsd_analyses'.format(firehose_dir)
    if version is 'Latest':
        version = sorted(os.listdir(path))[-1]
    run = pickle.load(open('{}/{}/RunObject.p'.format(path, version), 'rb'))
    return run
