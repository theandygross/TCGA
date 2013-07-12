'''
Created on Jul 2, 2013

@author: agross
'''

import pickle as pickle
from Data.Containers import Dataset
import Data.Firehose as FH



class CNDataset(Dataset):
    '''
    Inherits from Dataset class.  Adds some added processing for mutation
    data.
    '''
    def __init__(self, run, cancer, cn_type, patients=None):
        '''
        '''
        Dataset.__init__(self, cancer.path, cn_type, compressed=False)
        min_pat = run.parameters['min_patients']
        if cn_type == 'CN_broad':
            self.df = FH.get_gistic(run.data_path, cancer.name, 
                                    min_patients=min_pat)
            if patients is not None:
                self.df = self.df.ix[:, patients].dropna(1, how='all')
            self.features = self.df
            
                
def initialize_cn(cancer_type, report_path, cn_type, patients=None, save=True):
    '''
    Initialize mutation data for down-stream analysis.
    '''
    run = pickle.load(open(report_path + '/RunObject.p', 'rb'))
    cancer = run.load_cancer(cancer_type)

    data = CNDataset(run, cancer, cn_type, patients)
        
    if save is True:
        data.save()
    
    data.uncompress()
    return data