'''
Created on Oct 25, 2012

@author: agross
'''
from pandas import read_csv, read_table
import Data.Methylation as Methylation

PROBE_PATH = '/cellar/users/agross/Data/MethylationAge/Model_AllData.csv'

def get_age_signal(data_path, cancer):
    '''
    Goes through stddata directory, reads in the picked beta values
    for the reduced set of methylation probes used in Greg's aging
    model, then calculates the predicted age as well as the apearant
    methylomic aging rate (AMAR).
    '''
    clinical = read_csv(data_path + '/'.join(['ucsd_processing', cancer, 
                                      'Clinical','compiled.csv']), index_col=0)
    age = clinical.age
    
    folder = data_path + '/'.join(['ucsd_processing', cancer, 'methylation450',''])
    table = read_table(folder + 'beta_values_picked.txt', index_col=0)
    table = table.rename(columns=lambda s: s[:12] if s!=table.columns[0] 
                             else 'symbol')
    table = table.select(lambda s: s[:4] == 'TCGA', axis=1)
    good = table.index[(table > -1).sum(1) > 0]
    table = table.ix[good]
    
    probes = read_csv(PROBE_PATH, index_col=0)    
    c = probes.ix[good]['Coefficient']
    meth_age = table.apply(lambda s: s.dot(c))
    amar = meth_age.div(age)
    return meth_age, amar

def run_all_cancers(firehose_path, date, recalc=False):
    '''
    Pulls probes from big methylation files and puts them into
    beta_values_picked.txt files next to original data.
    '''
    probes = read_table(PROBE_PATH, index_col=0)
    Methylation.run_all_cancers(firehose_path, date, probeset=list(probes.index),
                                recalc=recalc)

