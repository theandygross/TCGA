'''
Created on Oct 25, 2012

@author: agross
'''
from pandas import read_csv, read_table

PROBE_PATH = '/cellar/users/menzies/Work/MethylCs/suppTable3_fromGH_20120517.txt'

def get_age_signal(stddata_path, clinical):
    folder = stddata_path + 'methylation__humanmethylation450__jhu_usc_edu__Level_3/'
    table = read_csv(folder + 'beta_values_picked.txt', index_col=0)
    table = table.select(lambda s: s[:4] == 'TCGA', axis=1)
    good = table.index[(table > -1).sum(1) > 0]
    table = table.ix[good]
    
    probes = read_table(PROBE_PATH, index_col=0)    
    c = probes.ix[good]['Coefficient']
    t = table.apply(lambda s: s.dot(c))
    amar = t.div(clinical.age)
    return t,amar