'''
Created on Oct 2, 2013

@author: agross
'''
import os as os
import pandas as pd
import numpy as np

from Data.ProcessClinical import try_float

def read_clinical_data(path, cancer):
    cancer = cancer.lower()
    na_vals = ['[Completed]', '[Not Available]', '[Not Applicable]', 'null']
    pat = pd.read_table(path + 'clinical_patient_{}.txt'.format(cancer),
                        index_col=0, na_values=na_vals)
    f = pat.dropna(axis=1, how='all')
    for fu in os.listdir(path):
        if 'clinical_follow_up' not in fu:
            continue
        followup = pd.read_table(path + fu, index_col=0, na_values=na_vals)
        f = pd.concat([f, followup])
    f.columns = f.columns.map(lambda s: s.replace('_', '').lower())
    
    time_vars = ['daystolastfollowup', 'daystolastknownalive',
                 'daystodeath']
    time_cols = list(f.columns.intersection(time_vars))
    
    # f['vitalstatus'] = f['vitalstatus'].map(lambda s: s in 
    #                                        ['DECEASED','Dead','deceased'], 
    #                                        na_action='skip')
    f['vitalstatus'] = (f['daystodeath'].isnull() == True)
    
    f = f.sort(columns=['vitalstatus'] + time_cols, ascending=True)
    f = f.groupby(lambda s: s[:12], axis=0).last()
    return f
    

def format_survival_data(timeline, clin):
    timeline = timeline.applymap(try_float)
    timeline['days'] = timeline.astype(float).max(1)
    
    # deceased = clin.vitalstatus.dropna().isin(['DECEASED','Dead','deceased', False])
    deceased = timeline.daystodeath.isnull() == False
    
    # timeline = events.combine_first(timeline).astype(float)
    timeline[timeline < -1] = np.nan
    survival = pd.concat([timeline.days, deceased], keys=['days', 'event'],
                         axis=1)
    survival = survival.dropna().stack().astype(float)
    
    f_vars = ['days', 'daystonewtumoreventafterinitialtreatment',
              'daystotumorprogression', 'daystotumorrecurrence',
              'daystonewtumoreventafterinitialtreatment',
              'daystonewtumoreventafterinitialtreatment']
    followup_cols = list(timeline.columns.intersection(f_vars))
    pfs = timeline[followup_cols].min(1)
    # pfs = pd.concat([timeline.days, timeline[followup_cols]], axis=1).min(1)
    event = ((pfs < timeline.days) + deceased) > 0
    pfs = pd.concat([pfs, event], keys=['days', 'event'], axis=1)
    pfs = pfs.dropna().stack().astype(float)
    
    survival = pd.concat([survival, pfs],
                         keys=['survival', 'event_free_survival'],
                         axis=1)
    return timeline, survival

def update_clinical_object(clinical, path):
    cancer = clinical.cancer
    f = read_clinical_data(path, cancer)
    fix_cols = lambda df: df.columns.map(lambda s: s.replace('_', '').lower())
    f.columns = fix_cols(f)
    clinical.clinical.columns = fix_cols(clinical.clinical)
    clinical.timeline.columns = fix_cols(clinical.timeline)
    clinical.stage.columns = fix_cols(clinical.stage)
    
    clin = f.combine_first(clinical.clinical).ix[:, clinical.clinical.columns]
    timeline = f.combine_first(clinical.timeline).ix[:, clinical.timeline.columns]
    stage = f.combine_first(clinical.stage).ix[:, clinical.stage.columns]    
    timeline, survival = format_survival_data(timeline, clin)
    timeline['age'] = f.ageatinitialpathologicdiagnosis.astype(float)
    
    clinical.clinical = clin
    clinical.stage = stage
    clinical.survival = survival
    clinical.timeline = timeline
    clinical.artificially_censor(5)
    return clinical

