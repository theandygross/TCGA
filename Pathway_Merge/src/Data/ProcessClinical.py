'''
Created on Jun 11, 2013

@author: agross
'''
import pandas as pd
import numpy as np

def to_date(s):
    try:
        return pd.datetime(int(s['yearofformcompletion']), int(s['monthofformcompletion']), 
                         int(s['dayofformcompletion']))
    except:
        return np.nan
def fix_date(df):
    df['form_completion'] = df.apply(to_date, 1)
    del df['yearofformcompletion']
    del df['monthofformcompletion']
    del df['dayofformcompletion']
    return df
    
def format_drugs(br):
    drugs = br.select(lambda s: s.startswith('patient.drugs.drug'))
    if len(drugs) == 0:
        return
    drug_annot = pd.MultiIndex.from_tuples(map(lambda s: tuple(s.split('.')[2:4]), drugs.index))
    drugs.index = drug_annot
    drugs = drugs.groupby(level=[0,1]).first() ##get rud of dups, need to fix
    drugs = drugs.stack().unstack(level=1)
    drugs.index = drugs.index.swaplevel(0,1)
    drugs = drugs.sort_index()
    drugs = fix_date(drugs)    
    return drugs

def format_followup1(br):
    followup = br.select(lambda s: s.startswith('patient.followups.followup'))
    if len(followup) == 0:
        return
    idx = pd.MultiIndex.from_tuples(map(lambda s: tuple(s.split('.')[2:5]), followup.index))
    followup.index = idx
    followup = followup.stack().unstack(level=2)
    followup.index = followup.index.reorder_levels([2,0,1])
    followup = followup.sort_index()
    followup = fix_date(followup)
    return followup

def format_followup(br):
    followup = br.select(lambda s: s.startswith('patient.followups.followup'))
    if len(followup) == 0:
        return
    followup['field'] = pd.Series(map(lambda s: s.split('.')[-1], followup.index), 
                                  followup.index)
    followup['label'] = pd.Series(map(lambda s: '.'.join(s.split('.')[:-1]), followup.index),
                                  followup.index)
    followup = followup.set_index(['field','label']).stack().unstack(level='field')
    followup.index = followup.index.reorder_levels([1,0])
    followup = followup.sort_index()
    followup = fix_date(followup)
    return followup

def format_radiation(br):
    followup = br.select(lambda s: s.startswith('patient.radiations.radiation'))
    if len(followup) == 0:
        return
    idx = pd.MultiIndex.from_tuples(map(lambda s: tuple(s.split('.')[2:4]), followup.index))
    followup.index = idx
    followup = followup.stack().unstack(level=1)
    followup.index = followup.index.reorder_levels([1,0])
    followup = followup.sort_index()
    followup = fix_date(followup)
    return followup

def get_meta_data(br, keeper_columns):  
    keeper_columns = map(lambda s: 'patient.' + s, keeper_columns)
    meta_data = br.ix[keeper_columns].rename(index=lambda s: s.split('.')[1])
    meta_data = meta_data.T
    meta_data['age'] = meta_data.daystobirth.astype(float) / -365.
    meta_data['deceased'] = (meta_data.vitalstatus.dropna() == 'deceased').astype(int)
    return meta_data

def format_cqcf(br):
    cqcf = br.select(lambda s: s.startswith('patient.clinicalcqcf'))
    cqcf.index = cqcf.index.map(lambda s: s.split('.')[-1])
    cqcf = cqcf.T.dropna(how='all', axis=1)
    return cqcf

def format_clinical_var(br):
    cl = [s for s in br.index if (s.count('.') == 1) and s.startswith('patient')]
    clinical = br.ix[cl]
    clinical.index = clinical.index.map(lambda s: s.split('.')[1])
    clinical = clinical.T.dropna(axis=1, how='all')
    return clinical

def format_survival(clin, followup):
    clin2 = clin.copy()
    clin2.index = pd.MultiIndex.from_tuples([(i,'surgery',0) for i in clin2.index])
    
    f = followup.append(clin2)
    time_cols = list(f.columns.intersection(['daystodeath', 'daystolastfollowup',
                                        'daystolastknownalive']))
    f = f.sort(columns=time_cols, ascending=False)
    
    f = f.groupby(lambda s: s[0][:12]).last()
    
    age = -1*f.daystobirth.dropna().astype(float) / 365.
    timeline = f[time_cols].dropna(how='all')
    timeline['days'] = timeline.astype(float).max(1)
    deceased = f.vitalstatus.dropna() == 'deceased'
    
    events = f.select(lambda s: 'days' in s, 1)
    timeline = events.combine_first(timeline).astype(float)
    timeline[timeline < -1] = np.nan
    survival = pd.concat([timeline.days, deceased], keys=['days','event'], axis=1)
    survival = survival.dropna().stack().astype(float)
    
    pfs = pd.concat([timeline.days, timeline.daystonewtumoreventafterinitialtreatment], axis=1).min(1)
    event = (pfs < timeline.days) | deceased
    pfs = pd.concat([pfs, event], keys=['days','event'], axis=1)
    pfs = pfs.dropna().stack().astype(float)
    
    survival = pd.concat([survival, pfs], keys=['survival','event_free_survival'], axis=1)
    return survival,timeline

def get_clinical(cancer, data_path, patients=None, filtered_patients=None, **params):
    f = '{}stddata/{}/Clinical/{}.clin.merged.txt'.format(data_path, cancer, cancer)
    tab = pd.read_table(f, index_col=0)
    tab.columns = tab.ix['patient.bcrpatientbarcode'].map(str.upper)
    tab = tab.T.sort_index().drop_duplicates().T
    
    drugs = format_drugs(tab)
    followup = format_followup(tab)
    radiation = format_radiation(tab)
    cqcf = format_cqcf(tab)
    clin = format_clinical_var(tab)
    survival,timeline = format_survival(clin, followup)
    
    return clin, drugs, followup, timeline, survival