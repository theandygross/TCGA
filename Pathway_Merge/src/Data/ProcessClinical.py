'''
Created on Jan 23, 2013

@author: agross
'''
import os
from pandas import Series, DataFrame
from pandas import read_table, datetime

from numpy import nan, array


def to_date(s):
    try:
        return datetime(int(s['yearofformcompletion']), int(s['monthofformcompletion']), 
                         int(s['dayofformcompletion']))
    except:
        return nan
    
def format_drugs(br):
    drugs = br.select(lambda s: s.startswith('patient.drugs.drug'))
    if len(drugs) == 0:
        return
    drug_annot = array(map(lambda s: (s.split('.')[2:4]), drugs.index))
    drug_annot = DataFrame(drug_annot, index=drugs.index, columns=['drug', 'info'])
    drug_annot['drug'] = drug_annot['drug'].map(lambda s: int(s.split('-')[-1]) if '-' 
                                                in s else 1)
    #barcodes = Series(br.ix['patient.bcrpatientbarcode'], name='patient')
    #drugs = drugs.rename(columns=barcodes).join(drug_annot)
    drugs = drugs.join(drug_annot)
    
    drugs = DataFrame({p: drugs.ix[:,[p,'drug','info']].dropna().groupby(['drug','info']).first().ix[:,0]
                       for p in drugs.columns[:-2]}).T
    drugs.columns = drugs.columns.reorder_levels([1,0])
    drugs = drugs.stack()
    
    drugs['form_completion'] = drugs.apply(to_date, 1)
    del drugs['yearofformcompletion']
    del drugs['monthofformcompletion']
    del drugs['dayofformcompletion']
    
    return drugs

def format_followup(br):
    followup = br.select(lambda s: s.startswith('patient.followups.followup'))
    if len(followup) == 0:
        return
    follow_annot = array(map(lambda s: ('.'.join(s.split('.')[2:4]), 
                                        s.split('.')[4]), followup.index))
    follow_annot = DataFrame(follow_annot, index=followup.index, 
                             columns=['followup', 'info'])
    #barcodes = Series(br.ix['patient.bcrpatientbarcode'].str.upper(), name='patient')
    #followup = followup.rename(columns=barcodes).join(follow_annot)
    followup = followup.join(follow_annot)
    followup = DataFrame({p: followup.ix[:,[p,'followup','info']].dropna().groupby(['followup','info']).first().ix[:,0]
                       for p in followup.columns[:-2]}).T
    followup.columns = followup.columns.reorder_levels([1,0])
    followup = followup.stack()
    
    followup['form_completion'] = followup.apply(to_date, 1)
    del followup['yearofformcompletion']
    del followup['monthofformcompletion']
    del followup['dayofformcompletion']
    return followup

def get_meta_data(br, keeper_columns):  
    keeper_columns = map(lambda s: 'patient.' + s, keeper_columns)
    #meta_data = br.ix[keeper_columns].rename(columns=br.ix['patient.bcrpatientbarcode'], 
    #                                         index=lambda s: s.split('.')[1])
    meta_data = br.ix[keeper_columns].rename(index=lambda s: s.split('.')[1])
    meta_data = meta_data.T
    meta_data['age'] = meta_data.daystobirth.astype(float) / -365.
    meta_data['deceased'] = (meta_data.vitalstatus.dropna() == 'deceased').astype(int)
    return meta_data

def get_timeline_o(followup, meta_data, event_cols):
    timeline = meta_data[['daystodeath', 'daystolastfollowup']].astype(float)
    days = Series(timeline[['daystolastfollowup','daystodeath']].min(1), name='days')
    
    def format_event_free_survival(timeline):
        df = followup.groupby(level=0).last()
        time_cols = [c for c in df.columns if c in event_cols]   
        if time_cols == []:
            return timeline     
        t = df[time_cols].astype(float)
        time_to_event = t[t > 0].min(1).dropna()
        time_to_event.name = 'days_to_event'
        new_tumor = Series(time_to_event > 0, name='new_tumor')
        event = (new_tumor + meta_data['deceased']).clip_upper(1).fillna(0)
        event  = Series(event, name='event').astype(float)
        timeline = timeline.join(df[time_cols]).astype(float)
        event_free_survival = Series(timeline.min(1), name='event_free_survival')
        timeline = timeline.join(event).join(meta_data['deceased'])
        timeline = timeline.join(event_free_survival)
        return timeline
        
    if followup is not None:
        timeline = format_event_free_survival(timeline)        
    timeline = timeline.join(meta_data['daystobirth']).join(days)
    return timeline

def get_timeline(followup, meta_data, event_cols):
    try:
        df = followup.groupby(level=0).last()
        time_cols = [c for c in df.columns if c in event_cols] 
        timeline = meta_data.join(df[time_cols], lsuffix='_f')
    except AttributeError:
        timeline = meta_data
    
    timeline = timeline.select(lambda s: 'days' in s and '_f' not in s, axis=1)       
    timeline = timeline.astype(float)
    timeline = timeline.ix[:,(timeline > 0).sum() > 0].dropna(how='all', axis=1)
    
    s_vars = ['daystodeath', 'daystolastfollowup']
    last_known = Series(timeline[s_vars].max(1), name='ln')
    pfs = timeline.drop(s_vars, 1).join(last_known).min(1)
    pfs = pfs[pfs > 30] #30 day buffer
    timeline['event_free_survival'] = pfs
    
    meta_data['days'] = meta_data[s_vars].astype(float).min(1)
    survival = meta_data[['days','deceased']].rename(columns={'deceased': 'event'})
    survival = survival.dropna().stack().astype(float)
    
    event = Series((timeline.event_free_survival != timeline.daystolastfollowup) +
                   ((timeline.event_free_survival == timeline.daystodeath)), name='event')
    s = timeline[['event_free_survival']].join(event).rename(columns={'event_free_survival': 
                                                                      'days'})
    evs = s.dropna().stack().astype(float)
    
    survival = DataFrame.from_dict({'survival':survival, 'event_free_survival': evs})

    return timeline, survival

def her2_combo(s):
    if s['her2_fish'] in ['positive','negative']:
        return s['her2_fish']
    elif s['her2_chem'] in ['positive','negative']:
        return s['her2_chem']
    elif s['her2_fish'] not in [nan,'performed but not available','not performed']:
        return s['her2_fish']
    else:
        return s['her2_chem']
    
def process_specific_features(c):
    '''
    Yes this is sort of hard coded but I don't want to pass functions around
    globally as I can't pickle them. 
    '''
    if 'tumor_stage' in c:
        c['tumor_stage'] = c['tumor_stage'].map(lambda s: s.replace('a','').replace('b','').replace('c',''),
                                                na_action='ignore')
    if 'tumor_grade' in c:
        c['tumor_t1t2']= c['tumor_grade'].map(lambda s: ('t1' in s) or ('t2' in s), na_action='ignore')
    if 'lymphnode_stage' in c:
        c['lymphnode_n0n1']= c.lymphnode_stage.map(lambda s: ('n0' in s) or ('n1' in s), 
                                                   na_action='ignore')
    if 'smoker' in c:
        c['smoker'] = c.smoking_status != 'lifelong non-smoker'
    if 'aml_cytogenics' in c:
        c['normal_cyto'] = c.aml_cytogenics == 'normal'
    if 'neo_depth' in c:
        c['neo_area'] = c[['neo_depth', 'neo_width', 'neo_length']].astype(float).product(1)
    if 'her2_fish' in c:
        c['her2'] = c.apply(her2_combo, 1)
        c['triple_neg'] = (c[['ER','PR','her2']] == 'negative').sum(1) == 3
    return c

def get_clinical(cancer, data_path, patients, filtered_patients, **params):
    p = params
    path = (data_path + '/'.join(['stddata', cancer,'Clinical', cancer]) + 
            '.clin.merged.txt')
    br = read_table(path, index_col=0)
    br = br.rename(columns=br.ix['patient.bcrpatientbarcode'].str.upper())
    
    if patients is not None:
        br = br[[pat for pat in br.columns if pat in patients]]
    elif filtered_patients is not None:
        br = br[[pat for pat in br.columns if pat not in filtered_patients]]  
          
    drugs = format_drugs(br)
    followup = format_followup(br)
    
    meta_data = get_meta_data(br, p['meta_cols'])
    timeline, survival = get_timeline(followup, meta_data, p['event_cols'])
    #generic_clinical = timeline.join(meta_data[p['meta_keep']])
    generic_clinical = meta_data[p['meta_keep']]
    
    br = br.rename(index=lambda s: '.'.join(s.split('.')[1:]))
    c = br.ix[p['feature_dict']].dropna(how='all')
    c = c.rename(index=p['feature_dict']).T
    c = process_specific_features(c)
    cl = c.join(generic_clinical)
    
    if followup is not None:
        df = followup.groupby(level=0).last().T
        followup_vars = df.ix[p['followup_dict']].dropna(how='all')
        followup_vars = followup_vars.rename(index=p['followup_dict']).T
        cl = cl.join(followup_vars)
        
    return cl, drugs, followup, timeline, survival