"""
Created on Jun 11, 2013

@author: agross
"""
import pandas as pd
import numpy as np


def to_date(s):
    """
    Pulls year, month, and day columns from clinical files and 
    formats into proper date-time field.
    """
    try:
        return pd.datetime(int(s['yearofformcompletion']),
                           int(s['monthofformcompletion']),
                           int(s['dayofformcompletion']))
    except:
        return np.nan


def fix_date(df):
    """
    Translate date to date-time, get rid of old columns.
    """
    try:
        df['form_completion'] = df.apply(to_date, 1)
        del df['yearofformcompletion']
        del df['monthofformcompletion']
        del df['dayofformcompletion']
    except:
        pass
    return df


def try_float(s):
    try:
        return float(s)
    except:
        return np.nan


def format_drugs(br):
    """
    Format drug rows in merged clinical file from Firehose.
    The data consists of one or more drug entries for each patient.
    Here we use a MultiIndex with Patient Barcode on level 0 and 
    the drug administration on level 1. 
    
    Input
        br: clinical DataFrame with patient bar-codes on the columns
    """
    drugs = br.select(lambda s: s.startswith('patient.drugs.drug'))
    if len(drugs) == 0:
        return
    ft = pd.MultiIndex.from_tuples  # long Pandas names
    drug_annot = ft(map(lambda s: tuple(s.split('.')[2:4]), drugs.index))
    drugs.index = drug_annot
    # drugs = drugs.groupby(level=[0,1]).first() ##get rud of dups, need to fix
    l = drugs.sort_index().groupby(level=[0, 1], axis=0)
    drugs = pd.DataFrame({i: s.ix[s.count(1).idxmax()] for i, s in l}).T
    drugs.index = ft(drugs.index)
    drugs = drugs.stack().unstack(level=1)
    drugs.index = drugs.index.swaplevel(0, 1)
    drugs = drugs.sort_index()
    drugs = fix_date(drugs)    
    return drugs


def format_followup(br):
    """
    Format follow-up rows in merged clinical file from Firehose.
    The data consists of one or more followup entries for each patient.
    Here we use a MultiIndex with Patient Barcode on level 0 and 
    the follow-up number on level 1. 
    
    Input
        br: clinical DataFrame with patient bar-codes on the columns
    """
    row_filter = lambda s: s.startswith('patient.followups.followup')
    followup = br.select(row_filter)
    if len(followup) == 0:
        return
    ft = pd.MultiIndex.from_tuples  # long Pandas names
    followup.index = ft([(s.split('.')[-1], '.'.join(s.split('.')[:-1])) 
                         for s in followup.index])
    followup = followup.stack().unstack(level=0)
    followup.index = followup.index.reorder_levels([1, 0])
    followup = followup.sort_index()
    followup = fix_date(followup)
    return followup


def format_stage(br):
    row_filter = lambda s: s.startswith('patient.stageevent')
    stage = br.select(row_filter)
    stage = stage.dropna(how='all')
    stage = stage.rename(index=lambda s: s.split('.')[-1])
    stage = stage.T
    return stage


def format_radiation(br):
    """
    Format radiation rows in merged clinical file from Firehose.
    The data consists of one or more entries for each patient.
    Here we use a MultiIndex with Patient Barcode on level 0 and 
    the treatment number on level 1. 
    
    Input
        br: clinical DataFrame with patient bar-codes on the columns
    """
    row_filter = lambda s: s.startswith('patient.radiations.radiation')
    followup = br.select(row_filter)
    if len(followup) == 0:
        return
    ft = pd.MultiIndex.from_tuples
    idx = ft([tuple(s.split('.')[2:4]) for s in followup.index])
    followup.index = idx
    followup = followup.stack().unstack(level=1)
    followup.index = followup.index.reorder_levels([1, 0])
    followup = followup.sort_index()
    followup = fix_date(followup)
    return followup


def format_clinical_var(br):
    """
    Format clinical variables that are not associated with drug, follow-up,
    or radiation.

    Input
        br: clinical DataFrame with patient bar-codes on the columns
    """
    cl = [s for s in br.index if (s.count('.') == 1)
                              and s.startswith('patient')]
    clinical = br.ix[cl]
    clinical.index = clinical.index.map(lambda s: s.split('.')[1])

    cl = [s for s in br.index if (s.count('.') == 2)
                              and s.startswith('patient.primarypathology')]
    clinical2 = br.ix[cl]
    clinical2.index = clinical2.index.map(lambda s: s.split('.')[2])

    clinical = clinical.append(clinical2)
    clinical = clinical.T.dropna(axis=1, how='all')


    clinical['age'] = clinical.ageatinitialpathologicdiagnosis.astype(float)
    del clinical['ageatinitialpathologicdiagnosis']
    return clinical

def format_survival(clin, followup):
    """
    Format survival for downstream analysis.
    For survival analysis we need to track the time to death/censoring
    as well as the censoring status (censored or deceased) for each patient.
    We use a MultiIndex with Patient Barcode on level 0 and ['days','event']
    on level 1, where days in the time variable and 'event' is the death 
    indicator. Here we extract the standard survival as well as event 
    free survival from the clinical information as well as the patient 
    followup. 
    
    Input
        br: clinical DataFrame with patient bar-codes on the columns
        
    Returns:
        survival: DataFrame consisting of event_free_survival, and survival
                  Series
        timeline: DataFrame of clinical variables related to patient cancer
                  timelines
    """
    clin2 = clin.copy()
    clin2.index = pd.MultiIndex.from_tuples([(i, 'surgery', 0) for i in clin2.index])
    if type(followup) == pd.DataFrame:
        f = followup.append(clin2)
    else:
        f = clin2
    time_vars = ['daystodeath', 'daystolastfollowup', 'daystolastknownalive',
                 'daystonewtumoreventafterinitialtreatment', 'daystotumorprogression',
                 'daystotumorrecurrence']
    time_cols = list(f.columns.intersection(time_vars))
    timeline = f[time_cols].dropna(how='all').astype(float)
    timeline['days'] = timeline.max(1)
    timeline = timeline.groupby(level=0).max()
    deceased = timeline.daystodeath.isnull() == False

    #days = timeline.days[((timeline.days > 7) | (deceased == False))]
    #days = days[days > 0]

    days = timeline.days[timeline.days >= 0]
    survival = pd.concat([days, deceased], keys=['days', 'event'], axis=1)
    survival = survival.dropna().stack().astype(float)

    pfs_var = 'daystonewtumoreventafterinitialtreatment'
    if (followup is not None) and (pfs_var in followup):
        new_tumor = followup[pfs_var].dropna().groupby(level=0).min()
        time_to_progression = pd.concat([new_tumor, timeline.days], 1).min(1)
        time_to_progression = time_to_progression[time_to_progression > 7]
        progression = (deceased | pd.Series(1, index=new_tumor.index))
        pfs = pd.concat([time_to_progression, progression], keys=['days', 'event'],
                        axis=1)
        pfs = pfs.dropna().stack().astype(float)
    else:
        pfs = survival
    survival = pd.concat([survival, pfs], keys=['survival', 'event_free_survival'],
                         axis=1)
    return survival, timeline


def get_clinical(cancer, data_path, patients=None, **params):
    """
    Reads in and formats clinical data for a given tumor type.  
    
    Returns
        clin:     clinical variables
        drugs:    drugs administered, Dataframe with with 
                  (patient, treatment_id) on index
        followup: patient followups, Dataframe with with 
                  (patient, follwup_id) on index
        timeline: patient cancer timeline variables
        survival: DataFrame consisting of event_free_survival, and 
                  survival with (patient, ['days','event']) on index.
        
    """
    f = '{}stddata/{}/Clinical/{}.clin.merged.txt'.format(data_path, cancer,
                                                          cancer)
    tab = pd.read_table(f, index_col=0, low_memory=False)
    tab.columns = tab.ix['patient.bcrpatientbarcode'].map(str.upper)
    tab = tab.T.sort_index().drop_duplicates().T
    
    drugs = format_drugs(tab)
    followup = format_followup(tab)
    stage = format_stage(tab)
    # radiation = format_radiation(tab)
    clin = format_clinical_var(tab)
    survival, timeline = format_survival(clin, followup)
    
    return clin, drugs, followup, stage, timeline, survival
