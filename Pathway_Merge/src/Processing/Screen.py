'''
Created on Aug 27, 2013

@author: agross
'''
import pandas as pd

import Processing.Helpers as H
import Stats.Survival as Surv
from Stats.Scipy import rev_kruskal, fisher_exact_test

def fc(hit_vec, response_vec):
    f = response_vec.groupby(hit_vec).median()
    return f[0] > f[1]

def mut_filter(df, rate, binary_cutoff=12):
    '''
    Filter out mutation features, ensuring that a feature
    is not entirely an artifact of mutation rate.
    '''
    df = df[df.sum(1) >= binary_cutoff]
    cc = H.screen_feature(rate, rev_kruskal, df)
    
    fc_apply = lambda s: fc(s, rate)
    direction = df.apply(fc_apply, axis=1)
    direction.name = 'direction'
    
    cc = cc.join(direction)
    cc = cc[cc.direction==False]
    
    df = df.ix[H.true_index(cc.p > .01)]
    df = df.dropna(axis=1)
    return df

def cn_filter(df, binary_cutoff=12):
    '''
    Extract copy number features.
    '''
    del_df = (df.ix['Deletion'].dropna(1) < 0).astype(int)
    del_df = del_df[del_df.sum(1) >= binary_cutoff]
    del_df.index = del_df.index.droplevel(1)
    del_df = del_df.T
    amp_df = (df.ix['Amplification'].dropna(1) > 0).astype(int)
    amp_df = amp_df[amp_df.sum(1) >= binary_cutoff]
    amp_df.index = amp_df.index.droplevel(1)
    amp_df = amp_df.T
    return amp_df, del_df

def process_real(df):
    '''
    Process real valued feature into binary feature.
    '''
    df_c = df.copy()
    df_c = df_c.apply(lambda s: H.to_quants(s, std=1), axis=1)
    df_c = df_c > 0
    if type(df.index) == pd.MultiIndex:
        df_c.index = map(lambda s: '_'.join(s), df_c.index)
    return df_c.T

def get_interactions(df, cov_df, surv, test):
    binary = df[df.T.describe().ix['unique'] == 2]
    n_tests = (len(binary) * (len(binary) - 1)) / 2
    s = pd.DataFrame({(a,b): Surv.interaction(v1,v2, surv) 
                          for a,v1 in binary.iterrows()
                          for b,v2 in binary.iterrows()
                          if (a < b)
                          and fisher_exact_test(v1,v2).ix['p'] < (.3 / n_tests)}).T
    int_pairs =  s.ix[s.p < 1].sort('p')
    
    int_associations = {}
    for p,vals in int_pairs.iterrows():
        combo = H.combine(binary.ix[p[0]], binary.ix[p[1]])
        vec = combo == vals['interaction']
        int_associations[p] = test(vec, surv, cov_df) 
    int_associations = pd.DataFrame(int_associations).T
    return s, int_associations

def binary_filter_fx(s):
    vc = s.value_counts()
    if len(vc) != 2:
        return -1
    else:
        return vc.min()

def filter_binary(df, cutoff):
    df = df.dropna(how='all')
    vc = df.apply(binary_filter_fx, axis=1)
    binary = df[vc > cutoff]
    return binary

