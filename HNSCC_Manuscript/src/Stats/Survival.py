'''
Created on Jun 30, 2013

@author: agross
'''
import pandas as pd
import numpy as np
import scipy.stats as stats

import rpy2.robjects as robjects
from pandas.rpy.common import convert_to_r_dataframe, convert_robj
from Processing.Helpers import get_vec_type, bhCorrection, powerset
from Processing.Helpers import match_series, combine
from Stats.Scipy import fisher_exact_test

survival = robjects.packages.importr('survival')
base = robjects.packages.importr('base')
mass = robjects.packages.importr('MASS')

robjects.r.options(warn=-1);
zz = robjects.r.file("all.Rout", open="wt")
robjects.r.sink(zz, type='message')

def log_rank(feature, surv):
    '''
    Perform log-rank test using r.survival.survdiff function.
    '''
    feature = sanitize_lr(feature)
    if type(feature) is type(None):
        return pd.Series(index=['chi2', 'p'])
    fmla = robjects.Formula('Surv(days, event) ~ feature')
    #use cox function to extract model
    m = get_cox_ph(surv, feature, formula=fmla)
    r_data = m.rx2('call')[2]
    s = survival.survdiff(fmla, r_data)
    p = stats.chi2.sf(s.rx2('chisq')[0], len(feature.unique()) - 1)
    return pd.Series({'chi2': s.rx2('chisq')[0], 'p': p})

def log_rank_more(feature, surv):
    '''
    Perform log-rank test using r.survival.survdiff function.
    '''
    feature = sanitize_lr(feature)
    if type(feature) is type(None):
        return pd.Series(index=['chi2', 'p'])
    fmla = robjects.Formula('Surv(days, event) ~ feature')
    #use cox function to extract model
    m = get_cox_ph(surv, feature, formula=fmla)
    r_data = m.rx2('call')[2]
    s = survival.coxph(fmla, r_data)
    
    b = base.summary(s)
    print b

    hazard = convert_robj(b.rx2('conf.int')).ix[0]
    stat = pd.Series(b.rx2('logtest'), index=['stat', 'df', 'p'])
    concordance = pd.Series(b.rx2('concordance'), index=['stat','se'])
    ret = pd.concat([hazard, stat, concordance], keys=['hazard','LR','concordance'])
    return ret


def test_model(p):
    interactions = [t for t in p if ':' in t] 
    for t in interactions:
        term = t.split(':')[-1]
        if term not in p:
            return False
        if 'feature' not in p:
            return False
    return True

def get_models(factors, interactions='just_feature'):
    if interactions == 'just_feature':
        cov = [c for c in factors if c != 'feature']
        models = ['feature:' + c for c in cov]
        models = [p for p in powerset(['feature'] + cov + models) if 
                  test_model(p)]
    elif interactions == True:
        int_terms = [':'.join([a,b]) for a in factors for b in factors 
                     if a<b]
        models = list(powerset(factors + int_terms))
    else:
        models = list(powerset(factors))
    models = map(lambda s: ' + '.join(s), models)
    models[0] = '1'
    return models
    
def cox_model_selection(surv, feature=None, covariates=None, interactions=True):
    df, factors = process_covariates(surv, feature, covariates)
    models = get_models(factors, interactions)      
    ll = {}
    for m in models:
        fmla = robjects.Formula('Surv(days, event) ~ ' + m)
        s = survival.coxph(fmla, df)
        ll[m] = max(s.rx2('loglik'))
    ll = pd.Series(ll)
    
    dof = pd.Series(ll.index, ll.index)
    dof = dof.apply(lambda s: s.count('+') + 1)
    q = 3
    AIC = (-2*ll) + (dof*q)
    best_model = AIC.idxmin()
    best_model = robjects.Formula('Surv(days, event) ~ ' + best_model)
    best_res = survival.coxph(best_model, df)
    return best_res

def LR_test(full, reduced):
    '''
    Perform Likelihood ratio test on two R models. 
    '''
    full_ll = list(full.rx2('loglik'))
    reduced_ll = list(reduced.rx2('loglik'))
    assert full_ll[0] == reduced_ll[0]
    if len(full_ll) == 1:
        return 1.
    full_df = len(full.rx2('coefficients')) 
    if len(reduced_ll) == 1:
        return stats.chi2.sf(2*full_ll[1] - 2*full_ll[0], full_df)
    reduced_df = len(reduced.rx2('coefficients'))
    df = max(full_df - reduced_df, 1)
    return stats.chi2.sf(2*full_ll[1] - 2*reduced_ll[1], df)

def sanitize_lr(feature):
    if feature is None:
        return feature
    if len(feature.unique()) <= 1:
        return feature.map(lambda s: np.nan)
    if feature.dtype not in ['str','object','bool']:
        return feature
    try:
        feature = feature.astype(float)
        return feature
    except:
        pass
    if len(feature.value_counts()) > 5:
        try:
            feature = feature.astype(float)
            return feature
        except:
            pass
    return feature
    

def get_formula(factors, get_interactions=True):
    if len(factors) > 1:
        interactions = ' + '.join(factors)
        if get_interactions:
            interactions += ' + '
        if get_interactions == 'just_first':
            interactions += ' + '.join((':'.join([factors[0],b]) for b in factors[1:]))
        else: #all pairs
            interactions += ' + '.join((':'.join([a,b]) for a in factors 
                                      for b in factors if a<b))
    elif len(factors) == 1: 
        interactions = factors[0]
    else:
        interactions = '1'
    fmla = 'Surv(days, event) ~ {}'.format(interactions)
    return fmla

def process_factors(clinical, hit_vec=None, covariates=[]):
    if not all([cov in clinical for cov in covariates]):
        covariates = [cov for cov in covariates if cov in clinical]
    if type(hit_vec) != type(None):
        factors = ['feature'] + covariates
        df = clinical.join(pd.Series(hit_vec, name='feature'))
    else:
        factors = covariates
        df = clinical
    return df, factors

def process_covariates(surv, feature=None, cov=None):
    '''
    Coerce covariates and feature into format suitable for R's
    survival functions. 
    '''
    if type(feature) is type(None):
        feature = pd.Series(index=surv.index.levels[0])
    if type(cov) is type(None):
        cov = pd.DataFrame(index=feature.index)
    if type(cov) == pd.Series:
        cov = pd.concat([cov], axis=1)
    elif type(cov) == list:
        assert map(type, cov) == ([pd.Series] * len(cov))
        cov = pd.concat(cov, axis=1)
    cov = cov.apply(sanitize_lr)
    feature = sanitize_lr(feature)
    c_real = cov.ix[:,cov.dtypes.isin([np.dtype(float), np.dtype(int)])]
    c_real = (c_real - c_real.mean()) / c_real.std()
    if c_real.shape[1] > 0:
        cov[c_real.columns] = c_real
    cov = cov.dropna(1, how='all')
    df = cov.join(surv.unstack()).dropna()
    df['days'] = df['days'] / 365.
    df = df.groupby(level=0).first()
    if len(feature.dropna()) == 0:
        feature = None 
    df, factors = process_factors(df, feature, list(cov.columns))
    df = df[factors + ['days','event']]
    df = df.dropna(axis=1, how='all')
    df = convert_to_r_dataframe(df)
    return df, factors

def get_cox_ph(surv, feature=None, covariates=None, formula=None, 
               interactions=True, get_model=True, print_desc=False):
    '''
    Fit a cox proportial hazzards model to the data.
    Returns a p-value on the hit_vec coefficient. 
    ---------------------------------------------------
    clinical: DataFrame of clinical variables
    hit_vec: vector of labels to test against
    covariates: names of covariates in the cox model,
                (must be columns in clinical DataFrame)
    '''
    if formula is None:
        s = cox_model_selection(surv, feature, covariates, interactions)
    else:
        df, _ = process_covariates(surv, feature, covariates)
        fmla = robjects.Formula(formula)
        s = survival.coxph(fmla, df)
    
    if print_desc:
        print '\n\n'.join(str(s).split('\n\n')[-2:])
        
    if get_model:
        return s
    
def get_cox_ph_ms(surv, feature=None, covariates=None, return_val='LR', 
                  null_model=None, formula=None, get_model=True,
                  interactions=True):
    '''
    Fit a cox proportial hazzards model to the data.
    Returns a p-value on the hit_vec coefficient. 
    ---------------------------------------------------
    clinical: DataFrame of clinical variables
    hit_vec: vector of labels to test against
    covariates: names of covariates in the cox model,
                (must be columns in clinical DataFrame)
    '''
    print_desc = return_val == 'model_desc'
    if covariates is None:
        covariates = pd.DataFrame(index=feature.index)
    elif type(covariates) == list:
        assert map(type, covariates) == ([pd.Series] * len(covariates))
        covariates = pd.concat(covariates, axis=1)
        
    s = get_cox_ph(surv, feature, covariates, formula, interactions,
                   get_model, print_desc)
    if s is None:
        return
    
    results = convert_robj(base.summary(s).rx2('coefficients'))
    
    def set_null_model(feature, covariates):
        null_int = False if interactions == 'just_feature' else interactions
        patients = covariates.index.intersection(feature.dropna().index)
        covariates = covariates.ix[patients]
        null_model = get_cox_ph(surv, covariates=covariates, 
                                interactions=null_int)
        get_cox_ph_ms.null_model = null_model
        get_cox_ph_ms.params = surv.name, str(covariates)
      
    # check if we need to recompute null model
    has_null = hasattr(get_cox_ph_ms, 'null_model') or (null_model)
    if has_null != True:
        set_null_model(feature, covariates)
        
    recalc = get_cox_ph_ms.params !=  surv.name, str(covariates)
    if recalc:
        set_null_model(feature, covariates)
    
    null_model = get_cox_ph_ms.null_model
    LR_p = LR_test(s, null_model)
    
    if type(results) == pd.DataFrame and 'feature' in results.index:
        coef_p = results.ix['feature','Pr(>|z|)']
        hazzard = results.ix['feature','exp(coef)']
    else:
        coef_p, hazzard = np.nan, np.nan
            
    if return_val == 'LR_p':
        return LR_p
    
    elif return_val == 'LR':
        f = str(s.rx2('formula'))
        results = pd.Series({'LR': LR_p, 
                             'feature_p': coef_p, 
                             'hazzard': hazzard, 
                             'fmla': f})
    return results

def get_surv_fit(surv, feature=None, covariates=None, interactions=None,
                 formula=None, time_cutoff=5):
    df, factors = process_covariates(surv, feature, covariates)
    if formula is None:
        fmla = get_formula(factors, interactions)
        fmla = robjects.Formula(fmla)
    else:
        fmla = robjects.Formula(formula)
    
    s = survival.survfit(fmla, df)
    summary = base.summary(s, times=robjects.r.c(time_cutoff))
    res =  convert_robj(summary.rx2('table'))
    
    if type(res) == list:
        r = summary.rx2('table')
        r = pd.Series(r, r.names)
        res = pd.DataFrame({'feature=all': r}).T
        
    res = res.rename(index=lambda idx: idx.split('=')[1])
    res = res[['records','events','median','0.95LCL','0.95UCL']]
    res.columns = pd.MultiIndex.from_tuples([('Stats','# Patients'), 
                                             ('Stats','# Events'), 
                                             ('Median Survival', 'Median'), 
                                             ('Median Survival', 'Lower'),
                                             ('Median Survival', 'Upper')])
    if feature is None:
        for f in ['surv','lower','upper']:
            res[(str(time_cutoff) + 'y Survival', 
                 f.capitalize())] = summary.rx2(f)
    else:
        idx = map(lambda s: s.replace('feature=',''), 
                  summary.rx2('strata').iter_labels())
        
        df = pd.DataFrame({d: list(summary.rx2(d)) for d in 
                           ['strata','surv','lower','upper']},
                          index=idx)
        for f in ['surv','lower','upper']:
            res[(str(time_cutoff) + 'y Survival',
                 f.capitalize())] = df[f]
            
    try:
        res.index = map(int, res.index)
    except:
        pass
    return res

def get_surv_fit_lr(surv, feature=None):
    t = get_surv_fit(surv, feature)
    s = log_rank(feature, surv)
    num_f = len(feature.dropna().unique())
    t[('Log-Rank','chi2')] = [''] * num_f
    t[('Log-Rank','p')] = [''] * num_f
    t = t.append(pd.Series([''] * (8) + [s['chi2'], s['p']], index=t.columns, name=''))
    t = t.sort([('Stats','# Patients')], ascending=False)
    return t
    
'''    
def get_cox_ph(clinical, hit_vec=None, covariates=[], time_var='days',
               event_var='censored'):
    p = get_cox_ph_ms(clinical, hit_vec, covariates, time_var, event_var, 
                      'LR_p')
    return p
'''
    
    
class SurvivalTest(object):
    def __init__(self, surv, covariates):
        self.statistic = ('Full', 'LR')
        self.surv = surv
        self.covariates = covariates
        self.first_pass = self.get_first_pass()
        self.full_test = self.get_full_test()
        
    def get_first_pass(self):
        '''
        Fist pass test for survival tests is basic test with no covariates.
        '''
        def test(feature):
            return get_cox_ph_ms(self.surv, feature, return_val='p_haz', 
                                 formula='Surv(days, event) ~ feature')
        return test
    
    def get_full_test(self):
        '''
        Run Cox-PH with full model.
        '''
        test = lambda feature: get_cox_ph_ms(self.surv, feature, 
                                             covariates=self.covariates, 
                                             return_val='LR')
        return test
    
    def check_feature(self, vec):
        vec_type = get_vec_type(vec)
        if vec_type == 'boolean':
            return vec.value_counts()[1] > 10
        elif vec_type == 'real':
            return vec.count() > 50
        else:
            return False
        
def run_feature_matrix(df, test, fp_cutoff=.5):
    df = df.ix[df.apply(test.check_feature, 1)]
    if hasattr(test, 'first_pass'):
        fp = df.apply(test.first_pass, 1)
        df = df[fp.p < fp_cutoff]
    full = df.apply(test.full_test, 1)
    res = pd.concat([full[['LR', 'fmla']], fp], keys=['Full', 'Univariate'], axis=1)
    if type(res.index[0]) == tuple: #pandas bug
                res.index = pd.MultiIndex.from_tuples(res.index, names=df.index.names) 
    res = res.join(pd.Series(bhCorrection(res[('Full', 'LR')], n=len(fp)), 
                          name=('Full','LR_q')))
    res = res.join(pd.Series(bhCorrection(res[('Univariate', 'p')], n=len(fp)), 
                          name=('Univariate','q')))
    return res.sort_index(axis=1).sort(columns=[('Full', 'LR')])

def stratified_cox(feature, surv, strata):
    fmla = 'Surv(days, event) ~ feature + strata({})'.format(strata.name)
    model = get_cox_ph(surv, feature, covariates=strata, formula=fmla)
    lr = model[3][0]
    p = stats.chi2.sf(lr, 1)
    return pd.Series({'LR': lr, 'p': p})

def cox(feature, surv):
    '''
    Perform log-rank test using r.survival.survdiff function.
    '''
    if feature.dtype in ['str','object','bool']:
        feature = sanitize_lr(feature)
    if type(feature) is type(None):
        return cox(feature)
    fmla = robjects.Formula('Surv(days, event) ~ feature')
    #use cox function to extract model
    s = get_cox_ph(surv, feature, formula=fmla)
    b = base.summary(s)

    hazard = convert_robj(b.rx2('conf.int')).ix[0]
    stat = pd.Series(b.rx2('logtest'), index=['stat', 'df', 'p'])
    concordance = pd.Series(b.rx2('concordance'), index=['stat','se'])
    ret = pd.concat([hazard, stat, concordance], keys=['hazard','LR','concordance'])
    return ret

def get_stats(s):
    b = base.summary(s)
    hazard = convert_robj(b.rx2('conf.int')).ix['feature']
    stat = pd.Series(b.rx2('logtest'), index=['stat', 'df', 'p'])
    concordance = pd.Series(b.rx2('concordance'), index=['stat','se'])
    ret = pd.concat([hazard, stat, concordance], keys=['hazard','LR','concordance'])
    return ret

def cox_screen(df, surv, axis=1):
    if axis==0:
        df = df.T
    c = df.apply(pd.value_counts, axis=1).count(1)
    df = df.ix[c[c>1].index]
    rr = df.apply(lambda s: cox(s.dropna(), surv), axis=1)
    
    rr[('LR', 'q')] = bhCorrection(rr['LR']['p'])
    rr = rr.sort([('LR','q')])
    rr = rr.sortlevel(0, axis=1)
    return rr

def lr_screen(df, surv):
    rr = df.astype(float).apply(log_rank, args=(surv,), axis=1)
    rr['q'] = bhCorrection(rr.p)
    rr = rr.sort('p')
    return rr

def _interaction(a,b, surv):
    a,b = a.copy(), b.copy()
    a.name, b.name = 'a','b'
    m1 = get_cox_ph(surv, covariates=[a,b], formula='Surv(days, event) ~ a + b')
    if fisher_exact_test(a,b)['odds_ratio'] > 1:
        int_direction = 'both'
    else:
        int_direction = 'neither'
        
    int_var = 1.*(combine(a,b)==int_direction)
    int_var.name = 'interaction'
    m2 = get_cox_ph(surv, int_var)
    return pd.Series({'interaction': int_direction, 'p': LR_test(m2, m1)})

def interaction(a,b, surv):
    try:
        return _interaction(a,b, surv)
    except:
        return pd.Series(index=['interaction','p'])

def extract_chi2(full, reduced):
    '''
    Extract chi2 statstic of likelihood ratio test 
    on two R models. 
    '''
    full_ll = list(full.rx2('loglik'))
    reduced_ll = list(reduced.rx2('loglik'))
    chi2 = 2*full_ll[1] - 2*reduced_ll[1]
    return chi2

def get_interaction_simple(a,b, surv, int_direction='both'):
    '''
    Get test statistic (chi2 distributed) of interaction between 
    two event vectors.  
    '''
    a,b = a.copy(), b.copy()
    a.name, b.name = 'a','b'
    m1 = get_cox_ph(surv, covariates=[a,b], 
                    formula='Surv(days, event) ~ a + b')

    int_var = 1.*(combine(a,b)==int_direction)
    int_var.name = 'interaction'
    m2 = get_cox_ph(surv, int_var)
    chi2 = extract_chi2(m2, m1)
    return chi2

def get_interaction(a,b, surv, int_direction='both'):
    '''
    Get test statistic (chi2 distributed) of interaction between 
    two event vectors.  
    
    We define 3 models: 
        1) a + b
        2) a:b
        3) a + b + a:b
        
    We return the improvement of fit from 2 to 1 minus the 
    improvement of fit from 3 to 2. That is we want to capture
    as much of the information in the interaction term as possible.
    '''
    a,b = a.copy(), b.copy()
    a.name, b.name = 'a','b'
    m1 = get_cox_ph(surv, covariates=[a,b], 
                    formula='Surv(days, event) ~ a + b')
    int_var = 1.*(combine(a,b)==int_direction)
    int_var.name = 'interaction'
    m2 = get_cox_ph(surv, int_var)
    
    m3 = get_cox_ph(surv, combine(a,b))
    
    chi2_a = extract_chi2(m2, m1)
    chi2_b = extract_chi2(m3, m2)
    return chi2_a - chi2_b

def interaction_empirical_p(a, b, surv, num_perm=101):
    '''
    Calculate an empirical p-value for an interaction by sampling
    with replacement.  
    
    We first test if there is an improvement in model fit by 
    considering the interaction of the two events.  If so, we 
    then derive an empirical p-value. 
    '''
    a,b = match_series(a,b)
    if fisher_exact_test(a,b)['odds_ratio'] > 1:
        int_direction = 'both'
    else:
        int_direction = 'neither'
    r = get_interaction(a, b, surv)
    mat = np.array([np.random.permutation(a.index) for i in range(num_perm)])
    
    vec = {}
    for i,idx in enumerate(mat):
        a_p = pd.Series(list(a.ix[idx]), range(len(idx)))
        b_p = pd.Series(list(b.ix[idx]), range(len(idx)))
        surv_p = pd.DataFrame(surv.unstack().ix[a.index].as_matrix(), 
                              index=range(len(idx)), 
                              columns=['days','event']).stack()
        vec[i] = get_interaction(a_p, b_p, surv_p, int_direction)
    vec = pd.Series(vec).dropna()
    empirical_p = 1.*(len(vec) - sum(vec <= r)) / len(vec)
    return pd.Series({'p': empirical_p, 'interaction': int_direction})

def interaction_empirical_p_resample(a, b, surv, num_perm=101, check_first=True):
    '''
    Calculate an empirical p-value for an interaction by sampling
    with replacement.  
    
    We first test if there is an improvement in model fit by 
    considering the interaction of the two events.  If so, we 
    then derive an empirical p-value. 
    '''
    a,b = match_series(a,b)
    if fisher_exact_test(a,b)['odds_ratio'] > 1:
        int_direction = 'both'
    else:
        int_direction = 'neither'
    r = get_interaction(a, b, surv)
    if (r < 0) and (check_first is True):
        return pd.Series({'p': 1, 'interaction': int_direction})
    
    mat = np.random.choice(a.index, size=(num_perm, len(a.index)))
    
    vec = {}
    for i,idx in enumerate(mat):
        a_p = pd.Series(list(a.ix[idx]), range(len(idx)))
        b_p = pd.Series(list(b.ix[idx]), range(len(idx)))
        surv_p = pd.DataFrame(surv.unstack().ix[a.index].as_matrix(), 
                              index=range(len(idx)), 
                              columns=['days','event']).stack()
        vec[i] = get_interaction(a_p, b_p, surv_p, int_direction)
    vec = pd.Series(vec)
    
    empirical_p = 1.*(len(vec) - sum(vec <= r)) / len(vec)
    return pd.Series({'p': empirical_p, 'interaction': int_direction})
'''
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
'''

def get_interactions(df, cov_df, surv, test):
    binary = df[df.T.describe().ix['unique'] == 2]
    
    '''drop redundant features within a data-type'''
    s = {b for i,(a,v1) in enumerate(binary.iterrows())
           for j,(b,v2) in enumerate(binary.iterrows())
           if (i < j)
           and a[0] == b[0]
           and np.log2(fisher_exact_test(v1,v2)['odds_ratio']) > 4}
    binary = binary.ix[binary.index.diff(s)]

    n_tests = (len(binary) * (len(binary) - 1)) / 2
    s = pd.DataFrame({(a,b): interaction_empirical_p(v1,v2, surv, num_perm=101) 
                          for a,v1 in binary.iterrows()
                          for b,v2 in binary.iterrows()
                          if (a < b)
                          and fisher_exact_test(v1,v2).ix['p'] < (.05 / n_tests)
                          and fisher_exact_test(v1,v2).ix['odds_ratio'] != np.inf
                          and a[0] != b[0]}).T
    int_pairs =  s.ix[s.p < .1].sort('p')
    
    int_associations = {}
    for p,vals in int_pairs.iterrows():
        combo = combine(binary.ix[p[0]], binary.ix[p[1]])
        vec = combo == vals['interaction']
        int_associations[p] = test(vec, surv, cov_df) 
    int_associations = pd.DataFrame(int_associations).T
    return s, int_associations
