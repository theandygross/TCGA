'''
Created on Jun 30, 2013

@author: agross
'''
import pandas as pd
import numpy as np
import scipy.stats as stats

import rpy2.robjects as robjects
from pandas.rpy.common import convert_to_r_dataframe, convert_robj
from Processing.Helpers import get_vec_type, bhCorrection

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

def cox_model_selection(fmla, df):
    '''
    Perform stepwise model selection on a cox model using r.step.
    '''
    s = survival.coxph(fmla, df)
    reduced =mass.stepAIC(s, direction='both', trace=0, k=2)
    fmla_reduced = reduced.rx2('formula')
    s_reduced = survival.coxph(fmla_reduced, df)
    return s_reduced

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
    if len(feature.unique()) == 1:
        return None
    if feature.dtype is np.dtype('bool'):
        feature = 1.*feature
    if len(feature.value_counts()) > 5:
        return None
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
    if type(cov) is type(None):
        cov = pd.DataFrame(index=feature.index)
    if type(cov) == pd.Series:
        cov = pd.concat([cov], axis=1)
    elif type(cov) == list:
        assert map(type, cov) == ([pd.Series] * len(cov))
        cov = pd.concat(cov, axis=1)
    c_real = cov.ix[:,cov.dtypes.isin([np.dtype(float), np.dtype(int)])]
    c_real = (c_real - c_real.mean()) / c_real.std()
    cov[c_real.columns] = c_real
    df = cov.join(surv.unstack()).dropna()
    df['days'] = df['days'] / 365.
    df = df.groupby(level=0).first()
    df, factors = process_factors(df, feature, list(cov.columns))
    df = df[factors + ['days','event']]
    df = df.dropna()
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
    df, factors = process_covariates(surv, feature, covariates)
    if formula is None:
        formula = get_formula(factors, interactions)
        fmla = robjects.Formula(formula)
        s = cox_model_selection(fmla, df)
    else:
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
    s = get_cox_ph(surv, feature, covariates, formula, interactions,
                   get_model, print_desc)
    if s is None:
        return
    
    results = convert_robj(base.summary(s).rx2('coefficients'))
    
    def set_null_model():
        #patients = feature.dropna().index
        null_model = get_cox_ph(surv, covariates=covariates, 
                                feature=feature.map(lambda s: 1), 
                                interactions=interactions)
        get_cox_ph_ms.null_model = null_model
        get_cox_ph_ms.params = surv.name, str(covariates)
      
    # check if we need to recompute null model
    has_null = hasattr(get_cox_ph_ms, 'null_model') or (null_model)
    if has_null != True:
        set_null_model()
        
    recalc = get_cox_ph_ms.params !=  surv.name, str(covariates)
    if recalc:
        set_null_model()
    
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
                 formula=None):
    df, factors = process_covariates(surv, feature, covariates)
    if formula is None:
        fmla = get_formula(factors, interactions)
        fmla = robjects.Formula(fmla)
    else:
        fmla = robjects.Formula(formula)
    
    s = survival.survfit(fmla, df)
    summary = base.summary(s, times=robjects.r.c(5))
    #return summary
    res =  convert_robj(summary.rx2('table'))
    res = res.rename(index=lambda idx: idx.split('=')[1])
    res = res[['records','events','median','0.95LCL','0.95UCL']]
    res.columns = pd.MultiIndex.from_tuples([('Stats','# Patients'), 
                                             ('Stats','# Events'), 
                                             ('Median Survival', 'Median'), 
                                             ('Median Survival', 'Lower'),
                                             ('Median Survival', 'Upper')])
    idx = map(lambda s: s.replace('feature=',''), 
              summary.rx2('strata').iter_labels())
    df = pd.DataFrame({d: list(summary.rx2(d)) for d in 
                       ['strata','surv','lower','upper']},
                      index=idx)
    res[('5y Survival', 'Surv')] = df['surv']
    res[('5y Survival', 'Lower')] = df['lower']
    res[('5y Survival', 'Upper')] = df['upper']
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
    Perform univariate Cox regression.
    '''
    fmla = 'Surv(days, event) ~ feature'
    model = get_cox_ph(surv, feature, formula=fmla)
    lr = model[3][0]
    p = stats.chi2.sf(lr, 1)
    return pd.Series({'LR': lr, 'p': p})

def cox_screen(df, surv):
    rr = df.apply(cox, args=(surv,), axis=1)
    rr['q'] = bhCorrection(rr.p)
    rr = rr.sort('p')
    return rr

def lr_screen(df, surv):
    rr = df.astype(float).apply(log_rank, args=(surv,), axis=1)
    rr['q'] = bhCorrection(rr.p)
    rr = rr.sort('p')
    return rr