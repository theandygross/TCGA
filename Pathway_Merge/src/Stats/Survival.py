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

robjects.r.options(warn=-1);
zz = robjects.r.file("all.Rout", open="wt")
robjects.r.sink(zz, type='message')


def cox_model_selection(fmla, df):
    '''
    Perform stepwise model selection on a cox model using r.step.
    '''
    s = survival.coxph(fmla, df)
    reduced = robjects.r.step(s, direction='both', trace=0)
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

def log_rank(feature, surv):
    fmla = robjects.Formula('Surv(days, event) ~ feature')
    m = get_cox_ph_ms(surv, feature, return_val='model', formula=fmla)
    r_data = m.rx2('call')[2]
    s = survival.survdiff(fmla, r_data)
    p = stats.chi2.sf(s.rx2('chisq')[0], 1)
    return pd.Series({'chi2': s.rx2('chisq')[0], 'p': p})

def get_formula(factors, get_interactions=True):
    if len(factors) > 1:
        interactions = ' + '.join(factors)
        if get_interactions:
            interactions += ' + '
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

def process_covariates(surv, feature=None, covariates=None):
    if covariates is None:
        covariates = pd.DataFrame(index=feature.index)
    c_real = covariates.ix[:,covariates.dtypes.isin([np.dtype(float), 
                                                     np.dtype(int)])]
    c_real = (c_real - c_real.mean()) / c_real.std()
    covariates[c_real.columns] = c_real
    df = covariates.join(surv.unstack()).dropna()
    df['days'] = df['days'] / 365.
    df = df.groupby(level=0).first()
    df, factors = process_factors(df, feature, list(covariates.columns))
    df = df[factors + ['days','event']]
    df = df.dropna()
    df = convert_to_r_dataframe(df)
    return df, factors

def get_cox_ph_ms(surv, feature=None, covariates=None, return_val='p', 
                  null_model=None, formula=None, get_model=False,
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
    df, factors = process_covariates(surv, feature, covariates)
    if formula is None:
        formula = get_formula(factors, interactions)
        fmla = robjects.Formula(formula)
        s = cox_model_selection(fmla, df)
    else:
        fmla = robjects.Formula(formula)
        s = survival.coxph(fmla, df)
  
    results = convert_robj(base.summary(s).rx2('coefficients'))
      
    if return_val == 'model':
        return s
    
    elif return_val == 'model_desc':
        desc = '\n\n'.join(str(s).split('\n\n')[-2:])
        print desc
        return desc
    
    elif return_val in ['LR', 'LR_p']:
        # check if we need to recompute null model
        has_null = hasattr(get_cox_ph_ms, 'null_model') or (null_model)
        recalc = get_cox_ph_ms.params !=  surv.name, list(covariates.columns)
        if (has_null == False) or (recalc == True): 
            patients = feature.dropna().index
            null_model = get_cox_ph_ms(surv, covariates=covariates.ix[patients], 
                                       return_val='model')
            get_cox_ph_ms.null_model = null_model
            get_cox_ph_ms.params = surv.name, list(covariates.columns)
        else:
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
        if get_model:
            results['model'] = s
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
    res[('5y Survival', 'Lower')] = df['lower']
    res[('5y Survival', 'Surv')] = df['surv']
    res[('5y Survival', 'Upper')] = df['upper']
    return res
    
def get_cox_ph(clinical, hit_vec=None, covariates=[], time_var='days',
               event_var='censored'):
    p = get_cox_ph_ms(clinical, hit_vec, covariates, time_var, event_var, 
                      'LR_p')
    return p
    
    
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