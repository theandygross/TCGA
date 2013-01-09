'''
Created on Jan 8, 2013

@author: agross
'''
from scipy.stats import chi2
from numpy import nan

import rpy2.robjects as robjects
from pandas.rpy.common import convert_to_r_dataframe, convert_robj
from pandas import DataFrame

survival = robjects.packages.importr('survival')
base = robjects.packages.importr('base')

robjects.r.options(warn=-1);
zz = robjects.r.file("all.Rout", open="wt")
robjects.r.sink(zz, type='message')

def cox_model_selection(fmla, df):
    s = survival.coxph(fmla, df)
    reduced = robjects.r.step(s, trace=0);
    fmla_reduced = reduced.rx2('formula')
    s_reduced = survival.coxph(fmla_reduced, df)
    return s_reduced

def LR_test(full, reduced):
    full_ll = list(full.rx2('loglik'))
    reduced_ll = list(reduced.rx2('loglik'))
    assert full_ll[0] == reduced_ll[0]
    if len(full_ll) == 1:
        return 1.
    full_df = len(full.rx2('coefficients')) 
    if len(reduced_ll) == 1:
        return chi2.sf(2*full_ll[1] - 2*full_ll[0], full_df)
    reduced_df = len(reduced.rx2('coefficients'))
    df = max(full_df - reduced_df, 1)
    return chi2.sf(2*full_ll[1] - 2*reduced_ll[1], df)

def process_factors(clinical, hit_vec=None, covariates=[]):
    if not all([cov in clinical for cov in covariates]):
        missing = [cov for cov in covariates if cov not in clinical]
        covariates = [cov for cov in covariates if cov in clinical]
    if type(hit_vec) != type(None):
        hit_vec.name = 'pathway'
        factors = ['pathway'] + covariates
        df = clinical.join(hit_vec)
    else:
        factors = covariates
        df = clinical
    return df, factors
                        
def get_cox_ph_ms(clinical, hit_vec=None, covariates=[], time_var='days',
                  event_var='censored', return_val='p', null_model=None):
    '''
    Fit a cox proportial hazzards model to the data.
    Returns a p-value on the hit_vec coefficient. 
    ---------------------------------------------------
    clinical: DataFrame of clinical variables
    hit_vec: vector of labels to test against
    covariates: names of covariates in the cox model,
                (must be columns in clinical DataFrame)
    '''
    df, factors = process_factors(clinical, hit_vec, ['age', 'rate'])
    df = df[factors + [time_var, event_var]]
    df = df.dropna()
    df[factors] = (df[factors] - df[factors].mean()) / df[factors].std()
    df = convert_to_r_dataframe(df)
    if len(factors) > 1:
        interactions = ' + '.join(('*'.join([a,b]) for a in factors 
                                   for b in factors if a<b))
    else: 
        interactions = factors[0]
    fmla = 'Surv(' + time_var + ', ' + event_var + ') ~ ' + interactions
    fmla = robjects.Formula(fmla)
    
    try:
        s = cox_model_selection(fmla, df)
        results = convert_robj(base.summary(s).rx2('coefficients'))
    except robjects.rinterface.RRuntimeError:
        return 1.23
    
    if return_val == 'p':
        return results.ix['pathway','Pr(>|z|)']
    elif return_val == 'coef':
        return results
    elif return_val == 'model':
        return s
    elif return_val in ['LR', 'LR_p']:
        if (((not hasattr(get_cox_ph_ms, 'null_model')) and (null_model == None))
            or (get_cox_ph_ms.params != time_var, event_var, covariates)):
            '''Save some time by non-recomputing null model'''
            patients = hit_vec.dropna().index
            null_model = get_cox_ph_ms(clinical.ix[patients], return_val='model',
                                       covariates=covariates, time_var=time_var, 
                                       event_var=event_var)
            get_cox_ph_ms.null_model = null_model
            get_cox_ph_ms.params = time_var, event_var, covariates
        else:
            null_model = get_cox_ph_ms.null_model
        LR_p = LR_test(s, null_model)
        if type(results) == DataFrame and 'pathway' in results.index:
            coef_p = results.ix['pathway','Pr(>|z|)']
            hazzard = results.ix['pathway','exp(coef)']
        else:
            coef_p, hazzard = nan, nan
    if return_val == 'LR':
        return {'LR': LR_p, 'pathway_p': coef_p, 'hazzard': hazzard, 'model': s}
    if return_val == 'LR_p':
        return coef_p if LR_p < .1 else nan
    
def get_cox_ph(clinical, hit_vec=None, covariates=[], time_var='days',
               event_var='censored'):
    p = get_cox_ph_ms(clinical, hit_vec, covariates, time_var, event_var, 'LR_p')
    return p