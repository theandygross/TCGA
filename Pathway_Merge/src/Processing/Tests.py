'''
Created on Jan 8, 2013

@author: agross
'''
import pickle as pickle
import pandas as pd
import rpy2.robjects as robjects

from scipy.stats import f_oneway, fisher_exact, pearsonr, chi2, kruskal
from numpy import nan, dtype

from pandas import Series, DataFrame
from pandas.rpy.common import convert_to_r_dataframe, convert_robj

from Processing.Helpers import match_series, get_vec_type
from Processing.Helpers import make_path_dump, bhCorrection

survival = robjects.packages.importr('survival')
base = robjects.packages.importr('base')

robjects.r.options(warn=-1);
zz = robjects.r.file("all.Rout", open="wt")
robjects.r.sink(zz, type='message')

def cox_model_selection(fmla, df):
    s = survival.coxph(fmla, df)
    reduced = robjects.r.step(s, direction='both', trace=0);
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
        covariates = [cov for cov in covariates if cov in clinical]
    if type(hit_vec) != type(None):
        hit_vec.name = 'feature'
        factors = ['feature'] + covariates
        df = clinical.join(hit_vec)
    else:
        factors = covariates
        df = clinical
    return df, factors

def get_formula(factors, get_interactions=True):
    if len(factors) > 1:
        if get_interactions:
            interactions = ' + '.join(('*'.join([a,b]) for a in factors 
                                      for b in factors if a<b))
        else:
            interactions = ' + '.join(factors)
    elif len(factors) == 1: 
        interactions = factors[0]
    else:
        interactions = '1'
    fmla = 'Surv(days, event) ~ {}'.format(interactions)
    return fmla
                        
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
    if covariates is None:
        covariates = DataFrame(index=feature.index)
    #covariates = (covariates - covariates.mean()) / covariates.std()
    c_real = covariates.ix[:,covariates.dtypes.isin([dtype(float), dtype(int)])]
    c_real = (c_real - c_real.mean()) / c_real.std()
    covariates[c_real.columns] = c_real
    clinical = covariates.join(surv.unstack()).dropna()
    clinical = clinical.groupby(level=0).first()
    
    df, factors = process_factors(clinical, feature, list(covariates.columns))
    df = df[factors + ['days','event']]
    df = df.dropna()
    
    df = convert_to_r_dataframe(df)
    if formula is None:
        fmla = get_formula(factors, interactions)
        fmla = robjects.Formula(fmla)
    else:
        fmla = robjects.Formula(formula)
    
    if formula is not None:
        s = survival.coxph(fmla, df)
    else:
        s = cox_model_selection(fmla, df)
    results = convert_robj(base.summary(s).rx2('coefficients'))
   
    if return_val == 'p':
        return results.ix['feature','Pr(>|z|)']
    
    elif return_val == 'p_haz':
        try:
            hazzard = results.ix['feature','exp(coef)']
            p = results.ix['feature','Pr(>|z|)']
        except:
            hazzard = nan
            p = nan
        return Series({'p': p, 'hazzard': hazzard})
    
    elif return_val == 'model':
        return s
    elif return_val in ['LR', 'LR_p']:
        if (((not hasattr(get_cox_ph_ms, 'null_model')) and (null_model == None))
            or (get_cox_ph_ms.params !=  surv.name, list(covariates.columns))):
            '''Save some time by non-recomputing null model'''
            patients = feature.dropna().index
            null_model = get_cox_ph_ms(surv, covariates=covariates.ix[patients], 
                                       return_val='model')
            get_cox_ph_ms.null_model = null_model
            get_cox_ph_ms.params = surv.name, list(covariates.columns)
        else:
            null_model = get_cox_ph_ms.null_model
        LR_p = LR_test(s, null_model)
        
        if type(results) == DataFrame and 'feature' in results.index:
            coef_p = results.ix['feature','Pr(>|z|)']
            hazzard = results.ix['feature','exp(coef)']
        else:
            coef_p, hazzard = nan, nan
    if return_val == 'LR':
        f = str(s.rx2('formula'))
        results = Series({'LR': LR_p, 'feature_p': coef_p, 'hazzard': hazzard, 
                          'fmla': f})
        if get_model:
            results['model'] = s
        return results
    if return_val == 'LR_p':
        return LR_p
    
def get_cox_ph(clinical, hit_vec=None, covariates=[], time_var='days',
               event_var='censored'):
    p = get_cox_ph_ms(clinical, hit_vec, covariates, time_var, event_var, 'LR_p')
    return p


def anova(hit_vec, response_vec, min_size=5):
    '''
    Wrapper to do a one way anova on pandas Series
    ------------------------------------------------
    hit_vec: Series of labels
    response_vec: Series of measurements
    '''
    if hit_vec.value_counts().min < min_size:
        return nan
    hit_vec, response_vec = match_series(hit_vec, response_vec)
    return f_oneway(*[response_vec[hit_vec == num] for num in 
                      hit_vec.unique()])[1]

def fisher_exact_test(hit_vec, response_vec):
    '''
    Wrapper to do a fischer's exact test on pandas Series
    ------------------------------------------------
    hit_vec: Series of labels (boolean, or (0,1))
    response_vec: Series of measurements (boolean, or (0,1))
    '''
    #assert ((len(hit_vec.unique()) <= 2) and 
    #        (len(response_vec.unique()) <= 2)) 
    cont_table = pd.crosstab(hit_vec, response_vec)
    if (cont_table.shape != (2,2)):
        return 1
    return fisher_exact(cont_table)[1]

def kruskal_p(hit_vec, response_vec, min_size=5):
    '''
    Wrapper to do a one way anova on pandas Series
    ------------------------------------------------
    hit_vec: Series of labels
    response_vec: Series of measurements
    '''
    try:
        hit_vec, response_vec = match_series(hit_vec, response_vec)
        return kruskal(*[response_vec[hit_vec == num] for num in 
                          hit_vec.unique()])[1]
    except:
        return nan
                      
def pearson_p(a,b):
    '''
    Find pearson's correlation and return p-value.
    ------------------------------------------------
    a, b: Series with continuous measurements
    '''
    a,b = match_series(a.dropna(), b.dropna())
    _,p = pearsonr(a,b)
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
        
class AnovaTest(object):
    def __init__(self, response_vec):
        self.response_vec = response_vec
        self.full_test = self.get_full_test()
    
    def get_full_test(self):
        def test(hit_vec):
            hit_vec, response_vec = match_series(hit_vec, self.response_vec)
            res =  f_oneway(*[response_vec[hit_vec == num] for num in 
                      hit_vec.unique()])
            return Series({'stat': res[0], 'p': res[1]})
        return test
    def check_feature(self, feature):
        return feature.value_counts().min > 5.
    
def run_feature_matrix(df, test, fp_cutoff=.5):
    df = df.ix[df.apply(test.check_feature, 1)]
    if hasattr(test, 'first_pass'):
        fp = df.apply(test.first_pass, 1)
        df = df[fp.p < fp_cutoff]
    full = df.apply(test.full_test, 1)
    res = pd.concat([full[['LR', 'fmla']], fp], keys=['Full', 'Univariate'], axis=1)
    if type(res.index[0]) == tuple: #pandas bug
                res.index = pd.MultiIndex.from_tuples(res.index, names=df.index.names) 
    res = res.join(Series(bhCorrection(res[('Full', 'LR')], n=len(fp)), 
                          name=('Full','LR_q')))
    res = res.join(Series(bhCorrection(res[('Univariate', 'p')], n=len(fp)), 
                          name=('Univariate','q')))
    return res.sort_index(axis=1).sort(columns=[('Full', 'LR')])

class TestResult(object):
    def __init__(self, test, dataset, results):
        self.test = test
        self.name = test.name
        self.results = results.sort(columns=[test.statistic])
        self.p_values = results[test.statistic].order()
        self.features_tested = dataset.features.index
        self._dataset_path = dataset.path
        self.path = '/'.join([self._dataset_path, 'Tests', self.name, ''])
    
    def get_dataset(self):
        obj = pickle.load(open(self._dataset_path + '/DataObject.p', 'rb'))
        return obj
    
    def save(self):
        if hasattr(self.test, 'first_pass'):
            del self.test.first_pass
        if hasattr(self.test, 'full_test'):
            del self.test.full_test
        make_path_dump(self, self.path + '/TestResult.p')
        self.results.to_csv(self.path + '/TestResults.csv')