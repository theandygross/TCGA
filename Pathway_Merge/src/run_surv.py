#! /usr/bin/env python

import sys as sys
import os as os

import pickle as pickle
import pandas as pd

from Processing.Tests import SurvivalTest, TestResult
from Processing.Tests import run_feature_matrix
from Reports.Figures import draw_survival_curves_model, draw_survival_curves

def run_surv(cancer, data, covariates, survival_test, dry_run=False):
    clinical = cancer.load_clinical()
    global_vars = cancer.load_global_vars()
    
    try:
        surv = clinical.survival[survival_test]
    except AttributeError:
        surv = clinical.survival.survival_5y
        
    cov_df = global_vars.join(clinical.clinical, how='outer')
    cov_df = cov_df[covariates]
    remerge = lambda s: '__'.join(s) if type(s) != str else s
    cov_df = cov_df.rename(columns=remerge)
    test = SurvivalTest(surv, cov_df)
    test.name = survival_test
    
    events = run_feature_matrix(data.features, test)
    r = TestResult(test, data, events)
    if dry_run:
        return r
    
    r.save()
    thresholded = r.results.ix[r.p_values < .01]
    fig_path = '/'.join([r.path, 'Figures',''])
    if not os.path.isdir(fig_path):
        os.makedirs(fig_path)
        
    print data.data_type
    for f in thresholded.index:
        feature = data.features.ix[f]
        if type(feature.name) != str:
            feature.name = str(feature.name)
        filename = fig_path + feature.name + '.png'
        if data.data_type.startswith('CN'):
            labels = pd.Series({-2: 'Homozygous Deletion', -1: 'Deletion', 
                             0: 'Normal', 1: 'Amp', 2: 'High Amp'})
            colors = pd.Series({-2: 'black', -1: 'purple', 0: 'blue', 
                             1: 'orange', 2: 'red'})
            vals = sorted(feature.unique())
            draw_survival_curves(feature, test.surv, colors=colors[vals].tolist(), 
                                 labels=labels[vals].tolist(), filename=filename,
                                 ann='p')
        elif data.data_type in ['mRNASeq', 'RPPA', 'Methylation']:
            vec = (feature - feature.mean()) / feature.std()
            vec = ((vec > vec.quantile(.75)).astype(int) - 
                   (vec < vec.quantile(.25)).astype(int)).astype(float)
            draw_survival_curves(vec, test.surv, colors=['blue','orange','red'], 
                           labels=['Bottom 25%', 'Normal', 'Top 25%'],
                           filename=filename)
                       
                
        else:
            draw_survival_curves_model(feature, r.test, ann='p', 
                                       filename=filename)
    return r

if __name__ == "__main__":
    report_path = sys.argv[1]
    cancer_type = sys.argv[2]
    data_type = sys.argv[3]
    covariates = sys.argv[4].split(',')
    parse = lambda c: tuple(c.split('__')) if '__' in c else c
    covariates = map(parse, covariates)
    survival_test = sys.argv[5]   
    
    run = pickle.load(open('/'.join([report_path, 'RunObject.p']), 'rb'))
    cancer = run.load_cancer(cancer_type) 
    data = cancer.load_data(data_type)
    try:
        data.uncompress()
    except:
        do_nothing = True
    run_surv(cancer, data, covariates, survival_test)
    sys.exit(0)
