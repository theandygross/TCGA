'''
Created on Aug 14, 2013

@author: agross
'''
from Processing.Helpers import to_quants, screen_feature
from Stats.Scipy import chi2_cont_test

import numpy as np
import pandas as pd

from sklearn.cross_validation import train_test_split
from sklearn.grid_search import GridSearchCV
from sklearn.metrics import auc_score, precision_score
from sklearn.svm import SVC

def SVC_fill_old(feature, df):
    gg = df.apply(lambda s: to_quants(s, std=1) > 0)
    diff = screen_feature(feature, chi2_cont_test, gg)
    dd = diff[diff.p < .05]
    
    pats = gg.columns.intersection(feature.index)
    mat = gg.ix[dd.index]
    X = mat.ix[:, pats].T.as_matrix()
    Y = (feature)
    Y = np.array(Y.ix[pats])
    
    X_train, X_test, y_train, y_test = train_test_split(
        X, Y, test_size=0.35, random_state=5796543)
    
    params = [{'kernel': ['rbf'], 'gamma': [0, .1, .05, .01, 1e-3, 1e-4, 1e-5],
               'C': [.1, 1, 10, 100, 1000], 'class_weight': ['auto']},
              {'kernel': ['linear'], 'C': [1, 10, 100, 1000],
               'class_weight': ['auto']},
              {'kernel': ['poly'], 'C': [1, 10, 100, 1000],
               'class_weight': ['auto']}]
    
    clf = GridSearchCV(SVC(C=1), params, score_func=auc_score)
    clf.fit(X_train, y_train, cv=5);
    best = clf.best_estimator_
    auc = clf.score(X, Y)
    
    mat_all = gg.ix[mat.index].T.as_matrix()
    inferred = best.predict(mat_all)
    inferred = pd.Series(inferred, index=gg.columns)
    fun = pd.Series(best.decision_function(mat_all)[:, 0], mat.columns)
    f = feature.copy()
    f = f.ix[inferred.index]
    f[f.isnull()] = inferred[f.isnull()]
    filled_feature = f.astype(float)
    return {'auc': auc, 'model': best, 'decision_function': fun,
            'inferred_values': inferred, 'filled_feature': filled_feature} 
    
def SVC_fill(feature, df, metric='auc', test_size=.3):
    '''
    ###Support Vector Inference 
    * Using [SVC](http://scikit-learn.org/stable/modules/generated/sklearn.svm.SVC.html) 
    function from [Scikit Learn Package](http://scikit-learn.org/stable/) 
    * Features are the binarized differential expression vectors  
      * Have high change in expression from tumor to normal 
      * Thresholded at 1 standard deviation over the mean to reduce overfitting
    * Parameters are fit using cross validation, optimizing for AUC score 
      * I try linear, RBF, and polynomial kernels under a variety of parameters 
      * The best model in cross validation is fit on the entire dataset 
    * Missing values are filled in based on the model prediction
    '''
    gg = df.apply(lambda s: to_quants(s, std=1) > 0, axis=1)
    mat = gg
    
    pats = gg.columns.intersection(feature.index)
    X = mat.ix[:, pats].T.as_matrix()
    Y = (feature)
    Y = np.array(Y.ix[pats])
    
    X_train, X_test, y_train, y_test = train_test_split(
        X, Y, test_size=test_size, random_state=5796503)
    
    tuned_parameters = [{'kernel': ['rbf'], 'gamma': [0, .1, .05, .01, 1e-3, 1e-4, 1e-5],
                         'C': [.1, 1, 10, 100, 1000], 'class_weight': ['auto']},
                        {'kernel': ['linear'], 'C': [1, 10, 100, 1000], 'class_weight': ['auto']},
                        {'kernel': ['poly'], 'C': [1, 10, 100, 1000], 'class_weight': ['auto']}]
    
    if metric == 'auc':
        metric = auc_score
    elif metric == 'precision':
        metric = precision_score
    
    clf = GridSearchCV(SVC(C=1), tuned_parameters, score_func=metric)
    clf.fit(X_train, y_train, cv=5);
    best = clf.best_estimator_
    auc = clf.score(X, Y)
    
    mat_all = gg.ix[mat.index].T.as_matrix()
    inferred = best.predict(mat_all)
    inferred = pd.Series(inferred, index=gg.columns)
    fun = pd.Series(best.decision_function(mat_all)[:, 0], mat.columns)
    f = feature.copy()
    f = f.ix[inferred.index]
    f[f.isnull()] = inferred[f.isnull()]
    filled_feature = f.astype(float)
    return {'auc': auc, 'model': best, 'decision_function': fun,
            'inferred_values': inferred, 'filled_feature': filled_feature} 
