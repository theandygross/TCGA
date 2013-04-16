'''
Created on Apr 7, 2013

@author: agross
'''
from Processing.Helpers import get_vec_type, to_quants
from Processing.Tests import get_cox_ph_ms

import pandas as pd
import matplotlib.pylab as plt
import pandas.rpy.common as com
import rpy2.robjects as robjects 

survival = robjects.packages.importr('survival')
base = robjects.packages.importr('base')

def get_markers(censoring, survival):
    '''
    Get locations for markers in KM plot.
    censoring is a list of censoring times.
    survival is a time-series of survival values
    '''
    markers = []
    for cc in censoring:
        d = (pd.Series(survival.index, survival.index, dtype=float) - cc)
        t = d[d<=0].idxmax()
        markers += [(cc, survival[t])]
    return markers

def draw_survival_curves_mpl(fit, ax=None, title=None, colors=None):
    '''
    Takes an R survfit.
    '''
    if ax is None:
        _, ax = plt.subplots(1,1, figsize=(4,3))
    s = base.summary(fit)
    tab = pd.DataFrame({v: s.rx2(v) for v in s.names 
                                    if len(s.rx2(v)) == len(s.rx2('time'))}, 
                       index=s.rx2('time'))
    call = com.convert_robj(fit.rx2('call')[2])
    groups = sorted(call.feature.unique())
    
    if len(tab.strata.unique()) == 1: #TODO: fix for more than one curve
        groups = call.feature.value_counts().index
    
    for i,group in enumerate(groups):
        censoring = call[(call.event==0) * (call.feature==group)].days
        surv = tab[tab.strata==(i+1)].surv
        surv = surv.set_value(0, 1).sort_index()
        if surv.index[-1] < censoring.max():
            surv = surv.set_value(censoring.max(), surv.iget(-1)).sort_index()

        censoring_pos = get_markers(censoring, surv)
        ax.step(surv.index, surv, lw=3, where='post', alpha=.7, label=group)
        if colors is not None:
            ax.lines[-1].set_color(colors[i])
        ax.scatter(*zip(*censoring_pos), marker='|', s=80, 
                   color=ax.lines[-1].get_color())
        
    ax.set_ylim(0,1.05)
    ax.set_xlim(0, max(surv.index)*1.05)
    ax.legend(loc='best')
    ax.set_ylabel('Survival')
    ax.set_xlabel('Years')
    if title:
        ax.set_title(title)
        
def process_feature(feature, q, std):
    if (get_vec_type(feature) == 'real') and (len(feature.unique()) > 10):
        feature = to_quants(feature, q=q, std=std, labels=True)
    return feature

def draw_survival_curve(feature, surv, q=.25, std=None, **args):
    feature = process_feature(feature, q, std)
    fmla = robjects.Formula('Surv(days, event) ~ feature')           
    m = get_cox_ph_ms(surv, feature, return_val='model', formula=fmla)
    r_data = m.rx2('call')[2]
    s = survival.survdiff(fmla, r_data)
    p = str(s).split('\n\n')[-1].strip().split(', ')[-1]
    draw_survival_curves_mpl(survival.survfit(fmla, r_data), **args)
    
def draw_survival_curves(feature, surv, assignment=None):
    if assignment is None:
        draw_survival_curve(feature, surv)
        return
    num_plots = len(assignment.unique())
    fig, axs = plt.subplots(1,num_plots, figsize=(num_plots*4,3), sharey=True)
    for i,(l,s) in enumerate(feature.groupby(assignment)):
        draw_survival_curve(s, surv, ax=axs[i], 
                            title='{} = {}'.format(assignment.name,l))