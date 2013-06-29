'''
Created on Jun 24, 2013

@author: agross
'''
import matplotlib.pylab as plt
import numpy as np
import pandas as pd

import rpy2.robjects as robjects

from pandas import Series, DataFrame
from pandas.rpy.common import convert_to_r_dataframe, convert_robj
from Processing.Tests import process_factors, get_formula

survival = robjects.packages.importr('survival')
base = robjects.packages.importr('base')

#TODO: fix, see HNSC_now


def get_surv_fit(surv, feature=None, covariates=None, interactions=None,
                 formula=None):
    if covariates is None:
        covariates = DataFrame(index=feature.index)
    c_real = covariates.ix[:,covariates.dtypes.isin([np.dtype(float), np.dtype(int)])]
    c_real = (c_real - c_real.mean()) / c_real.std()
    covariates[c_real.columns] = c_real
    clinical = covariates.join(surv.unstack()).dropna()
    clinical['days'] = clinical['days'] / 365.
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
    
    s = survival.survfit(fmla, df)
    summary = base.summary(s, times=robjects.r.c(5))
    #return summary
    res =  convert_robj(summary.rx2('table'))
    res = res.rename(index=lambda idx: idx.split('=')[1])
    res = res[['records','events','median','0.95LCL','0.95UCL']]
    res.columns = pd.MultiIndex.from_tuples([('Stats','# Patients'), ('Stats','# Events'), 
                           ('Median Survival', 'Median'), ('Median Survival', 'Lower'),
                           ('Median Survival', 'Upper')])
    idx = map(lambda s: s.replace('feature=',''), summary.rx2('strata').iter_labels())
    df = pd.DataFrame({d: list(summary.rx2(d)) for d in ['strata','surv','lower','upper']},
                      index=idx)
    res[('5y Survival', 'Lower')] = df['lower']
    res[('5y Survival', 'Surv')] = df['surv']
    res[('5y Survival', 'Upper')] = df['upper']
    return res

def survival_stat_plot(t, axs=None):
    if axs is None:
        fig = plt.figure(figsize=(6,1.5))
        ax = plt.subplot2grid((1,3), (0,0), colspan=2)
        ax2 = plt.subplot2grid((1,3), (0,2))
    else:
        ax, ax2 = axs
        fig = plt.gcf()
    for i,(idx,v) in enumerate(t.iterrows()):
        conf_int = v['Median Survival']
        median_surv = v[('Median Survival','Median')]
        if (v['Stats']['# Events']  / v['Stats']['# Patients']) < .5:
            median_surv = np.nanmin([median_surv, 6])
            conf_int['Upper'] = np.nanmin([conf_int['Upper'], 6])
        l = ax.plot(*zip(*[[conf_int['Lower'],i], [median_surv,i], [conf_int['Upper'],i]]), lw=3, ls='--', 
                    marker='o', dash_joinstyle='bevel')
        ax.scatter(median_surv,i, marker='s', s=100, color=l[0].get_color(), edgecolors=['black'], zorder=10,
                   label=idx)
    ax.set_yticks(range(len(t)))
    ax.set_yticklabels(['{} ({})'.format(idx, int(t.ix[idx]['Stats']['# Patients'])) 
                        for idx in t.index])
    ax.set_ylim(-.5, i+.5)
    ax.set_xlim(0,5)
    ax.set_xlabel('Median Survival (Years)')
    
    tt = t['5y Survival']
    b = (tt['Surv']).plot(kind='barh', ax=ax2, 
         color=[l.get_color() for l in ax.lines],
         xerr=[tt.Surv-tt.Lower, tt.Upper-tt.Surv], ecolor='black')
    ax2.set_xlabel('5Y Survival')
    ax2.set_xticks([0, .5, 1.])
    ax2.set_yticks([])
    fig.tight_layout()
    
def plot_5y(t, ax):
    tt = t['5y Survival']
    b = ax.bar(range(len(tt)), tt['Surv'],
         color=plt.rcParams['axes.color_cycle'],
         yerr=[tt.Surv-tt.Lower, tt.Upper-tt.Surv], ecolor='black')
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_position(('data', 0))
    #adjust_spines(ax,['left','bottom'])
    ax.spines['left'].set_position(('data', 0))
    
    ax.set_xticks(np.arange(len(tt))+.4)
    ax.set_xticklabels(['n={}'.format(int(t.ix[idx]['Stats']['# Patients']))
                        for idx in t.index], rotation=0)
    
    ax.set_ylabel('5Y Survival', position=(0.5,.66))
    ax.set_yticks([0, .5, 1.])
    #ax.yaxis.set_data_interval(0,1)
    ax.spines['left'].set_bounds(0, 1)
    ax.spines['bottom'].set_bounds(0, len(tt))
    
    for i,idx in enumerate(t.index):
        ax.annotate(idx[1], (i + .4, -.2), ha='center')
    
    ax.hlines(-.25, 0, 1.9)  
    ax.annotate(t.index.names[0], (1, -.35), ha='center')
    ax.annotate(t.index.names[1], (0, -.2), ha='right')
    ax.set_ylim(-.5,1)
    ax.set_xlim(-1,len(tt))
    ax.set_xlim(-2.,len(tt))
    
def plot_5y2(t, ax):
    tt = t['5y Survival']
    b = ax.bar(range(len(tt)), tt['Surv'],
         color=plt.rcParams['axes.color_cycle'][1],
         yerr=[tt.Surv-tt.Lower, tt.Upper-tt.Surv], ecolor='black')
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_position(('data', 0))
    #adjust_spines(ax,['left','bottom'])
    ax.spines['left'].set_position(('data', 0))
    
    ax.set_xticks(np.arange(len(tt))+.4)
    ax.set_xticklabels([])
    #ax.set_xticklabels(['{}'.format(int(t.ix[idx]['Stats']['# Patients']))
    #                    for idx in t.index], rotation=0)
    
    ax.set_ylabel('5Y Survival', position=(0.5,.66))
    ax.set_yticks([0, .5, 1.])
    #ax.yaxis.set_data_interval(0,1)
    ax.spines['left'].set_bounds(0, 1)
    ax.spines['bottom'].set_bounds(0, len(tt))
    
    for i,idx in enumerate(t.index):
        ax.annotate(str(int(t.ix[idx]['Stats']['# Patients'])), (i + .4, -.45), ha='center')
        ax.annotate({'Mut': 'X', 'WT': ''}[idx[1]], (i + .4, -.15), ha='center')
        ax.annotate({1: 'X', 0: ''}[idx[0]], (i + .4, -.3), ha='center')
    
    #ax.hlines(-.3, 0, 1.9)  
    ax.annotate('n= ', (0, -.45), ha='right')
    ax.annotate(t.index.names[0], (0, -.3), ha='right')
    ax.annotate(t.index.names[1], (0, -.15), ha='right')
    ax.set_ylim(-.5,1)
    ax.set_xlim(-2.5,len(tt))