'''
Created on Apr 7, 2013

@author: agross
'''
from Processing.Helpers import get_vec_type, to_quants
from Stats.Survival import get_cox_ph

import pandas as pd
import matplotlib.pylab as plt
import pandas.rpy.common as com
import rpy2.robjects as robjects 
from Stats.Survival import get_surv_fit

import numpy as np

survival = robjects.packages.importr('survival')
base = robjects.packages.importr('base')

colors_global = plt.rcParams['axes.color_cycle']*10

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

def draw_survival_curves_mpl(fit, ax=None, title=None, colors=None, ms=80, alpha=1):
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
    
    groups = robjects.r.sort(robjects.r.c(*call.feature.unique()))
    
    if 'strata' not in tab:
        groups = [0]
        tab['strata'] = 1
    elif len(tab.strata.unique()) != len(groups):
        gg = list(call[call.event > 0].feature.unique())
        gg = [g for g in groups if g in gg]
        bg = [g for g in groups if g not in gg]
        groups = gg + bg
           
    for i,group in enumerate(groups):
        censoring = call[(call.event==0) * (call.feature==group)].days
        surv = tab[tab.strata==(i+1)].surv
        surv = surv.set_value(0., 1.).sort_index()  #Pandas bug, needs to be float
        if surv.index[-1] < censoring.max():
            surv = surv.set_value(censoring.max(), surv.iget(-1)).sort_index()

        censoring_pos = get_markers(censoring, surv)
        ax.step(surv.index, surv, lw=3, where='post', alpha=alpha, label=group)
        if colors is not None:
            try:
                '''fix for R-Python str-to-int conversion'''
                color = colors[group]
            except:
                color = colors[i]
            ax.lines[-1].set_color(color)
        if len(censoring_pos) > 0:
            ax.scatter(*zip(*censoring_pos), marker='|', s=ms, 
                       color=ax.lines[-1].get_color())
        
    ax.set_ylim(0,1.05)
    #ax.set_xlim(0, max(surv.index)*1.05)
    ax.set_xlim(0, max(call.days)*1.05)
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
    m = get_cox_ph(surv, feature)
    r_data = m.rx2('call')[2]
    #s = survival.survdiff(fmla, r_data)
    #p = str(s).split('\n\n')[-1].strip().split(', ')[-1]
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
        
def survival_stat_plot(t, upper_lim=5, axs=None, colors=None):
    '''
    t is the DataFrame returned from a get_surv_fit call.
    '''
    if axs is None:
        fig = plt.figure(figsize=(6,1.5))
        ax = plt.subplot2grid((1,3), (0,0), colspan=2)
        ax2 = plt.subplot2grid((1,3), (0,2))
    else:
        ax, ax2 = axs
        fig = plt.gcf()
    if colors is None:
        colors = colors_global
    for i,(idx,v) in enumerate(t.iterrows()):
        conf_int = v['Median Survival']
        median_surv = v[('Median Survival','Median')]
        if (v['Stats']['# Events']  / v['Stats']['# Patients']) < .5:
            median_surv = np.nanmin([median_surv, 20])
            conf_int['Upper'] = np.nanmin([conf_int['Upper'], 20])
        l = ax.plot(*zip(*[[conf_int['Lower'],i], [median_surv,i], [conf_int['Upper'],i]]), lw=3, ls='--', 
                    marker='o', dash_joinstyle='bevel', color=colors[i])
        ax.scatter(median_surv,i, marker='s', s=100, color=l[0].get_color(), edgecolors=['black'], zorder=10,
                   label=idx)
    ax.set_yticks(range(len(t)))
    ax.set_yticklabels(['{} ({})'.format(idx, int(t.ix[idx]['Stats']['# Patients'])) 
                        for idx in t.index])
    ax.set_ylim(-.5, i+.5)
    ax.set_xlim(0,upper_lim)
    ax.set_xlabel('Median Survival (Years)')
    
    tt = t['5y Survival']
    b = (tt['Surv']).plot(kind='barh', ax=ax2, 
         color=[l.get_color() for l in ax.lines],
         xerr=[tt.Surv-tt.Lower, tt.Upper-tt.Surv], ecolor='black')
    ax2.set_xlabel('5Y Survival')
    ax2.set_xticks([0, .5, 1.])
    ax2.set_yticks([])
    fig.tight_layout()
    
def survival_and_stats(feature, surv, upper_lim=5, axs=None, figsize=(7,5), title=None, 
                       order=None, colors=None, **args):
    if axs is None:
        fig = plt.figure(figsize=figsize)
        ax1 = plt.subplot2grid((3,3), (0,0), colspan=3, rowspan=2)
        ax2 = plt.subplot2grid((3,3), (2,0), colspan=2)
        ax3 = plt.subplot2grid((3,3), (2,2))
    else:
        ax1, ax2, ax3 = axs
        fig = plt.gcf()
    if feature.dtype != str:
        feature = feature.astype(str)
    if colors is None:
        colors = colors_global
    
    t = get_surv_fit(surv, feature)
    if order is None:
        t = t.sort([('5y Survival','Surv')], ascending=True)
    else:
        t = t.ix[order]
    survival_stat_plot(t, axs=[ax2, ax3], upper_lim=upper_lim, colors=colors)
    r = pd.Series({s:i for i,s in enumerate(t.index)})
    color_lookup = {c: colors[i % len(colors)] for i,c in enumerate(t.index)}
    
    draw_survival_curve(feature, surv, ax=ax1, colors=color_lookup, **args)
    ax1.legend().set_visible(False)
    if title:
        ax1.set_title(title)
    
    fig.tight_layout()

