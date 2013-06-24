'''
Created on Apr 24, 2013

@author: agross
'''
from scipy import stats
import numpy as np
import pandas as pd
import matplotlib.pylab as plt
from Processing.Tests import anova

from Processing.Helpers import match_series

def ttest_rel(a,b):
    a,b = match_series(a,b)
    z, p = stats.ttest_rel(a, b)
    return pd.Series({'t': z, 'p': p})

def paired_boxplot(boxes):
    fig = plt.figure(figsize=(len(boxes)/2.5,4))
    ax1 = fig.add_subplot(111)
    plt.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)
    bp = ax1.boxplot(boxes, notch=0, positions=np.arange(len(boxes)) + 
                     1.5*(np.arange(len(boxes)) / 2), patch_artist=True);
    [p.set_color('r') for p in bp['boxes'][::2]]
    [p.set_color('black') for p in bp['whiskers']]
    [p.set_color('black') for p in bp['fliers']]
    [p.set_alpha(.4) for p in bp['fliers']]
    [p.set_alpha(.6) for p in bp['boxes']];
    [p.set_edgecolor('black') for p in bp['boxes']];
    ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
                  alpha=0.5)
    
    # Hide these grid behind plot objects
    ax1.set_axisbelow(True)
    ax1.set_ylabel('$Log_{2}$ RNA Expression')
    ax1.set_xticks(3.5*np.arange(len(boxes)/2) + .5);
    return ax1, bp

def paired_boxplot2(b):
    n = b.groupby(level=0).size()==2
    b = b.ix[n[n].index]
    o = b.xs('11', level=1).median().order().index
    b = b[o[::-1]]
    l1, l2 = list(b.xs('01', level=1).as_matrix().T), list(b.xs('11', level=1).as_matrix().T)
    boxes = [x for t in zip(l1, l2) for x in t]
    ax1, bp = paired_boxplot(boxes)
    
    test = lambda v: ttest_rel(v.unstack()['01'], v.unstack()['11'])
    res = b.apply(test).T
    p = res.p
    
    pts = [(i*3.5 +.5,18) for i,n in enumerate(p) if n < .00001]
    if len(pts) > 0:
        s1 = ax1.scatter(*zip(*pts), marker='$**$', label='$p<10^{-5}$', s=200)
    else:
        s1 = None
    pts = [(i*3.5 +.5,18) for i,n in enumerate(p) if (n < .01) and (n > .00001)]
    if len(pts) > 0:
        s2 = ax1.scatter(*zip(*pts), marker='$*$', label='$p<10^{-2}$', s=30)
    else:
        s2 = None
    ax1.set_xticklabels(b.columns);
    ax1.legend(bp['boxes'][:2] + [s2,s1], ('Tumor','Normal', '$p<10^{-2}$', '$p<10^{-5}$'), loc='best',
               scatterpoints=1);
    
def boxplot_pannel(hit_vec, response_vec):
    
    b = response_vec.copy()
    b.columns = pd.MultiIndex.from_arrays([b.columns, hit_vec.ix[b.columns]])
    b = b.T
    
    v1, v2 = hit_vec.unique()
    test = lambda v: anova(v.reset_index(level=1)[v.index.names[1]], v.reset_index(level=1)[v.name])
    res = b.apply(test).T
    p = res.p.order()
    b = b.ix[:,p.index]
    
    l1, l2 = list(b.xs(v1, level=1).as_matrix().T), list(b.xs(v2, level=1).as_matrix().T)
    boxes = [x for t in zip(l1, l2) for x in t]
    ax1, bp = paired_boxplot(boxes)
        
    pts = [(i*3.5 +.5,18) for i,n in enumerate(p) if n < .00001]
    if len(pts) > 0:
        s1 = ax1.scatter(*zip(*pts), marker='$**$', label='$p<10^{-5}$', s=200)
    else:
        s1 = None
    pts = [(i*3.5 +.5,18) for i,n in enumerate(p) if (n < .01) and (n > .00001)]
    if len(pts) > 0:
        s2 = ax1.scatter(*zip(*pts), marker='$*$', label='$p<10^{-2}$', s=30)
    else:
        s2 = None
    ax1.set_xticklabels(b.columns);
    ax1.legend(bp['boxes'][:2] + [s2,s1], (v1, v2, '$p<10^{-2}$', '$p<10^{-5}$'), loc='best',
                   scatterpoints=1);