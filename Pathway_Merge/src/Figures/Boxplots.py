'''
Created on Apr 24, 2013

@author: agross
'''
import numpy as np
import pandas as pd
import matplotlib.pylab as plt

import Stats.Scipy as Stats
from Figures.Helpers import latex_float, init_ax
from Processing.Helpers import match_series
colors = plt.rcParams['axes.color_cycle']*10
        
def _violin_plot(ax, data,pos=[], bp=False):
    '''
    http://pyinsci.blogspot.com/2009/09/violin-plot-with-matplotlib.html
    
    Create violin plots on an axis.  Internal to module as it does not 
    use Pandas data-structures.  This is split off due to it's being a 
    reuse of the code from the blog-post linked above, and I wanted to keep
    the original code untouched. 
    '''
    from scipy.stats import gaussian_kde
    from numpy import arange
    
    #dist = max(pos)-min(pos)
    dist = len(pos)
    w = min(0.25*max(dist,1.0),0.5)
    for p,d in enumerate(data):
        try:
            k = gaussian_kde(d) #calculates the kernel density
            m = k.dataset.min() #lower bound of violin
            M = k.dataset.max() #upper bound of violin
            x = arange(m,M,(M-m)/100.) # support for violin
            v = k.evaluate(x) #violin profile (density curve)
            v = v/v.max()*w #scaling the violin to the available space
            ax.fill_betweenx(x,p,v+p,facecolor='y',alpha=0.1)
            ax.fill_betweenx(x,p,-v+p,facecolor='y',alpha=0.1)
        except:
            pass
    if bp:
        boxPlot = ax.boxplot(data,notch=1,positions=range(len(pos)),vert=1, 
                             widths=.25)
        return boxPlot
    
def box_plot_pandas(bin_vec, real_vec, ax=None):
    '''
    Wrapper around matplotlib's boxplot function.
    
    Inputs
        bin_vec: Series of labels
        real_vec: Series of measurements to be grouped according to bin_vec
    '''
    _, ax = init_ax(ax)
    bin_vec, real_vec = match_series(bin_vec, real_vec)
    categories = bin_vec.value_counts().index
    data = [real_vec[bin_vec==num] for num in categories]
    bp = ax.boxplot(data, positions=range(len(categories)), widths=.3,
                    patch_artist=True);
    if real_vec.name:
        ax.set_ylabel(real_vec.name)
    if bin_vec.name:
        ax.set_xlabel(bin_vec.name)
    [p.set_visible(False) for p in bp['fliers']]
    [p.set_visible(False) for p in bp['caps']]
    [p.set_visible(False) for p in bp['whiskers']]
    for p in bp['medians']:
        p.set_color(colors[0])
        p.set_lw(3)
        p.set_alpha(.8)
    for p in bp['boxes']:
        p.set_color('grey')
        p.set_lw(3)
        p.set_alpha(.7)
        
def violin_plot_pandas(bin_vec, real_vec, ann='p', order=None, ax=None, 
                       filename=None):
    '''
    http://pyinsci.blogspot.com/2009/09/violin-plot-with-matplotlib.html
    Wrapper around matplotlib's boxplot function to add violin profile.
    
    Inputs
        bin_vec: Series of labels
        real_vec: Series of measurements to be grouped according to bin_vec
    '''   
    fig, ax = init_ax(ax)
    ax.set_ylabel(real_vec.name)
    ax.set_xlabel(bin_vec.name)
    bin_vec, real_vec = match_series(bin_vec, real_vec)
    try:
        if order is None:
            categories = bin_vec.value_counts().index
        else:
            categories = order
        _violin_plot(ax, [real_vec[bin_vec==num] for num in categories], 
                     pos=categories, bp=True)
        ax.set_xticklabels([str(c) +'\n(n=%i)'%sum(bin_vec==c) 
                            for c in categories])
    except:
        box_plot_pandas(bin_vec, real_vec, ax=ax)
        
    if type(bin_vec.name) == str:
        ax.set_title(str(bin_vec.name) +' x '+ str(real_vec.name))
        
    p_value = Stats.kruskal_pandas(bin_vec, real_vec)['p']
    if ann == 'p_fancy':
        ax.annotate('$p = {}$'.format(latex_float(p_value)), (.95, -.02),
                    xycoords='axes fraction', ha='right', va='bottom', size=14)
    if ann == 'p':
        ax.annotate('p = {0:.1e}'.format(p_value), (.95, .02),
                    xycoords='axes fraction', ha='right', va='bottom', size=12)
    elif ann != None:
        ax.annotate(ann, (.95, .02), xycoords='axes fraction', ha='right', 
                    va='bottom', size=12)
    if filename is not None:
        fig.savefig(filename)
    return

def violin_plot_series(s, **kw_args):
    '''
    Wrapper for drawing a violin plot on a series with a multi-index.
    The second level of the index is used as the binning variable. 
    '''
    assert s.index.levshape[1] > 1
    violin_plot_pandas(pd.Series(s.index.get_level_values(1), s.index), s,
                        **kw_args)


def paired_boxplot_o(boxes):
    '''
    Wrapper around plt.boxplot to draw paired boxplots
    for a set of boxes. 
    
    Input is the same as plt.boxplot:
        Array or a sequence of vectors.
    '''
    fig = plt.figure(figsize=(len(boxes)/2.5,4))
    ax1 = fig.add_subplot(111)
    plt.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)
    bp = ax1.boxplot(boxes, notch=0, positions=np.arange(len(boxes)) + 
                     1.5*(np.arange(len(boxes)) / 2), patch_artist=True);
    [p.set_color(colors[0]) for p in bp['boxes'][::2]]
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

def paired_boxplot(boxes, ax1=None):
    if not ax1:
        fig = plt.figure(figsize=(len(boxes)/2.5,4))
        ax1 = fig.add_subplot(111)
    plt.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)
    bp = ax1.boxplot(boxes, notch=0, positions=np.arange(len(boxes)) + 
                     1.5*(np.arange(len(boxes)) / 2), patch_artist=True);
    [p.set_color(colors[0]) for p in bp['boxes'][::2]]
    [p.set_color(colors[1]) for p in bp['boxes'][1::2]]
    [p.set_color('black') for p in bp['whiskers']]
    [p.set_color('black') for p in bp['fliers']]
    [p.set_alpha(.4) for p in bp['fliers']]
    [p.set_alpha(.8) for p in bp['boxes']];
    [p.set_edgecolor('black') for p in bp['boxes']];
    ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
                  alpha=0.5)
    
    # Hide these grid behind plot objects
    ax1.set_axisbelow(True)
    ax1.set_ylabel('$Log_{2}$ RNA Expression')
    ax1.set_xticks(3.5*np.arange(len(boxes)/2) + .5);
    return ax1, bp

def paired_boxplot_tumor_normal(df, sig=True, cutoffs=[.01, .00001], 
                                order=None, ax=None):
    '''
    Draws a paired boxplot given a DataFrame with both tumor and normal
    samples on the index. '01' and '11' are hard-coded as the ids for
    tumor/normal. 
    '''
    n = df.groupby(level=0).size()==2
    df = df.ix[n[n].index]
    if order is None:
        o = df.xs('11', level=1).median().order().index
        df = df[o[::-1]]
    else:
        df = df[order]
    l1 = list(df.xs('01', level=1).as_matrix().T)
    l2 = list(df.xs('11', level=1).as_matrix().T)
    boxes = [x for t in zip(l1, l2) for x in t]
    ax1, bp = paired_boxplot(boxes, ax)
    
    test = lambda v: Stats.ttest_rel(v.unstack()['01'], v.unstack()['11'])
    res = df.apply(test).T
    p = res.p
    
    if sig:
        pts = [(i*3.5 +.5,18) for i,n in enumerate(p) if n < cutoffs[1]]
        if len(pts) > 0:
            s1 = ax1.scatter(*zip(*pts), marker='$**$', label='$p<10^{-5}$', s=200)
        else:
            s1 = None
        pts = [(i*3.5 +.5,18) for i,n in enumerate(p) 
                              if (n < cutoffs[0]) and (n > cutoffs[1])]
        if len(pts) > 0:
            s2 = ax1.scatter(*zip(*pts), marker='$*$', label='$p<10^{-2}$', s=30)
        else:
            s2 = None
        ax1.legend(bp['boxes'][:2] + [s2,s1], 
                   ('Tumor','Normal', '$p<10^{-2}$', '$p<10^{-5}$'), 
                   loc='best', scatterpoints=1);
    else:
        ax1.legend(bp['boxes'][:2], ('Tumor','Normal'), loc='best');
    ax1.set_xticklabels(df.columns);

    
def boxplot_pannel(hit_vec, response_df):
    '''
    Draws a series of paired boxplots with the rows of the response_df
    split according to hit_vec.  
    '''
    b = response_df.copy()
    b.columns = pd.MultiIndex.from_arrays([b.columns, hit_vec.ix[b.columns]])
    b = b.T
    v1, v2 = hit_vec.unique()
    test = lambda v: Stats.anova(v.reset_index(level=1)[v.index.names[1]], 
                                 v.reset_index(level=1)[v.name])
    res = b.apply(test).T
    p = res.p.order()
    b = b.ix[:,p.index]
    
    l1 = list(b.xs(v1, level=1).as_matrix().T)
    l2 = list(b.xs(v2, level=1).as_matrix().T)

    boxes = [x for t in zip(l1, l2) for x in t]
    ax1, bp = paired_boxplot(boxes)
    
    y_lim = (response_df.T.quantile(.9).max()) * 1.2
    pts = [(i*3.5 +.5, y_lim) for i,n in enumerate(p) if n < .00001]
    if len(pts) > 0:
        s1 = ax1.scatter(*zip(*pts), marker='$**$', label='$p<10^{-5}$', s=200)
    else:
        s1 = None
    pts = [(i*3.5 +.5, y_lim) for i,n in enumerate(p) 
                              if (n < .01) and (n > .00001)]
    if len(pts) > 0:
        s2 = ax1.scatter(*zip(*pts), marker='$*$', label='$p<10^{-2}$', s=30)
    else:
        s2 = None
    ax1.set_xticklabels(b.columns);
    ax1.legend(bp['boxes'][:2] + [s2,s1], 
               (v1, v2, '$p<10^{-2}$', '$p<10^{-5}$'), 
               loc='best', scatterpoints=1);