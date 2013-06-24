'''
Created on Jun 12, 2013

@author: agross
'''
import Processing.Tests as Tests
import matplotlib.pylab as plt
from Processing.Helpers import *
from Reports.Figures import violin_plot, box_plot_pandas
import pandas as pd

colors = plt.rcParams['axes.color_cycle']

def latex_float(f):
    '''http://stackoverflow.com/questions/13490292/format-number-using-latex-notation-in-python'''
    float_str = "{0:.2g}".format(f)
    if "e" in float_str:
        base, exponent = float_str.split("e")
        return r"{0} \times 10^{{{1}}}".format(base, int(exponent))
    else:
        return float_str

def series_scatter(s1, s2, ax=None, ann='p', filename=None, **plot_args):
    if ax is None:
        fig, ax = plt.subplots(1,1, figsize=(6,4))
    if 's' not in plot_args:
        plot_args['s'] = 75
    if 'alpha' not in plot_args:
        plot_args['alpha'] = .5
    ax.scatter(*match_series(s1, s2), **plot_args)
    ax.set_xlabel(s1.name)
    ax.set_ylabel(s2.name)
    if ann == 'p':
        ax.annotate('p = {0:.1e}'.format(Tests.spearman_pandas(s1, s2)['p']), (.95, -.02),
                    xycoords='axes fraction', ha='right',va='bottom', size=14)
    if ann == 'fancy_p':
        ax.annotate('$p = {}$'.format(latex_float(Tests.spearman_pandas(s1, s2)['p'])), (.95, -.02),
                    xycoords='axes fraction', ha='right',va='bottom', size=14)
    if filename is not None:
        fig.savefig(filename)
        
def violin_plot_pandas(bin_vec, real_vec, ann='p', ax=None, filename=None):
    '''
    Wrapper around matplotlib's boxplot function to add violin profile.
    '''   
    bin_vec, real_vec = match_series(bin_vec, real_vec)
    if ax is None:
        fig, ax = plt.subplots(1,1)
    else:
        fig = plt.gcf()
    try:
        categories = bin_vec.value_counts().index
        violin_plot(ax, [real_vec[bin_vec==num] for num in categories], 
                    pos=categories, bp=True)
        ax.set_xticklabels([str(c) +'\n(n=%i)'%sum(bin_vec==c) 
                            for c in categories])
    except:
        box_plot_pandas(bin_vec, real_vec, ax=ax)
    ax.set_ylabel(real_vec.name)
    ax.set_xlabel(bin_vec.name)
    if type(bin_vec.name) == str:
        ax.set_title(str(bin_vec.name) +' x '+ str(real_vec.name))
    if ann == 'p_fancy':
        ax.annotate('$p = {}$'.format(latex_float(Tests.kruskal_pandas(bin_vec, real_vec)['p'])), 
                    (.95, -.02),
                    xycoords='axes fraction', ha='right',va='bottom', size=14)
    if ann == 'p':
        ax.annotate('p = {0:.1e}'.format(Tests.kruskal_pandas(bin_vec, real_vec).ix['p']), 
                    (.95, .02),
                    xycoords='axes fraction', ha='right',va='bottom', size=12)
    elif ann != None:
        ax.annotate(ann, (.95, .02),
                    xycoords='axes fraction', ha='right',va='bottom', size=12)
    if filename is not None:
        fig.savefig(filename)
    return

def violin_plot_series(s, **kw_args):
    '''
    Wrapper for drawing a violin plot on a series with a multi-index.
    The second level of the index is used as the binning variable. 
    '''
    assert s.index.levshape[1] == 2
    violin_plot_pandas(pd.Series(s.index.get_level_values(1), s.index), s, **kw_args)