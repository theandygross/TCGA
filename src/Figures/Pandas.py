'''
Created on Jun 12, 2013

@author: agross
'''
import Stats.Scipy as Tests
from Processing.Helpers import match_series, split_a_by_b
from Figures.Helpers import init_ax, latex_float
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy as sp

colors = plt.rcParams['axes.color_cycle']

def series_scatter(s1, s2, ax=None, ann='p', filename=None, **plot_args):
    fig, ax = init_ax(ax, figsize=(6, 4))
    if 's' not in plot_args:
        plot_args['s'] = 75
    if 'alpha' not in plot_args:
        plot_args['alpha'] = .5
    ax.scatter(*match_series(s1, s2), **plot_args)
    ax.set_xlabel(s1.name)
    ax.set_ylabel(s2.name)
    if ann == 'p':
        ax.annotate('p = {0:.1e}'.format(Tests.spearman_pandas(s1, s2)['p']), (.95, -.02),
                    xycoords='axes fraction', ha='right', va='bottom', size=14)
    if ann == 'fancy_p':
        ax.annotate('$p = {}$'.format(latex_float(Tests.spearman_pandas(s1, s2)['p'])), (.95, -.02),
                    xycoords='axes fraction', ha='right', va='bottom', size=14)
    if filename is not None:
        fig.savefig(filename)
    
def fischer_bar_chart(bin_vec, response_vec, ax=None, filename=None):
    fig, ax = init_ax(ax)
    t = pd.crosstab(bin_vec, response_vec)
    t.plot(kind='bar', ax=ax)
    if filename is not None:
        fig.savefig(filename)
    return fig     
    
def histo_compare(hit_vec, response_vec, ax=None):
    '''
    Split response_vec by hit_vec and compared histograms.  
    Also plots the kde of the whole response_vec.
    '''
    fig, ax = init_ax(ax)
    kde1 = sp.stats.gaussian_kde(response_vec)
    x_eval = np.linspace(min(response_vec), max(response_vec), num=200)
    ax.plot(x_eval, kde1(x_eval), 'k-')
    miss, hit = split_a_by_b(response_vec, hit_vec)
    ax.hist(miss, bins=20, normed=True, alpha=.2, label='WT');
    ax.hist(hit, bins=10, normed=True, alpha=.5, label='Mut');
    ax.legend()
    return fig

def fancy_raster(df, cluster=False, cmap=plt.cm.get_cmap('Spectral'),
                 norm=None, ax=None):
    if cluster:
        d = sp.spatial.distance.pdist(df)
        D = sp.spatial.distance.squareform(d)
        Y = sp.cluster.hierarchy.linkage(D)
        Z = sp.cluster.hierarchy.dendrogram(Y, no_plot=True)
        order = Z['leaves']
        df = df.ix[order, order]
        
    _, ax = init_ax(ax, figsize=(12, 8))
    img = ax.imshow(df, interpolation='Nearest', cmap=cmap, norm=norm)
    ax.set_yticks(range(len(df.index)))
    ax.set_yticklabels(df.index)
    ax.set_xticks(np.arange(len(df.columns)))
    ax.set_xticklabels(df.columns, rotation=360 - 90, ha='center');
    ax.hlines(np.arange(len(df.index) - 1) + .5, -.5, len(df.columns) - .5,
              color='white', lw=6)
    ax.vlines(np.arange(len(df.columns) - 1) + .5, -.5, len(df.index) - .5,
              color='white', lw=6)
    
    if cluster:
        icoord = np.array(Z['icoord']) - np.array(Z['icoord']).min()
        icoord = icoord * ((len(Z['leaves']) - 1) / icoord.max())
    
        dcoord = -1 * np.array(Z['dcoord']) - .7 
        for i, z, c in zip(icoord, dcoord, Z['color_list']):
            ax.plot(i, z, color=c, lw=2, alpha=.8)
            
        ax.tick_params(axis='x', top='off')
        ax.set_frame_on(False)
    return img

def count_plot(vec, name=None, ax=None):
    _, ax = init_ax(ax)
    vec.value_counts().sort_index().plot(kind='bar', ax=ax)
    ax.set_ylabel('# of Patients')
    ax.set_xlabel(name if name is not None else vec.name)

def venn_pandas_o(a, b):
    from matplotlib_venn import venn2
    
    colors = plt.rcParams['axes.color_cycle']
    gc = pd.concat([a, b], axis=1).dropna().astype(int).astype(str).apply(lambda s: ''.join(s), axis=1)
    v = venn2(gc.value_counts().sort_index()[1:], set_labels=[b.name, a.name], normalize_to=1.)
    v.patches[0].set_facecolor(colors[0])
    v.patches[1].set_facecolor(colors[2])
    v.patches[2].set_facecolor(colors[4])
    v.patches[0].set_alpha(.7)
    v.patches[1].set_alpha(.7)
    v.patches[2].set_alpha(.7)
    
def venn_pandas(a, b, colors=None, alpha=.7):
    from matplotlib_venn import venn2
    
    if colors is None:
        colors = np.array(plt.rcParams['axes.color_cycle'])[[0, 2, 4]]
    gc = pd.concat([a, b], axis=1).dropna().astype(int).astype(str).apply(lambda s: ''.join(s), axis=1)
    v = venn2(gc.value_counts().sort_index()[1:], set_labels=[b.name, a.name], normalize_to=1.0)
    v.patches[0].set_facecolor(colors[0])
    v.patches[1].set_facecolor(colors[1])
    v.patches[2].set_facecolor(colors[2])
    for p in v.patches:
        p.set_alpha(alpha)
        p.set_lw(2)
    for l in v.subset_labels:
        l.set_fontsize(12)
    return v
