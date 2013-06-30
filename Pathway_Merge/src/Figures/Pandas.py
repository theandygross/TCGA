'''
Created on Jun 12, 2013

@author: agross
'''
import Processing.Tests as Tests
from Processing.Helpers import match_series, split_a_by_b
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy as sp

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
        
def violin_plot(ax,data,pos=[], bp=False):
    '''
    http://pyinsci.blogspot.com/2009/09/violin-plot-with-matplotlib.html
    create violin plots on an axis
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
            do_nothing = True
    if bp:
        boxPlot = ax.boxplot(data,notch=1,positions=range(len(pos)),vert=1, widths=.25)
        return boxPlot
    
def box_plot_pandas(hitVec, expVec, ax='None'):
    '''
    http://pyinsci.blogspot.com/2009/09/violin-plot-with-matplotlib.html
    Wrapper around matplotlib's boxplot function for KW eQTLs
    '''
    if ax=='None':
        fig, ax = plt.subplots(1,1)
    ax.boxplot([expVec[hitVec==num] for num in set(hitVec)], 
           positions=list(set(hitVec)));
    ax.set_ylabel('Sub-Cohort Gene Expression')
    ax.set_xlabel('Number of Mutations')
    if type(hitVec.name) == str:
        ax.set_title(hitVec.name +' x '+ expVec.name)
        
def violin_plot_pandas(bin_vec, real_vec, ann='p', ax=None, filename=None):
    '''
    http://pyinsci.blogspot.com/2009/09/violin-plot-with-matplotlib.html
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
    
def fischer_bar_chart(bin_vec, response_vec, ax=None, filename=None):
    if ax is None:
        fig, ax = plt.subplots(1,1)
    t = pd.crosstab(bin_vec, response_vec)
    t.plot(kind='bar', ax=ax)
    if filename is not None:
        fig.savefig(filename)
    return fig 

def draw_pathway_count_bar(p, cancer, file_name='tmp.svg', colors='red'):
    fig, ax = plt.subplots(1,1, figsize=(3+len(cancer.gene_sets[p])/30.,2))
    if 'deceased' in cancer.clinical:
        patients = cancer.clinical[['age', 'deceased']].dropna().index
    else:
        patients = cancer.clinical.index
    m = cancer.hit_matrix.ix[list(cancer.gene_sets[p]), patients] > 0
    m = np.sort(m[m.sum(1) > 0].sum(1))
    #colors = ['orange']*(len(m))
    m.plot(kind='bar', ax=ax, alpha=.7, color=colors);
    ax.set_yticks(range(m.max()+1))
    ax.set_yticklabels(range(m.max()+1), size=14)
    ax.set_xticklabels(m.index, ha='center', va='bottom', position=(0,.1), size=20)
    ax.set_ylabel('# Patients', size=16)
    
    fig.tight_layout()
    fig.savefig(file_name)
    
def draw_pathway_eig_bar(U, file_name='tmp.svg'):
    if sp.rank(U)  == 2:
        U = U[0]
    fig, ax = plt.subplots(1,1, figsize=(3+len(U)/15.,2.5))
    np.sort(U).plot(kind='bar', ax=ax, alpha=.4)
    ax.set_ylabel('Loading')
    fig.tight_layout()
    fig.savefig(file_name)
    
def draw_pathway_age_scatter(p, cancer, file_name='tmp.svg'):
    fig, ax = plt.subplots(1,1, figsize=(6,4))
    ax.scatter(*match_series(cancer.pc.ix[p], cancer.clinical.age), alpha=.5, s=75)
    ax.set_xlabel('Principal Component Loading')
    ax.set_ylabel('Age')
    fig.savefig(file_name) 
    
def histo_compare(hit_vec, response_vec, ax=None):
    '''
    Split response_vec by hit_vec and compared histograms.  
    Also plots the kde of the whole response_vec.
    '''
    if ax is None:
        fig, ax = plt.subplots(1,1)
    else:
        fig = fig.gcf()
    kde1 = sp.stats.gaussian_kde(response_vec)
    x_eval = np.linspace(min(response_vec), max(response_vec), num=200)
    ax.plot(x_eval, kde1(x_eval), 'k-')
    miss, hit = split_a_by_b(response_vec, hit_vec)
    ax.hist(miss, bins=20, normed=True, alpha=.2, label='WT');
    ax.hist(hit, bins=10, normed=True, alpha=.5, label='Mut');
    ax.legend()
    return fig

def mut_module_raster(cluster_num, mut, ax=None):
    assert hasattr(mut, 'clustered')
    if ax is None:
        fig, ax = plt.subplots(1,1, figsize=(10,8))
    else:
        fig = fig.gcf()
    p = mut.assignments[mut.assignments == cluster_num].index
    s = mut.meta_matrix.ix[p].sum()
    x_order = np.sort(-s[s > 0]).index
    y_order = np.sort(mut.meta_matrix.ix[p].sum(1)).index
    plt.imshow(mut.meta_matrix.ix[y_order, x_order], aspect=8, 
               interpolation='Nearest', cmap=plt.cm.bone_r) #@UndefinedVariable
    ax.set_yticks(range(len(p)))
    ax.set_yticklabels(p);
    ax.hlines(np.arange(len(p)) + .5, 0, len(x_order), colors='grey', alpha=.3, lw=4)
    ax.set_xlabel('Patients');
    ax.axvline(x=(mut.clustered.ix[cluster_num] == 1).sum(), ymin=0, 
               ymax=len(y_order), color='r', alpha=.7, lw=5)
    ax.set_xbound(0, len(x_order))
    fig.tight_layout()
    return fig

def mut_box(x, y, width=.7, height=.9, small=False):
    '''
    Takes a coordinate and returns a mutation box rectangle.
    Small makes it a little one. 
    '''
    from matplotlib.patches import FancyBboxPatch
    
    if small:
        height = height / 2
        y = y - height / 2
    width = width - .2
    height = height - .2
    rect = FancyBboxPatch([x - width / 2, y - height / 2], width, height,
                          facecolor='green', edgecolor='black', alpha=.8,
                          lw=2, mutation_scale=.33)
    return rect
    
def null_box(x, y, width=.7, height=.9):
    '''
    Grey null box for MeMO plots. 
    '''
    from matplotlib.patches import FancyBboxPatch
    
    width = width - .2
    height = height - .2
    rect = FancyBboxPatch([x - width / 2, y - height / 2], width, height,
                     facecolor='grey', edgecolor='black', alpha=.5,
                     mutation_scale=.33)
    return rect

def memo_plot(df, ax=None):
    '''
    MeMO style plot from a DataFrame.
    '''
    if ax == None:
        ax = plt.gca()
    ax.patch.set_facecolor('white')
    ax.set_aspect(3, 'box')
    ax.set_xticks([])
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df.index[::-1])
    
    W = df.as_matrix().T[:, ::-1]
    for (x,y),w in np.ndenumerate(W):
        if w == 0: 
            rect = null_box(x, y)
        else:     
            rect = null_box(x, y)
            ax.add_patch(rect)
            rect = mut_box(x, y) 
        ax.add_patch(rect)
    df
    ax.set_xbound(-.5, len(W)-.5)
    ax.set_ybound(-.5, len(W[0])-.5)
    
def pathway_plot(df, ax=None):
    df = df.ix[df.sum(1) > 0, df.sum() > 0]
    df = df.ix[df.sum(1).order(ascending=False).index]
    o = np.sort(df.apply(lambda s: ''.join(map(str, s)))).index[::-1]
    df = df[o]
    
    if df.shape[0] > 20:
        rest = pd.Series(df.ix[10:].sum().clip_upper(1.), name='rest')
        df = df.ix[:10]
        df = df.append(rest)
    if ax is None:
        fig, ax = plt.subplots(figsize=(df.shape[1]*.2,df.shape[0]*.5))
    else:
        fig = ax.get_figure()
    memo_plot(df, ax=ax)
    ax.bar(np.arange(len(df.columns)) - .3, df.sum() / df.sum().max(), bottom=-1.5, 
           width=.6, alpha=.5)
    counts = df.sum(1)[::-1]
    width = df.shape[1]
    ax.barh(np.arange(len(counts)) - .3, (counts / counts.max())*width*.25, left=width - .2, 
            height=.6, alpha=.5)
    ax.set_frame_on(False)
    ax.tick_params(right='off')
    fig.tight_layout()
    
def draw_pathway_overlaps(mat, bars, filename=None):  
    fig, ax = plt.subplots(figsize=(25,5))   
    memo_plot(mat, ax=ax)
    ax.bar(np.arange(len(mat.columns)) - .3, bars, bottom = -2, width=.6, alpha=.5)
    ax.hlines(-1, -.5, len(mat.columns))
    ax.annotate('Days to Death', (-.5,-1.5), ha='right', va='center', size=15)
    fig.tight_layout()
    if filename is not None:
        fig.savefig(filename)
    return fig



def fancy_raster(df, cluster=False, cmap=plt.cm.get_cmap('Spectral'), norm=None):
    if cluster:
        d = sp.spatial.distance.pdist(df)
        D = sp.spatial.distance.squareform(d)
        Y = sp.cluster.hierarchy.linkage(D)
        Z = sp.cluster.hierarchy.dendrogram(Y, no_plot=True)
        order = Z['leaves']
        df = df.ix[order, order]
        
    fig, ax = plt.subplots(1,1, figsize=(12,8))
    img = ax.imshow(df, interpolation='Nearest', cmap=cmap, norm=norm)
    ax.set_yticks(range(len(df.index)))
    ax.set_yticklabels(df.index)
    ax.set_xticks(np.arange(len(df.columns)))
    ax.set_xticklabels(df.columns, rotation=360-90, ha='center');
    ax.hlines(np.arange(len(df.index)-1)+.5, -.5, len(df.columns)-.5, 
              color='white', lw=6)
    ax.vlines(np.arange(len(df.columns)-1)+.5, -.5, len(df.index)-.5, 
              color='white', lw=6)
    
    if cluster:
        icoord = np.array(Z['icoord']) - np.array(Z['icoord']).min()
        icoord = icoord * ((len(Z['leaves']) - 1) / icoord.max())
    
        dcoord = -1*np.array(Z['dcoord']) - .7 
        for i,z,c in zip(icoord, dcoord, Z['color_list']):
            ax.plot(i,z,color=c, lw=2, alpha=.8)
            
        ax.tick_params(axis='x', top='off')
        ax.set_frame_on(False)
    return img

def count_plot(vec, name=None, ax=None):
    if ax is None:
        ax = plt.gca()
    vec.value_counts().sort_index().plot(kind='bar', ax=ax)
    ax.set_ylabel('# of Patients')
    ax.set_xlabel(name if name is not None else vec.name)

def venn_pandas(a,b):
    from matplotlib_venn import venn2
    
    colors = plt.rcParams['axes.color_cycle']
    gc = pd.concat([a,b], axis=1).dropna().astype(int).astype(str).apply(lambda s: ''.join(s), axis=1)
    v = venn2(gc.value_counts().sort_index()[1:], set_labels=[b.name, a.name], normalize_to=1.)
    v.patches[0].set_facecolor(colors[0])
    v.patches[1].set_facecolor(colors[2])
    v.patches[2].set_facecolor(colors[4])
    v.patches[0].set_alpha(.7)
    v.patches[1].set_alpha(.7)
    v.patches[2].set_alpha(.7)