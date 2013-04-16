'''
Created on Oct 15, 2012

@author: agross
'''

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from pylab import cm
from numpy import nanmax, sort, linspace, arange, rank, ndenumerate, array
import numpy as np
from scipy.stats import gaussian_kde

from Processing.Helpers import match_series, split_a_by_b
from Processing.Helpers import get_vec_type ,to_quants
from Processing.Tests import get_cox_ph_ms, pearson_p, anova

import rpy2.robjects as robjects 
from rpy2.robjects import r
from pandas import crosstab, Series, DataFrame
import pandas.rpy.common as com 
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as dist

import base64
from IPython.display import Image

survival = robjects.packages.importr('survival')
base = robjects.packages.importr('base')


def violin_plot(ax,data,pos=[], bp=False):
    '''
    create violin plots on an axis
    '''
    from scipy.stats import gaussian_kde
    from numpy import arange
    
    #dist = max(pos)-min(pos)
    dist = len(pos)
    w = min(0.25*max(dist,1.0),0.5)
    for p,d in enumerate(data):
        k = gaussian_kde(d) #calculates the kernel density
        m = k.dataset.min() #lower bound of violin
        M = k.dataset.max() #upper bound of violin
        x = arange(m,M,(M-m)/100.) # support for violin
        v = k.evaluate(x) #violin profile (density curve)
        v = v/v.max()*w #scaling the violin to the available space
        ax.fill_betweenx(x,p,v+p,facecolor='y',alpha=0.1)
        ax.fill_betweenx(x,p,-v+p,facecolor='y',alpha=0.1)
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
    if ann == 'p':
        ax.annotate('p = {0:.2e}'.format(anova(bin_vec, real_vec)), (.95, .02),
                    xycoords='axes fraction', ha='right',va='bottom', size=12)
    elif ann != None:
        ax.annotate(ann, (.95, .02),
                    xycoords='axes fraction', ha='right',va='bottom', size=12)
    if filename is not None:
        fig.savefig(filename)
    return fig
    
def fischer_bar_chart(bin_vec, response_vec, ax=None, filename=None):
    if ax is None:
        fig, ax = plt.subplots(1,1)
    t = crosstab(bin_vec, response_vec)
    t.plot(kind='bar', ax=ax)
    if filename is not None:
        fig.savefig(filename)
    return fig    
    
        
def draw_survival_curves_KM(clinical, hit_vec, time_var='days', event_var='deceased',
                            filename='tmp.png', show=False, ax=None, title=True, 
                            labels=['No Mutation', 'Mutation'], colors=['blue','red'], 
                            ann=None):
    if not hasattr(hit_vec, 'name'):
        hit_vec.name = 'pathway'
    df = clinical[[time_var, event_var]].join(hit_vec)
    df = df.dropna()
    df_r = com.convert_to_r_dataframe(df) #@UndefinedVariable
    fmla = robjects.Formula('Surv(' + time_var + ', ' + event_var + ') ~ ' + 
                            hit_vec.name)
    fit = survival.survfit(fmla, df_r)
    ls = r('2:' + str(len(set(hit_vec))+1)) #R line styles
    if filename.endswith('.png'):
        r.png(filename=filename, width=300, height=250, res=100, pointsize=8)
    else:
        r.pdf(filename, width=4.5, height=3.75)
    s = survival.survdiff(fmla, df_r)
    p = str(s).split('\n\n')[-1].strip().split(', ')[-1]
    ls = r.c(*colors)
    r.plot(survival.survfit(fmla, df_r), lty=1, col=ls, lwd=3, cex=1.5, 
                            xlab='Days to Event', ylab= 'Survival');
    if title:
        r('title')(hit_vec.name)
    r('legend')(nanmax(df[time_var]) * .5,.9, labels, lty=1, col=ls, lwd=3)
    if ann=='p':
        r.text(0, labels='logrank ' + p, pos=4)
    elif ann != None:
        r.text(0, labels=ann, pos=4)
    r('dev.off()')
        
    
def draw_survival_curves(feature, surv, assignment=None, filename='tmp.png', show=False, 
                               title=True, labels=None, 
                               colors=['blue','red'], ann=None, show_legend=True, q=.25, std=None):
    if assignment is None:
        num_panels = 1
        assignment = np.ones_like(feature)
        name = lambda v: str(feature.name) if feature.name != None else ''
    else:
        num_panels = len(assignment.unique())
        name = lambda v: str(assignment.name) + ' = ' + str(v)
    if (labels is None) and (feature.nunique() < 10):
        labels = r.sort(r.c(*feature.unique()))  #R sorts bad
        colors = ['blue','green','red','cyan','magenta','yellow','black']
    if feature.dtype == 'bool':
        feature = feature.map({True: 'True', False: 'False'})
        
    r.png(filename=filename, width=200*(num_panels+1), height=300, res=75)
        
    fmla = robjects.Formula('Surv(days, event) ~ feature')
    r.par(mfrow=r.c(1, num_panels))
    r.par(mar=r.c(4,5,4,1))
    r.par(xpd=True)
    
    if (get_vec_type(feature) == 'real') and (len(feature.unique()) > 10):
        colors=['blue','orange','red']
        if q == .5:
            labels=['Bottom 50%', 'Top 50%']
        else:
            labels=['Bottom {}%'.format(int(q*100)), 'Normal', 'Top {}%'.format(int(q*100))]
            
    ls = r.c(*colors)
    
    def plot_me(sub_f, label):
        if (get_vec_type(sub_f) == 'real') and (len(sub_f.unique()) > 10):
            sub_f = to_quants(sub_f, q=q, std=std)
        m = get_cox_ph_ms(surv, sub_f, return_val='model', formula=fmla)
        r_data = m.rx2('call')[2]
        s = survival.survdiff(fmla, r_data)
        p = str(s).split('\n\n')[-1].strip().split(', ')[-1]
        ls = r.c(*colors)
        
        r.plot(survival.survfit(fmla, r_data), lty=1, col=ls, lwd=4, cex=1.25, 
                                xlab='Years to Event', ylab='Survival');
        r.title(label, cex=3.)
        if ann=='p':
            r.text(3, 0, labels='logrank ' + p, pos=4)
        elif ann != None:
            r.text(0, labels=ann, pos=4)

    if show_legend == 'out':  
        r.par(xpd=True, mar=r.c(4,5,5,8))
    for value in sorted(assignment.ix[feature.index].dropna().unique()):
        plot_me(feature.ix[assignment[assignment==value].index], name(value))

    if show_legend == True:
        mean_s = surv.ix[:,'event'].ix[assignment[assignment==value].index].mean()
        if mean_s < .5:
            r.legend(surv.ix[:,'days'].max() * .05 / 365., .45, labels, 
                     lty=1, col=ls, lwd=3, bty='o')
        else:
            r.legend(surv.ix[:,'days'].max() * .4 / 365, .9, labels, 
                     lty=1, col=ls, lwd=3, bty='o')
    elif show_legend == 'out':
        r.legend(surv.ix[:,'days'].max() * 1.1  / 365, .9, labels, 
                     lty=1, col=ls, lwd=3, bty='o')
    r('dev.off()')
    if show:
        return Show(filename)
    
class Show(object):
    def __init__(self, filename):
        self.filename = filename
        self.image = Image(filename=self.filename)
    def _repr_html_(self):
        return ('''<img src='data:image/png;base64,''' + 
                base64.standard_b64encode(self.image.data) + '\'>')
    
def draw_survival_curves_model(feature, test, filename='tmp.png', show=False, 
                               title=True, labels=['No Mutation', 'Mutation'], 
                               colors=['blue','red'], ann=None):
    name = feature.name
    fmla = robjects.Formula('Surv(days, event) ~ feature')
    m = get_cox_ph_ms(test.surv, feature, test.covariates, 
                      return_val='model', formula=fmla)
    ls = r('2:' + str(len(set(feature))+1)) #R line styles
    if filename.endswith('.png'):
        r.png(filename=filename, width=300, height=250, res=100, pointsize=8)
    else:
        r.pdf(filename, width=4.5, height=3.75)
    
    r_data = m.rx2('call')[2]
    s = survival.survdiff(fmla, r_data)
    p = str(s).split('\n\n')[-1].strip().split(', ')[-1]
    ls = r.c(*colors)
    r.plot(survival.survfit(fmla, r_data), lty=1, col=ls, lwd=3, cex=1.15, 
                            xlab='Days to Event', ylab= 'Survival');
    if title:
        r.title(name)
        
    mean_s = test.surv.ix[:,'event'].mean()
    if mean_s < .5:
        r.legend(test.surv.ix[:,'days'].max() * .5, .4, labels, 
                 lty=1, col=ls, lwd=3, bty='o')
    else:
        r.legend(test.surv.ix[:,'days'].max() * .5, .9, labels, 
                 lty=1, col=ls, lwd=3, bty='o')
    r.par(xpd=True)
    if ann=='p':
        r.text(0, labels='logrank ' + p, pos=4)
    elif ann != None:
        r.text(0, labels=ann, pos=4)
    r('dev.off()')
    

def draw_survival_curves_split(feature, assignment, surv, filename='tmp.png', show=False, 
                               title=True, labels=['No Mutation', 'Mutation'], 
                               colors=['blue','red'], ann=None, show_legend=True, q=.25):
    if filename.endswith('.png'):
        r.png(filename=filename, width=200*(len(assignment.unique())+1), height=300, res=75)
    else:
        r.pdf(filename, width=9.5, height=3.75)
        
    fmla = robjects.Formula('Surv(days, event) ~ feature')
    r.par(mfrow=r.c(1,len(assignment.unique())))
    r.par(mar=r.c(4,5,4,1))
    
    if (get_vec_type(feature) == 'real') and (len(feature.unique()) > 5):
        colors=['blue','orange','red']
        if q == .5:
            labels=['Bottom 25%', 'Top 50%']
        else:
            labels=['Bottom {}%'.format(int(q*100)), 'Normal', 'Top {}%'.format(int(q*100))]
            
    ls = r.c(*colors)
    
    def plot_me(sub_f, label):
        if (get_vec_type(sub_f) == 'real') and (len(sub_f.unique()) > 5):
            sub_f = to_quants(sub_f, q=q)
        m = get_cox_ph_ms(surv, sub_f, return_val='model', formula=fmla)
        r_data = m.rx2('call')[2]
        s = survival.survdiff(fmla, r_data)
        p = str(s).split('\n\n')[-1].strip().split(', ')[-1]
        ls = r.c(*colors)
        
        
        r.plot(survival.survfit(fmla, r_data), lty=1, col=ls, lwd=4, cex=1.25, 
                                xlab='Days to Event', ylab='Survival');
        r.title(label, cex=3.)
        if ann=='p':
            r.text(0, labels='logrank ' + p, pos=4)
        elif ann != None:
            r.text(0, labels=ann, pos=4)
            
    for value in sorted(assignment.ix[feature.index].dropna().unique()):
        plot_me(feature.ix[assignment[assignment==value].index], 
                str(assignment.name) + ' = ' + str(value))

    if show_legend:
        mean_s = surv.ix[assignment[assignment==value].index].ix[:,'event'].mean()
        if mean_s < .4:
            r.legend(surv.ix[:,'days'].max() * .05, .4, labels, 
                     lty=1, col=ls, lwd=3, bty='o')
        else:
            r.legend(surv.ix[:,'days'].max() * .4, .9, labels, 
                     lty=1, col=ls, lwd=3, bty='o')
    
    r.par(xpd=True)
    r('dev.off()')
    if show:
        return Show(filename)

    
def draw_pathway_count_bar(p, cancer, file_name='tmp.svg', colors='red'):
    fig, ax = plt.subplots(1,1, figsize=(3+len(cancer.gene_sets[p])/30.,2))
    if 'deceased' in cancer.clinical:
        patients = cancer.clinical[['age', 'deceased']].dropna().index
    else:
        patients = cancer.clinical.index
    m = cancer.hit_matrix.ix[list(cancer.gene_sets[p]), patients] > 0
    m = sort(m[m.sum(1) > 0].sum(1))
    #colors = ['orange']*(len(m))
    m.plot(kind='bar', ax=ax, alpha=.7, color=colors);
    ax.set_yticks(range(m.max()+1))
    ax.set_yticklabels(range(m.max()+1), size=14)
    ax.set_xticklabels(m.index, ha='center', va='bottom', position=(0,.1), size=20)
    ax.set_ylabel('# Patients', size=16)
    
    fig.tight_layout()
    fig.savefig(file_name)
    
def draw_pathway_eig_bar(U, file_name='tmp.svg'):
    if rank(U)  == 2:
        U = U[0]
    fig, ax = plt.subplots(1,1, figsize=(3+len(U)/15.,2.5))
    sort(U).plot(kind='bar', ax=ax, alpha=.4)
    ax.set_ylabel('Loading')
    fig.tight_layout()
    fig.savefig(file_name)
    
def draw_pathway_age_scatter(p, cancer, file_name='tmp.svg'):
    fig, ax = plt.subplots(1,1, figsize=(6,4))
    ax.scatter(*match_series(cancer.pc.ix[p], cancer.clinical.age), alpha=.5, s=75)
    ax.set_xlabel('Principal Component Loading')
    ax.set_ylabel('Age')
    fig.savefig(file_name)
    
def series_scatter(s1, s2, ax=None, ann='p', filename=None, **plot_args):
    if ax is None:
        fig, ax = plt.subplots(1,1, figsize=(6,4))
    ax.scatter(*match_series(s1, s2), alpha=.5, s=75, **plot_args)
    ax.set_xlabel(s1.name)
    ax.set_ylabel(s2.name)
    if ann == 'p':
        ax.annotate('p = {0:.2e}'.format(pearson_p(s1, s2)), (.95, .02),
                    xycoords='axes fraction', ha='right',va='bottom', size=12)
    if filename is not None:
        fig.savefig(filename)    
    
def histo_compare(hit_vec, response_vec, ax=None):
    '''
    Split response_vec by hit_vec and compared histograms.  
    Also plots the kde of the whole response_vec.
    '''
    if ax is None:
        fig, ax = plt.subplots(1,1)
    else:
        fig = fig.gcf()
    kde1 = gaussian_kde(response_vec)
    x_eval = linspace(min(response_vec), max(response_vec), num=200)
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
    x_order = sort(-s[s > 0]).index
    y_order = sort(mut.meta_matrix.ix[p].sum(1)).index
    plt.imshow(mut.meta_matrix.ix[y_order, x_order], aspect=8, 
               interpolation='Nearest', cmap=plt.cm.bone_r) #@UndefinedVariable
    ax.set_yticks(range(len(p)))
    ax.set_yticklabels(p);
    ax.hlines(arange(len(p)) + .5, 0, len(x_order), colors='grey', alpha=.3, lw=4)
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
    for (x,y),w in ndenumerate(W):
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
    o = sort(df.apply(lambda s: ''.join(map(str, s)))).index[::-1]
    df = df[o]
    
    if df.shape[0] > 20:
        rest = Series(df.ix[10:].sum().clip_upper(1.), name='rest')
        df = df.ix[:10]
        df = df.append(rest)
    if ax is None:
        fig, ax = plt.subplots(figsize=(df.shape[1]*.2,df.shape[0]*.5))
    else:
        fig = ax.get_figure()
    memo_plot(df, ax=ax)
    ax.bar(arange(len(df.columns)) - .3, df.sum() / df.sum().max(), bottom=-1.5, 
           width=.6, alpha=.5)
    counts = df.sum(1)[::-1]
    width = df.shape[1]
    ax.barh(arange(len(counts)) - .3, (counts / counts.max())*width*.25, left=width - .2, 
            height=.6, alpha=.5)
    ax.set_frame_on(False)
    ax.tick_params(right='off')
    fig.tight_layout()
    
def draw_pathway_overlaps(mat, bars, filename=None):  
    fig, ax = plt.subplots(figsize=(25,5))   
    memo_plot(mat, ax=ax)
    ax.bar(arange(len(mat.columns)) - .3, bars, bottom = -2, width=.6, alpha=.5)
    ax.hlines(-1, -.5, len(mat.columns))
    ax.annotate('Days to Death', (-.5,-1.5), ha='right', va='center', size=15)
    fig.tight_layout()
    if filename is not None:
        fig.savefig(filename)
    return fig



def fancy_raster(df, cluster=False, cmap=cm.get_cmap('Spectral'), norm=None):
    if cluster:
        d = dist.pdist(df)
        D = dist.squareform(d)
        Y = sch.linkage(D)
        Z = sch.dendrogram(Y, no_plot=True)
        order = Z['leaves']
        df = df.ix[order, order]
        
    fig, ax = plt.subplots(1,1, figsize=(12,8))
    img = ax.imshow(df, interpolation='Nearest', cmap=cmap, norm=norm)
    ax.set_yticks(range(len(df.index)))
    ax.set_yticklabels(df.index)
    ax.set_xticks(arange(len(df.columns)))
    ax.set_xticklabels(df.columns, rotation=360-90, ha='center');
    ax.hlines(arange(len(df.index)-1)+.5, -.5, len(df.columns)-.5, 
              color='white', lw=6)
    ax.vlines(arange(len(df.columns)-1)+.5, -.5, len(df.index)-.5, 
              color='white', lw=6)
    
    if cluster:
        icoord = array(Z['icoord']) - array(Z['icoord']).min()
        icoord = icoord * ((len(Z['leaves']) - 1) / icoord.max())
    
        dcoord = -1*array(Z['dcoord']) - .7 
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

    

