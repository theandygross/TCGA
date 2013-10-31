'''
Created on Jun 30, 2013

@author: agross
'''
import numpy as np
import pandas as pd

import matplotlib.pylab as plt
from matplotlib.patches import FancyBboxPatch

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
    for (x, y), w in np.ndenumerate(W):
        if w == 0: 
            rect = null_box(x, y)
        else:     
            rect = null_box(x, y)
            ax.add_patch(rect)
            rect = mut_box(x, y) 
        ax.add_patch(rect)
    df
    ax.set_xbound(-.5, len(W) - .5)
    ax.set_ybound(-.5, len(W[0]) - .5)
    
def pathway_plot(df, ax=None, bar='both'):
    df = df.ix[df.sum(1) > 0, df.sum() > 0]
    df = df.ix[df.sum(1).order(ascending=False).index]
    o = df.astype(int).apply(lambda s: ''.join(map(str, s))).order()
    df = df[o.index[::-1]]
    
    if df.shape[0] > 20:
        rest = pd.Series(df.ix[10:].sum().clip_upper(1.), name='rest')
        df = df.ix[:10]
        df = df.append(rest)
    if ax is None:
        fig, ax = plt.subplots(figsize=(df.shape[1] * .2, df.shape[0] * .5))
    else:
        fig = ax.get_figure()
    memo_plot(df, ax=ax)
    if bar in ['x', 'both']:
        ax.bar(np.arange(len(df.columns)) - .3, 1.*df.sum() / df.sum().max(),
               bottom= -1.5, width=.6, alpha=.5)
    if bar in ['y', 'both']:
        counts = df.sum(1)[::-1]
        width = df.shape[1]
        ax.barh(np.arange(len(counts)) - .3, (1.*counts / counts.max()) * width * .25,
                left=width - .2, height=.6, alpha=.5)
        
    ax.set_frame_on(False)
    ax.tick_params(right='off')
    fig.tight_layout()
    return df
    
def draw_pathway_overlaps(mat, bars, filename=None):  
    fig, ax = plt.subplots(figsize=(25, 5))   
    memo_plot(mat, ax=ax)
    ax.bar(np.arange(len(mat.columns)) - .3, bars, bottom= -2, width=.6,
           alpha=.5)
    ax.hlines(-1, -.5, len(mat.columns))
    ax.annotate('Days to Death', (-.5, -1.5), ha='right', va='center', size=15)
    fig.tight_layout()
    if filename is not None:
        fig.savefig(filename)
    return fig
