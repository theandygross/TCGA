'''
Created on Jun 12, 2013

@author: agross
'''
import matplotlib.pyplot as plt

def init_ax(ax, figsize=None):
    '''
    Helper to initialize an axis.  If the axis is not specified (eg is None),
    this fill create a new figure. 
    
    ax: matplotlib axis object, or None
    figsize: size of figure if we have to create a new axis.
    '''
    if ax is None:
        if figsize is None:
            fig, ax = plt.subplots(1,1)
        else:
            fig, ax = plt.subplots(1,1, figsize=figsize)
    else:
        fig = plt.gcf()
    return fig, ax

def prettify_ax(ax):
    ax.grid(b=False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
def latex_float(f):
    '''
    http://stackoverflow.com/questions/13490292/
    format-number-using-latex-notation-in-python
    '''
    float_str = "{0:.2g}".format(f)
    if "e" in float_str:
        base, exponent = float_str.split("e")
        return r"{0} \times 10^{{{1}}}".format(base, int(exponent))
    else:
        return float_str
