'''
Created on Jun 30, 2013

@author: agross
'''
import numpy as np
import rpy2.robjects as robjects
from rpy2.robjects import r
from Processing.Helpers import get_vec_type ,to_quants
from Stats.Survival import get_cox_ph, log_rank, survival
from Reports.NotebookTools import Show
    
def draw_survival_curves(feature, surv, assignment=None, filename='tmp.png', show=False, 
                        title=True, labels=None, colors=['blue','red'], ann=None, 
                        show_legend=True, q=.25, std=None):
    if assignment is None:
        num_panels = 1
        assignment = feature.map(lambda s: 1)
        name = lambda v: str(feature.name) if feature.name != None else ''
    else:
        num_panels = len(assignment.unique())
        name = lambda v: str(assignment.name) + ' = ' + str(v)
    if (labels is None) and ((len(feature) / feature.nunique()) > 10):
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
            
        m = get_cox_ph(surv, sub_f, formula=fmla)
        r_data = m.rx2('call')[2]
        p = log_rank(sub_f, surv)['p']
        ls = r.c(*colors)
        
        r.plot(survival.survfit(fmla, r_data), lty=1, col=ls, lwd=4, cex=1.25, 
                                xlab='Years to Event', ylab='Survival');
        r.title(label, cex=3.)
        if ann=='p':
            r.text(.2, 0, labels='logrank p = {0:.1e}'.format(p), pos=4)
        elif ann != None:
            r.text(0, labels=ann, pos=4)

    if show_legend == 'out':  
        r.par(xpd=True, mar=r.c(4,5,5,8))
    for value in sorted(assignment.ix[feature.index].dropna().unique()):
        f = feature.ix[assignment[assignment==value].index]
        if len(f.unique()) > 1:
            plot_me(f, name(value))

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