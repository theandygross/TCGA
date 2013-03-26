'''
Created on Mar 6, 2013

@author: agross
'''
import base64
from IPython.display import Image

class side_by_side(object):
    def __init__(self, fnames):
        self.images = [Image(filename=f) for f in fnames]
    def _repr_html_(self):
        imgs = ['''<img src='data:image/png;base64,''' + base64.standard_b64encode(i.data) + '\'>' for i in self.images]
        return ''.join(imgs)
    
class stack(object):
    def __init__(self, objs):
        self.objs = objs
    def _repr_html_(self):
        s = '<div >'
        for obj in self.objs:
            s += '<p>{}</p>'.format(obj._repr_html_())
        return s
    
class fig_tab(object):
    '''Temporary fix'''
    def __init__(self, fig, tab):
        self.fig = fig
        self.tab = tab
    def _repr_html_(self):
        return ('<div style="float:left;">' + self.fig._repr_html_() + '</div>' + 
                '<div style="float:left; padding-top:50px">' + 
                self.tab._repr_html_()  + '</div>')