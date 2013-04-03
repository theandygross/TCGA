'''
Created on Mar 6, 2013

@author: agross
'''
import base64
from IPython.display import Image

class stack(object):
    def __init__(self, objs):
        self.objs = objs
    def _repr_html_(self):
        s = '''<div style='border:2px solid #a1a1a1;padding:5px 20px;border-radius:10px;float:left'>'''
        for obj in self.objs:
            s += '<p>{}</p>'.format(obj._repr_html_())
        s += '</div>'
        return s
    
class side_by_side2(object):
    def __init__(self, objs):
        self.objs = objs
    def _repr_html_(self):
        s = '<div >'
        for obj in self.objs:
            s += obj._repr_html_()
        s += '</div>'
        return s
    
class side_by_side(object):
    def __init__(self, objs):
        self.images = [Image(filename=f) if (type(f) == str) else f for f in objs]
    def _repr_html_(self):
        imgs = ['''<img src='data:image/png;base64,''' + base64.standard_b64encode(i.data) + '\'>' for i in self.images]
        return ''.join(imgs)
    
class fig_tab(object):
    '''Temporary fix'''
    def __init__(self, fig, tab):
        self.fig = fig
        self.tab = tab
    def _repr_html_(self):
        return ('<div style="float:left;">' + self.fig._repr_html_() + '</div>' + 
                '<div style="float:left; padding-top:50px">' + 
                self.tab._repr_html_()  + '</div>')
        
class RepHTML(object):
    def __init__(self, obj):
        self.obj = obj
    def _repr_html_(self):
        return ('''<p style='align: center'><h3>{}</h3></p>'''.format(str(self.obj)))
    
class Show(object):
    def __init__(self, filename):
        self.filename = filename
        self.image = Image(filename=self.filename)
    def _repr_html_(self):
        return ('''<img src='data:image/png;base64,''' + 
                base64.standard_b64encode(self.image.data) + '\'>')