__author__ = 'agross'

import pandas as pd

from pandas.core.frame import com, StringIO
from IPython.core.display import Javascript, display

old_rep = pd.DataFrame._repr_html_


def hover_table(self):
    """
    Return a html representation for a particular DataFrame.
    Mainly for IPython notebook.
    """
    # ipnb in html repr mode allows scrolling
    # users strongly prefer to h-scroll a wide HTML table in the browser
    # then to get a summary view. GH3541, GH3573
    ipnbh = com.in_ipnb() and pd.get_option('display.notebook_repr_html')

    # qtconsole doesn't report it's line width, and also
    # behaves badly when outputting an HTML table
    # that doesn't fit the window, so disable it.
    if com.in_qtconsole():
        raise NotImplementedError('HTML output is disabled in QtConsole')

    if self._info_repr():
        buf = StringIO(u(""))
        self.info(buf=buf)
        return '<pre>' + buf.getvalue() + '</pre>'

    if pd.get_option("display.notebook_repr_html"):
        max_rows = pd.get_option("display.max_rows")
        max_cols = pd.get_option("display.max_columns")
        html = self.to_html(max_rows=max_rows,
                            max_cols=max_cols,
                            show_dimensions=False,
                            classes='table table-hover')

        text =  '<div style="max-height:1000px; max-width:900px;overflow:auto;">\n'
        text += html
        text += '\n</div>'
        text = text.replace('dataframe ','')
        return text
    else:
        return None


def fancy_table(df, div_id):
    def _repr_html_(self):
        if self._info_repr():
            buf = StringIO(u(""))
            self.info(buf=buf)
            return '<pre>' + buf.getvalue() + '</pre>'


        max_rows = pd.get_option("display.max_rows")
        max_cols = pd.get_option("display.max_columns")
        html = self.to_html(max_rows=max_rows,
                            max_cols=max_cols,
                            show_dimensions=False,
                            classes='table table-bordered table-striped')

        text =  '<div>\n'
        text += html
        text += '\n</div>'
        #text = text.replace('border="1" ','border="2" ')
        text = text.replace('dataframe ','')
        text = text.replace('<table ','<table id="{}" '.format(div_id))
        return text

    pd.DataFrame._repr_html_ = _repr_html_
    display(df)
    pd.DataFrame._repr_html_ = old_rep
    #display(Javascript("$('#{}').dataTable();".format(div_id)))

