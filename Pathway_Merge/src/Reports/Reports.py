'''
Created on Oct 30, 2012

@author: agross
'''
import os as os
import pandas.rpy.common as com
from pandas import Series, DataFrame
from numpy import nan, roll
from rpy2.robjects.packages import importr

from Figures import *

from Processing.Helpers import frame_svd, match_series

nz = importr('Nozzle.R1')
bool_ = {True: 'TRUE', False: 'FALSE'}
FIG_EXT = 'clinical_figures/'
SORT_ORDER = ['event_free_survival','AMAR','age','gender','radiation','therapy']

roll_df = lambda df, num: df.ix[:,roll(range(df.shape[1]),num)]

def create_figure(cancer, type, vec, file_name, overwrite=False):
    if (overwrite is False) and os.path.isfile(file_name):
        return
    if type in cancer.survival_tests:
        draw_survival_curves(cancer.clinical, vec.clip_upper(1.), 
                             filename=file_name, **cancer.survival_tests[type])
    elif type in cancer.real_variables:
        violin_plot_pandas(vec > 0, cancer.clinical[type].astype(float), 
                           filename=file_name)
    elif type in cancer.binary_variables:
        fischer_bar_chart(vec > 0, cancer.clinical[type], filename=file_name)
        
    elif type == 'pathway_bar':
        draw_pathway_count_bar(vec.name, cancer, cancer.gene_sets, file_name)

class NozzleTable(object):
    def __init__(self, table, path, table_file_name, caption, cutoff=.25,
                 idx_name='', tests_run=[]):
        self.path = path
        self.cutoff = cutoff
        
        table = table.sort(columns=[s for s in SORT_ORDER if s in table
                                    and len(table[s].unique()) > 1])
        table.to_csv(path + table_file_name)
        keepers = table[(table[tests_run] < cutoff).sum(1) > 0].index
        table = table.ix[keepers].head(20)
        
        index_series = Series(dict((i,i) for i in table.index), name=idx_name)
        self.table = table.join(index_series)
        self.table = roll_df(self.table,1)
        self.tests_run = tests_run
        self.row_pos = dict((g,i+1) for i,g in enumerate(self.table.index))
        self.col_pos = dict((c,i+1) for i,c in enumerate(self.table.columns))
        
        r_table = self.table.copy().fillna(1.)
        r_table = com.convert_to_r_dataframe(r_table) #@UndefinedVariable
        self.nozzle_table = nz.newTable(r_table, caption, file=table_file_name, 
                                        significantDigits=2)
        
    def add_to_table(self, fig_files, row, col, captions='', title='',
                     isSignificant=True):
        '''
        Add a figure from a file to the Nozzle table.
        This is just a wrapper to put a figure(s) in a result section.
        '''
        if type(fig_files) == str:
            fig_files = [fig_files]
            captions = [captions]
        elif len(captions) != len(fig_files):
            captions = [captions]*len(fig_files)
        figs = [nz.newFigure(file, caption) for file, caption 
                in zip(fig_files, captions)]
        result = nz.addTo(nz.newResult('', isSignificant=bool_[isSignificant]),
                          nz.addTo(nz.newSubSection(title), *figs))
        self.nozzle_table = nz.addTo(self.nozzle_table, result, 
                                     row=self.row_pos[row], 
                                     column=self.col_pos[col])                
 
class PathwayTable(NozzleTable):
    def __init__(self, cancer, cutoff=.25):
        path = cancer.report_folder + '/'
        caption = ('Association of pathway level ' + cancer.data_type + 
                   ' patients with patient clinical features')
        tests_run = list(cancer.q_pathways.columns)
        
        annotation = DataFrame(dict((p,get_pathway_annotation_vec(p, cancer)) 
                       for p in cancer.q_pathways.index)).T
        pathway_table = annotation.join(cancer.q_pathways)
        pathway_table = pathway_table[pathway_table['n mut pat'] > 3]
        
        NozzleTable.__init__(self, pathway_table, path, 'pathway_table.csv', 
                             caption, cutoff, 'pathways', tests_run)
        
    def fill_in_table(self, cancer, data_frame):
        for p,vec in self.table.iterrows():
            file_name = self.path + FIG_EXT + p + '.png'
            create_figure(cancer, 'pathway_bar', vec, file_name)
            
        for test,vals in self.table[self.tests_run].iteritems():
            t_caption  = (test.capitalize() + ' association for patients ' + 
                        'with or without ' + cancer.data_type + ' to genes ' + 
                        'in pathway.')
            b_caption = ('Distribution of ' + cancer.data_type + 's to genes' + 
                         ' in the pathway.')
            for g in vals[vals < self.cutoff].index: #@NoEffect, does have effect
                title = test + ' test for ' + g
                fig_file = self.path + FIG_EXT + g + '_' + test + '.png'
                create_figure(cancer, test, data_frame.ix[g], fig_file)
                path_file = self.path + FIG_EXT + g + '.png'          
                self.add_to_table([fig_file, path_file], g, test, 
                                  [t_caption, b_caption], title)
        
        
class HitTable(NozzleTable):
    def __init__(self, cancer, cutoff=.25):
        path = cancer.report_folder + '/'
        self.marker = ('region' if cancer.data_type == 'amplification' 
                       else 'gene')
        caption = ('Association of single ' + self.marker + ' ' + 
                   cancer.data_type + ' with patient clinical features.')
        
        hit_matrix = (cancer.lesion_matrix if (cancer.data_type == 
                      'amplification') else cancer.hit_matrix)
        hit_matrix = hit_matrix.groupby(level=0).first() #Make index unique
        tests_run = list(cancer.q_genes.columns)
        counts = Series((hit_matrix.ix[:,cancer.patients] > 0).sum(1), 
                        name='n_patients')
        gene_table = cancer.q_genes.join(counts)
        gene_table = roll_df(gene_table, 1)
        
        NozzleTable.__init__(self, gene_table, path, 'gene_table.csv', caption, 
                             cutoff, self.marker, tests_run)
        
    def fill_in_table(self, cancer, data_frame):
        for test,vals in self.table[self.tests_run].iteritems():
            caption  = (test.capitalize() + ' association for patients ' + 
                        'with or without ' + cancer.data_type + ' to ' + 
                        self.marker + '.') 
            for g in vals[vals < self.cutoff].index: #@NoEffect, does have effect
                fig_file = self.path + FIG_EXT + g + '_' + test + '.png'
                title = test + ' test for ' + g
                create_figure(cancer, test, data_frame.ix[g], fig_file)
                self.add_to_table(fig_file, g, test, caption, title)
        
def get_pathway_annotation_vec(p, cancer):
    genes = cancer.gene_sets[p].intersection(set(cancer.hit_matrix.index))
    sub_matrix = cancer.hit_matrix.ix[genes, cancer.patients] > 0
    annot_vec = Series({'n mut genes' :  (sub_matrix.sum(1) > 0).sum(),
                        'n mut pat' : (sub_matrix.sum() > 0).sum(),
                        'n genes' : len(cancer.gene_sets[p])})
    return annot_vec

def generic_header(report, cancer, prev_cancer, next_cancer):
    report = nz.setMaintainerName(report, 'Andrew Gross')
    report = nz.setMaintainerEmail(report, "agross@ucsd.edu" );
    report = nz.setMaintainerAffiliation(report, 'UCSD- Bioinf. and '
                                                + 'Systems Biology' );
    next_file = cancer.report_folder.replace(cancer.cancer, next_cancer)
    prev_file = cancer.report_folder.replace(cancer.cancer, prev_cancer)                                           
    report = nz.setPreviousReport(report, prev_file  + '/index.html')
    report = nz.setNextReport(report, next_file + '/index.html')
    return report