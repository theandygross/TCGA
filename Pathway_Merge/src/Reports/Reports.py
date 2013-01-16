'''
Created on Oct 30, 2012

@author: agross
'''
import os as os
from pandas import Series
from numpy import roll
from rpy2.robjects.packages import importr

from Figures import *
from Processing.Helpers import frame_svd


nz = importr('Nozzle.R1')
bool_ = {True: 'TRUE', False: 'FALSE'}
FIG_EXT = 'clinical_figures/'
SORT_ORDER = ['event_free_survival_5y', 'survival_5y', 'tumor_t1t2', 'lymphnode_n0n1',
              'metastatic_recurrence', 'organ_subdivision', 'AMAR','age',
              'gender','radiation', 'therapy', 'rate_non']
DATA_TYPE = {'mutation': 'bool', 'amplification': 'bool', 'deletion': 'bool',
             'methylation': 'real', 'expression': 'real', 
             'expression_array': 'real'}

roll_df = lambda df, num: df.ix[:,roll(range(df.shape[1]),num)]

def generic_header(report, cancer, prev_cancer, next_cancer, 
                   address='index.html'):
    report = nz.setMaintainerName(report, 'Andrew Gross')
    report = nz.setMaintainerEmail(report, "agross@ucsd.edu" );
    report = nz.setMaintainerAffiliation(report, 'UCSD- Bioinf. and '
                                                + 'Systems Biology' ); 
    '''Little wonky due to relative paths but works'''
    page_ext = cancer.report_folder.split('/')[2:-1] +  [address]                                       
    next_page = '../../../' + '/'.join([next_cancer] + page_ext)  
    prev_page = '../../../' + '/'.join([prev_cancer] + page_ext)                    
    report = nz.setPreviousReport(report, next_page, next_cancer)
    report = nz.setNextReport(report, prev_page, prev_cancer)
    return report

def draw_clinical_figs(cancer):
    fig_dir = cancer.data_path + cancer.report_folder + FIG_EXT
    for var, vals in cancer.clinical.iteritems():
        fig_file = fig_dir + var + '.png'
        if os.path.isfile(fig_file) or (vals.notnull().sum() < 5):
            continue
        fig, ax = plt.subplots()
        if (vals.dtype in ['object', 'bool']) or (len(set(vals)) < 3):
            vals.value_counts().plot(kind='bar', title=var, ax=ax)
            ax.set_ylabel('Number of Patients')
        elif vals.dtype in ['float', 'float64', 'int', 'int64']:
            vals.hist(ax=ax)
            ax.set_xlabel(var)
        fig.savefig(fig_file)
        
def create_clinical_overview(cancer):
    fig_dir = cancer.data_path + cancer.report_folder + FIG_EXT
    if not os.path.isdir(fig_dir):
        os.makedirs(fig_dir)
    results = []
    for var in cancer.clinical:
        f = fig_dir + var + '.png'
        res = nz.newResult(var, isSignificant=False)
        res = nz.addTo(res, nz.addTo(nz.newSection(var), nz.newFigure(f)))
        results.append(res)
    l = [['\t',nz.asSummary(r)] for r in results]
    l = [b for a in l for b in a]
    par = nz.newParagraph(*l)
    clin = nz.addTo(nz.newSubSection('Clinical Variables'), par)
    return clin

def desc(s):
    return Series({'FDR.25' : sum(s < .25),
                   'Best' : s.idxmin(),
                   'Best q-value' : s.min()})

def make_nz_table(table, idx_name, caption):
    index_series = Series(dict((i,i) for i in table.index), name=idx_name)
    table = table.join(index_series)
    table = roll_df(table,1)
    r_table = table.copy().fillna(1.)
    r_table = com.convert_to_r_dataframe(r_table)
    nozzle_table = nz.newTable(r_table, caption, significantDigits=2)
    return nozzle_table

def get_association_overview(cancer):
    association_overview = nz.newSubSection('Association Overview')
    if len(cancer.q_genes) > 0:
        gene_summary = cancer.q_genes.apply(desc).T
        caption = ('Summary of clinical associations with ' + cancer.data_type)
        gene_summary = make_nz_table(gene_summary, 'test', caption)
        association_overview = nz.addTo(association_overview, gene_summary)
    
    path_summary = cancer.q_pathways.apply(desc).T
    caption = ('Summary of clinical associations with ' + cancer.data_type)
    path_summary = make_nz_table(path_summary, 'test', caption)
    association_overview = nz.addTo(association_overview, path_summary)
    return association_overview

def create_figure(cancer, fig_type, vec, file_name, overwrite=False):
    if (overwrite is False) and os.path.isfile(file_name):
        return 
    if DATA_TYPE[cancer.data_type] == 'bool':
        create_figure_bool(cancer, fig_type, vec, file_name)
    else:
        create_figure_real(cancer, fig_type, vec, file_name)
    
def create_figure_bool(cancer, fig_type, vec, file_name):
    if fig_type in cancer.survival_tests:
        draw_survival_curves(cancer.clinical, vec.clip_upper(1.), 
                             filename=file_name, **cancer.survival_tests[fig_type])
    elif fig_type in cancer.real_variables:
        violin_plot_pandas(vec > 0, cancer.clinical[fig_type].astype(float), 
                           filename=file_name)
    elif fig_type in cancer.binary_variables:
        fischer_bar_chart(vec > 0, cancer.clinical[fig_type], filename=file_name)
    elif fig_type == 'pathway_bar':
        draw_pathway_count_bar(vec.name, cancer, file_name)
        
def create_figure_real(cancer, fig_type, vec, file_name):
    if fig_type in cancer.survival_tests:
        hit_vec = -1*(vec < -1) + (vec > 1)
        draw_survival_curves(cancer.clinical, hit_vec, filename=file_name, 
                             labels=['low','normal','high'],
                             **cancer.survival_tests[fig_type])
    elif fig_type in cancer.real_variables:
        series_scatter(vec, cancer.clinical[fig_type].astype(float), 
                       filename=file_name)
    elif fig_type in cancer.binary_variables:
        violin_plot_pandas(cancer.clinical[fig_type], vec, filename=file_name)
        
    elif fig_type == 'pathway_bar':
        genes = cancer.gene_sets[vec.name]
        U,S,vH = frame_svd(cancer.data_matrix.ix[genes].dropna())
        draw_pathway_eig_bar(U, file_name)

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
        path = cancer.data_path + cancer.report_folder + '/'
        caption = ('Association of pathway level ' + cancer.data_type + 
                   ' patients with patient clinical features')
        tests_run = list(cancer.q_pathways.columns)
        
        if DATA_TYPE[cancer.data_type] == 'bool':
            annotation = DataFrame(dict((p,get_pathway_annotation_vec(p, cancer)) 
                           for p in cancer.q_pathways.index)).T
            pathway_table = annotation.join(cancer.q_pathways)
            pathway_table = pathway_table[pathway_table['n mut pat'] > 3]
        else:
            pathway_size = Series({p: len(cancer.gene_sets[p]) for p 
                              in cancer.q_pathways.index},name='n genes')
            pathway_table = DataFrame(pathway_size).join(cancer.q_pathways)
        
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
        path = cancer.data_path + cancer.report_folder + '/'
        self.marker = ('region' if cancer.data_type == 'amplification' 
                       else 'gene')
        caption = ('Association of single ' + self.marker + ' ' + 
                   cancer.data_type + ' with patient clinical features.')
        
        hit_matrix = (cancer.lesion_matrix if (cancer.data_type == 
                      'amplification') else cancer.hit_matrix)
        hit_matrix = hit_matrix.groupby(level=0).first()
        tests_run = list(cancer.q_genes.columns)
        if 'deceased' in cancer.clinical:
            patients = cancer.clinical[['age', 'deceased']].dropna().index
        else:
            patients = cancer.clinical.columns
        counts = Series((hit_matrix.ix[:, patients] > 0).sum(1), 
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
    if 'deceased' in cancer.clinical:
        patients = cancer.clinical[['age', 'deceased']].dropna().index
    else:
        patients = cancer.clinical.index
    sub_matrix = cancer.hit_matrix.ix[genes, patients] > 0
    annot_vec = Series({'n mut genes' :  (sub_matrix.sum(1) > 0).sum(),
                        'n mut pat' : (sub_matrix.sum() > 0).sum(),
                        'n genes' : len(cancer.gene_sets[p])})
    return annot_vec

def create_clinical_report(cancer, next_cancer, prev_cancer):
    if not os.path.isdir(cancer.report_folder + '/' + FIG_EXT):
        os.makedirs(cancer.report_folder + '/' + FIG_EXT)
    report = nz.newReport('Report for ' + cancer.cancer)
    report = generic_header(report, cancer, next_cancer, prev_cancer)
    
    clin = create_clinical_overview(cancer)
    draw_clinical_figs(cancer)
    no_association = (cancer.q_genes < .5).sum().add(
                     (cancer.q_pathways < .5).sum(), fill_value=0) == 0
    no_a = nz.newParagraph(('No associations obtained for: \n' + 
                 ', '.join(no_association[no_association].index)))
    report = nz.addToSummary(report, no_a)
    report = nz.addToSummary(report, clin)
    
    cancer.q_genes = cancer.q_genes.ix[:,no_association==False]
    cancer.p_genes = cancer.p_genes.ix[:,no_association==False]
    cancer.q_pathways = cancer.q_pathways.ix[:,no_association==False]
    cancer.p_pathways = cancer.p_pathways.ix[:,no_association==False]
    report = nz.addToSummary(report, get_association_overview(cancer))
    
    if hasattr(cancer, 'q_genes') and len(cancer.q_genes.dropna()) > 0:
        gene_table = HitTable(cancer)
        df = (cancer.hit_matrix if gene_table.marker == 'gene' else
              cancer.lesion_matrix)
        df = df.groupby(level=0).first()
        #df = cancer.hit_matrix.groupby(level=0).first()
        gene_table.fill_in_table(cancer, df)
        section = nz.addTo(nz.newSubSection(gene_table.marker.capitalize() + ' ' + 
                           cancer.data_type), gene_table.nozzle_table)
        report = nz.addToResults(report, section)
    
    pathway_table = PathwayTable(cancer)
    if len(pathway_table.table) > 0:
        if hasattr(cancer, 'meta_matrix'):
            meta_matrix = cancer.meta_matrix
        else:
            meta_matrix = cancer.pc
        pathway_table.fill_in_table(cancer, meta_matrix)
        section = nz.addTo(nz.newSubSection('Pathway ' + cancer.data_type), 
                           pathway_table.nozzle_table)
        report = nz.addToResults(report, section)
        
    nz.writeReport(report, filename=cancer.data_path + cancer.report_folder + 'index')