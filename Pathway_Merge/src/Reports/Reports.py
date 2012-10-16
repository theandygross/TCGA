'''
Created on Sep 6, 2012

@author: agross
'''
import os as os
import matplotlib.pyplot as plt
from pandas import Series
from numpy import nan
from rpy2.robjects.packages import importr
import pandas.rpy.common as com
from Figures import violin_plot_pandas
from Processing.Helpers import frame_svd, match_series

from Figures import draw_survival_curves, draw_pathway_count_bar
from Figures import draw_pathway_age_scatter, draw_pathway_eig_bar

nz = importr('Nozzle.R1')
bool_ = {True: 'TRUE', False: 'FALSE'}
FIG_EXT = 'clinical_figures/'

def generic_header(report, cancer, prev_cancer, next_cancer, folder):
    report = nz.setMaintainerName(report, 'Andrew Gross')
    report = nz.setMaintainerEmail(report, "agross@ucsd.edu" );
    report = nz.setMaintainerAffiliation(report, 'UCSD- Bioinf. and '
                                                + 'Systems Biology' );
    report = nz.setPreviousReport(report, prev_cancer.data_path + 
                                  folder + '/index.html')
    report = nz.setNextReport(report, next_cancer.data_path + 
                              folder + '/index.html')
    return report

def add_violin_plot(vec, cancer, table1, pos, fig_path):
    fig_file = fig_path + vec.name + '_age.png'
    fig, ax = plt.subplots(1,1)
    violin_plot_pandas(vec > 0, cancer.clinical.age.astype(float), ax=ax)
    try:
        fig.savefig(fig_file)
    except:
        fig
    age_fig1 = nz.newFigure(fig_file, 'Age of patients with or without ' + 
                                      'mutation to gene.')
    result1 = nz.addTo(nz.newResult('', isSignificant='TRUE'),
                       nz.addTo(nz.newSection(vec.name), age_fig1))
    table1 = nz.addTo(table1, result1, row=pos[0], column=pos[1])   
    return table1
    
def add_survival_curve(vec, cancer, table1, pos, fig_path):
    fig_file = fig_path + vec.name + '_survival.png'
    draw_survival_curves(cancer.clinical, vec.clip_upper(1.), filename=fig_file, title=None)
    sv_fig1 = nz.newFigure(fig_file, 'Survival of patients with or ' + 
                                      'without mutation to gene.')
    result1 = nz.addTo(nz.newResult('', isSignificant='TRUE'),
                       nz.addTo(nz.newSection(vec.name), sv_fig1))
    table1 = nz.addTo(table1, result1, row=pos[0], column=pos[1])
    return table1
    
def single_gene_section(cancer, hit_matrix, cutoff=.25, folder='/'):
    #Format data for report
    path = cancer.data_path + folder + '/'
    gene_table_file = path + 'gene_table.csv'
    hit_matrix = hit_matrix.groupby(level=0).first() #Make index unique
    counts = (hit_matrix.ix[:,cancer.patients] > 0).sum(1)
    counts.name = 'n_patients'
    genes = Series(dict((i,i) for i in cancer.q_genes.index), name='gene')
    gene_table = cancer.q_genes.join(counts).join(genes)
    gene_table = gene_table.ix[:,::-1]
    if 'survival' in gene_table:
        gene_table = gene_table.sort(columns='survival')
    gene_table.to_csv(gene_table_file)
    genes_to_show = cancer.q_genes[(cancer.q_genes < .2).sum(1) > 0].index
    gene_table = gene_table.ix[genes_to_show]
    if 'survival' in gene_table:
        gene_table = gene_table.sort(columns='survival')
    gene_table = gene_table.head(20)
    gene_table_r = com.convert_to_r_dataframe(gene_table) #@UndefinedVariable
    
    if len(gene_table) == 0:
        return nz.addTo(nz.newSubSection('Gene Mutations'), nz.newParagraph(''))
    
    #Overview
    tableCaption1 = "Association of gene mutations with patient clinical features."
    table1 = nz.newTable(gene_table_r, tableCaption1, file=gene_table_file, 
                         significantDigits=2);
    #Fill in the details
    gene_pos = dict((g,i+1) for i,g in enumerate(gene_table.index))
    col_pos = dict((c,i+1) for i,c in enumerate(gene_table.columns))
    
    #age violin plots
    if 'age' in gene_table:
        for g,val in gene_table['age'].iteritems():
            num_genes = (match_series(hit_matrix.ix[g], cancer.clinical.age)[0] > 0).sum()
            if val < cutoff and num_genes > 2:
                table1 = add_violin_plot(hit_matrix.ix[g], cancer, table1, 
                                         (gene_pos[g], col_pos['age']),
                                         path + FIG_EXT)
        
    #survival curves
    if 'survival' in gene_table:
        for g,val in gene_table['survival'].iteritems():
            if val < cutoff:
                table1 = add_survival_curve(hit_matrix.ix[g], cancer, table1, (gene_pos[g], 
                                            col_pos['survival']), path + FIG_EXT) 
    
    section = nz.addTo(nz.newSubSection('Gene Mutations'), table1)
    return section

def format_pathway_table(cancer, gene_sets):
    pathways = Series(dict((i,i) for i in cancer.q_pathways.index), name='pathway')
    n_mut_pat = Series((cancer.meta_matrix > 0).sum(1), name='n mut patients')
    p_genes = lambda p: ((cancer.hit_matrix.ix[gene_sets[p], 
                                              cancer.patients] > 0).sum(1) > 0).sum()
    n_mut_genes = Series(dict((p, p_genes(p)) for p in cancer.meta_matrix.index),
                         name='n mut genes')
    n_genes = Series(dict((p, len(gene_sets[p])) for p in cancer.meta_matrix.index),
                     name='n genes')
    pathway_table = cancer.q_pathways.join(n_genes).join(n_mut_genes).join(n_mut_pat)
    pathway_table = pathway_table.join(pathways)
    pathway_table = pathway_table.ix[:, ::-1]
    if 'survival' in pathway_table:
        pathway_table.sort(columns='survival')
    return pathway_table

def add_survival_curve_pathway(vec, cancer, table1, pos, fig_path):
    fig_file = fig_path + vec.name + '_survival.png'
    draw_survival_curves(cancer.clinical, vec.clip_upper(1), filename=fig_file, title=None)
    sv_fig1 = nz.newFigure(fig_file, 'Survival of patients with or ' + 
                                      'without mutation to gene.')
    fig_file2 = fig_path + vec.name + '.svg'
    draw_pathway_count_bar(vec.name, cancer, cancer.gene_sets, fig_file2)
    sv_fig_2 = nz.newFigure(fig_file2, 'Gene mutation frequencies in pathway.')
    result1 = nz.addTo(nz.newResult('', isSignificant='TRUE'),
                       nz.addTo(nz.newSection(vec.name), sv_fig1, sv_fig_2))
    table1 = nz.addTo(table1, result1, row=pos[0], column=pos[1])
    return table1
    
def pathway_mutation_section(cancer, gene_sets, cutoff=.25, folder='/'):
    #Format data for report
    path = cancer.data_path + folder + '/'
    pathway_table_file = path + 'pathway_table.csv'
    pathway_table = format_pathway_table(cancer, gene_sets)    
    if 'survival' in pathway_table:
        pathway_table = pathway_table.sort(columns='survival')
    pathway_table.to_csv(pathway_table_file)
    keepers = cancer.q_pathways[(cancer.q_pathways < .25).sum(1) > 0].index
    pathway_table = pathway_table.ix[keepers]
    if 'survival' in pathway_table:
        pathway_table = pathway_table.sort(columns='survival')
    pathway_table = pathway_table.head(20)
    pathway_table_r = com.convert_to_r_dataframe(pathway_table.replace(nan, 1.23)) #@UndefinedVariable
    if len(pathway_table) == 0:
        return nz.addTo(nz.newSubSection('Pathway Mutations'), nz.newParagraph(''))
    
    #Overview
    tableCaption1 = ('Association of pathway level mutations with patient' + 
                     'clinical features.')
    table1 = nz.newTable(pathway_table_r, tableCaption1, file=pathway_table_file, 
                             significantDigits=2);                      
   
    #Fill in the details
    pathway_pos = dict((p,i+1) for i,p in enumerate(pathway_table.index))
    col_pos = dict((c,i+1) for i,c in enumerate(pathway_table.columns))
    
    #age violin plots
    if 'age' in pathway_table:
        for g,val in pathway_table['age'].iteritems():
            num_patients = (match_series(cancer.meta_matrix.ix[g], cancer.clinical.age)[0] > 0).sum()
            if val < cutoff and num_patients > 2:
                table1 = add_violin_plot(cancer.meta_matrix.ix[g], cancer, table1, 
                                         (pathway_pos[g], col_pos['age']),
                                         path + FIG_EXT)        
    
    #survival curves
    if 'survival' in pathway_table:
        for g,val in pathway_table['survival'].iteritems():
            if val < cutoff:
                table1 = add_survival_curve_pathway(cancer.meta_matrix.ix[g], cancer, table1, 
                            (pathway_pos[g], col_pos['survival']), path + FIG_EXT) 
                
    section = nz.addTo(nz.newSubSection('Pathway Mutations'), table1)
    return section

def format_pathway_table_exp(cancer, gene_sets):
    pathways = Series(dict((i,i) for i in cancer.q_pathways.index), name='pathway')
    n_genes = Series(dict((p, len(gene_sets[p])) for p in cancer.pc.index),
                     name='n genes')
    pathway_table = cancer.q_pathways.join(n_genes)
    pathway_table = pathway_table.join(pathways)
    pathway_table = pathway_table.ix[:, range(pathway_table.shape[1]-1, -1, -1)]
    return pathway_table

    
def pathway_mutation_section_exp(cancer, gene_sets, cutoff=.25, folder=''):
    #Format data for report
    path = cancer.data_path + folder + '/'
    pathway_table_file = path + 'pathway_table.csv'
    pathway_table = format_pathway_table_exp(cancer, gene_sets) 
    if 'survival' in pathway_table:
        pathway_table.sort(columns='survival')
    pathway_table.to_csv(pathway_table_file)
    keepers = cancer.q_pathways[(cancer.q_pathways < .25).sum(1) > 0].index
    pathway_table = pathway_table.ix[keepers]
    if 'survival' in pathway_table:
        pathway_table = pathway_table.sort(columns='survival')
    pathway_table = pathway_table.head(20)
    pathway_table_r = com.convert_to_r_dataframe(pathway_table.replace(nan, 1.23)) #@UndefinedVariable
    if len(pathway_table) == 0:
        return nz.addTo(nz.newSubSection('Expressed Pathways'), nz.newParagraph(''))
    
    #Overview
    tableCaption1 = ('Association of pathway level expression patterns with patient' + 
                     'clinical features.')
    table1 = nz.newTable(pathway_table_r, tableCaption1, file=pathway_table_file, 
                             significantDigits=2);                      
   
    #Fill in the details
    pathway_pos = dict((p,i) for i,p in enumerate(pathway_table.index))
    col_pos = dict((c,i) for i,c in enumerate(pathway_table.columns))
    
    #age scatter plots
    for p in (pathway_table['age'][pathway_table['age'] < cutoff]).index:
        fig_file = path + FIG_EXT + p + '_age.png'
        draw_pathway_age_scatter(p, cancer, fig_file)
        age_fig1 = nz.newFigure(fig_file, 'Age of patients with or without' +
                                           'mutation to pathway.')
        result1 = nz.addTo(nz.newResult('', isSignificant='TRUE'),
                           nz.addTo(nz.newSection(p), age_fig1))
        table1 = nz.addTo(table1, result1, row=pathway_pos[p]+1, 
                          column=col_pos['age']+1)
        
    #survival curves
    for p in (pathway_table['survival'][pathway_table['survival'] < cutoff]).index:
        fig_file = path + FIG_EXT + p + '_survival.png'
        data_frame = cancer.data_matrix.ix[gene_sets[p]].dropna()
        U,S,vH = frame_svd(((data_frame.T - data_frame.mean(1)) / data_frame.std(1)).T)
        
        strat = (vH[0] > vH[0].std()).astype(int) - (vH[0] < -vH[0].std()) + 1
        draw_survival_curves(cancer.clinical, Series(strat, name='pc'), 
                             labels=['low','mid','high'], filename=fig_file)
        sv_fig1 = nz.newFigure(fig_file, 'Survival of patients with ' + 
                                          'varying levels of pathway expression.')
        fig_file2 = path + FIG_EXT + p + '.svg'
        draw_pathway_eig_bar(U, fig_file2)
        sv_fig_2 = nz.newFigure(fig_file2, 'Loading for first eigen-patient.')
        result1 = nz.addTo(nz.newResult('', isSignificant='TRUE'),
                           nz.addTo(nz.newSection(p), sv_fig1, sv_fig_2))
        table1 = nz.addTo(table1, result1, row=pathway_pos[p]+1, 
                          column=col_pos['survival']+1)
        
    section = nz.addTo(nz.newSubSection('Pathway Mutations'), table1)
    return section

        
def create_clinical_report(cancer, next_cancer, prev_cancer, gene_sets, 
                           folder='Pathway_Mutations'):
    if not os.path.isdir(cancer.data_path + folder + '/' + FIG_EXT):
        os.makedirs(cancer.data_path + folder + '/' + FIG_EXT)
    report = nz.newReport('Report for ' + cancer.cancer)
    report = generic_header(report, cancer, next_cancer, prev_cancer, folder)
    if cancer.data_type not in ['expression', 'methylation']:
        hit_matrix = (cancer.hit_matrix if cancer.data_type == 'mutation' else 
              cancer.lesion_matrix)
        report = nz.addToResults(report, single_gene_section(cancer, hit_matrix, 
                                                             folder=folder))
        report = nz.addToResults(report, pathway_mutation_section(cancer, gene_sets, 
                                                                  folder=folder))
    else:
        report = nz.addToResults(report, pathway_mutation_section_exp(cancer, gene_sets, 
                                                                  folder=folder))
    nz.writeReport(report, filename=cancer.data_path + folder + '/index')