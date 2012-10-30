'''
Created on Sep 6, 2012

@author: agross
'''
import os as os
import pandas.rpy.common as com
from pandas import Series, DataFrame
from numpy import nan
from rpy2.robjects.packages import importr


from Processing.Helpers import frame_svd, match_series
from Figures import violin_plot_pandas
from Figures import draw_survival_curves, draw_pathway_count_bar
from Figures import draw_pathway_age_scatter, draw_pathway_eig_bar
from Figures import histo_compare, mut_module_raster

nz = importr('Nozzle.R1')
bool_ = {True: 'TRUE', False: 'FALSE'}
FIG_EXT = 'clinical_figures/'
SORT_ORDER = ['survival','AMAR','age','gender','radiation','therapy']

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

def add_violin_plot(vec, cancer, table1, pos, fig_path):
    fig_file = fig_path + vec.name + '_age.png'
    fig = violin_plot_pandas(vec > 0, cancer.clinical.age.astype(float))
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

def add_eig_bar(pathway, cancer, table, pos, fig_path):
    fig_file = cancer.report_folder + '/' + FIG_EXT + pathway + + '.svg'
    if os.path.isfile(fig_file):
        data_frame = cancer.data_matrix.ix[cancer.gene_sets[pathway]].dropna()
        U,S,vH = frame_svd(((data_frame.T - data_frame.mean(1)) / data_frame.std(1)).T)
        draw_pathway_eig_bar(U, fig_file)
    sv_fig = nz.newFigure(fig_file, 'Loading for first eigen-patient.')
    result1 = nz.addTo(nz.newResult('', isSignificant='TRUE'),
                       nz.addTo(nz.newSection(pathway), sv_fig))
    table = nz.addTo(table, result1, row=pos[0], column=pos[1])
    
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

def add_histo_compare(hit_vec, response_vec, table, pos, fig_path, redraw=False):
    fig_file = (fig_path + str(hit_vec.name) + '_' + str(response_vec.name) + 
                '_histo_compare.png')
    if not os.path.isfile(fig_file) or redraw:
        fig = histo_compare(hit_vec, response_vec)
        fig.savefig(fig_file)
    histo_compare_fig = nz.newFigure(fig_file, 'Comparison of pathway level ' + 
                                               'distributions in patients with ' + 
                                               'and without mutation to pathway.')
    result1 = nz.addTo(nz.newResult('', isSignificant='TRUE'),
                       nz.addTo(nz.newSection(str(hit_vec.name) + ' vs. ' + 
                                str(response_vec.name)), 
                       histo_compare_fig))
    table = nz.addTo(table, result1, row=pos[0], column=pos[1])  
    return table
    
def add_mut_module_raster(cluster_num, mut, table, pos, fig_path, redraw=False):
    fig_file = fig_path + str(cluster_num) + '_mut_module_raster.png'
    if not os.path.isfile(fig_file) or redraw:
        fig = mut_module_raster(cluster_num, mut)
        fig.savefig(fig_file)
    fig = nz.newFigure(fig_file, 'Breakdown of patients covered by pathways in ' + 
                                 'the cluster.')
    p = mut.assignments[mut.assignments == cluster_num].index
    result1 = nz.addTo(nz.newResult('', isSignificant='FALSE'),
                       nz.addTo(nz.newSection('Mutation Pathway Cluster ' + 
                                              str(cluster_num)), 
                                nz.newParagraph(', '.join(p)), fig))
    table = nz.addTo(table, result1, row=pos[0], column=pos[1])  
    return table

def single_gene_section(cancer, hit_matrix, cutoff=.25):
    #Format data for report
    path = cancer.report_folder + '/'
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


    
def pathway_mutation_section(cancer, gene_sets, cutoff=.25):
    #Format data for report
    path = cancer.report_folder + '/'
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

    
def pathway_mutation_section_exp(cancer, gene_sets, cutoff=.25):
    #Format data for report
    path = cancer.report_folder + '/'
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



def format_interaction_table(path, table, mut_obj, exp_vec):
    interaction_table_file = path + '/interaction_table.csv'
    table.to_csv(interaction_table_file)
    table = table[['pathway','p_val','q_val','survival_p', 'age_p']].sort('p_val').head(20)
    table_r = com.convert_to_r_dataframe(table) #@UndefinedVariable
    tableCaption1 = "Association of pathway expression with mutated pathways."
    table1 = nz.newTable(table_r, tableCaption1, file=interaction_table_file, 
                         significantDigits=2);
    section = nz.addTo(nz.newSubSection('Association with Mutations'), table1)
    
    #Fill in the details
    gene_pos = dict((g,i+1) for i,g in enumerate(table.index))
    col_pos = dict((c,i+1) for i,c in enumerate(table.columns))
    
    for g,vec in table.iterrows():
        table1 = add_histo_compare(mut_obj.meta_matrix.ix[g], exp_vec, table1, 
                                   (gene_pos[g], col_pos['q_val']),
                                   path + '/' + FIG_EXT)
    return table1



def format_pathway_interaction_table(path, table, mut, exp):
    interaction_table_file = path + '/interaction_table.csv'
    table.to_csv(interaction_table_file)
    table_r = com.convert_to_r_dataframe(table) #@UndefinedVariable
    tableCaption1 = "Association of pathway expression with mutated pathways."
    table1 = nz.newTable(table_r, tableCaption1, file=interaction_table_file, 
                         significantDigits=2);
    #Fill in the details
    gene_pos = dict((g,i+1) for i,g in enumerate(table.index))
    col_pos = dict((c,i+1) for i,c in enumerate(table.columns))    
    for g,vec in table.iterrows():
        table1 = add_histo_compare(mut.clustered.ix[vec['Mut Group']], 
                                   exp.clustered.ix[vec['Exp Group']], 
                                   table1, (gene_pos[g], col_pos['q_val']),
                                   path + '/' + FIG_EXT)
        table1 = add_mut_module_raster(vec['Mut Group'], mut, table1, 
                                       (gene_pos[g], col_pos['Mut Group']), 
                                       path + '/' + FIG_EXT)
    return table1
        
def create_clinical_report(cancer, next_cancer, prev_cancer, gene_sets):
    if not os.path.isdir(cancer.report_folder + '/' + FIG_EXT):
        os.makedirs(cancer.report_folder + '/' + FIG_EXT)
    report = nz.newReport('Report for ' + cancer.cancer)
    report = generic_header(report, cancer, next_cancer, prev_cancer)
    if cancer.data_type not in ['expression', 'methylation']:
        hit_matrix = (cancer.hit_matrix if cancer.data_type != 'amplification' else 
              cancer.lesion_matrix)
        report = nz.addToResults(report, single_gene_section(cancer, hit_matrix))
        report = nz.addToResults(report, pathway_mutation_section(cancer, gene_sets))
    else:
        report = nz.addToResults(report, pathway_mutation_section_exp(cancer, gene_sets))
    nz.writeReport(report, filename=cancer.report_folder + '/index')