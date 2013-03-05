'''
Created on Sep 4, 2012

Set of functions to read in data that has been formated in the 
BROAD GDAC Firehose data processing pipeline. 

Nothing in here should depend on any other modules.

@author: agross
'''
import numpy as np
import os as os
from pandas import DataFrame, Series, read_table, read_csv
import pandas as pd

from Processing.Tests import kruskal_p
from Processing.Helpers import bhCorrection

def read_clinical(data_path):
    clinical = read_table(data_path + 'From_Stddata/' + 
                                 'clin.merged.picked.txt', index_col=0, 
                                 header=0, skiprows=[1])
    clinical = clinical.rename(columns=str.upper).T
    age = clinical['yearstobirth'].dropna().astype(int)
    days = clinical[['daystolastfollowup','daystodeath']]
    days = days.apply(np.nanmax,1).dropna().astype(int)
    days = days[days > 0]
    censored = clinical.vitalstatus.dropna().astype(int)
    gender = clinical.gender
    check_feature = lambda f: clinical[f] if f in clinical else Series()
    therapy = check_feature('neoadjuvanttherapy')
    rrr = check_feature('radiations.radiation.regimenindication')
    clinical = DataFrame({'age': age, 'days': days, 'censored': censored, 
                         'gender': gender, 'therapy': therapy, 'radiation': rrr})
    clinical = clinical.replace('yes',True).replace('no',False)
    clinical = clinical.replace('male',True).replace('female',False)
    return clinical

def get_mutation_matrix_o(cancer, data_path, q_cutoff=.25, get_mutsig=False):
    f = data_path + '/'.join(['broad_analyses', cancer, 'MutSigRun2',''])
    if os.path.isdir(f):
        mutation_matrix = read_table(f + cancer + '.per_gene.mutation_counts.txt', 
                                     index_col=0)
        mutation_matrix = mutation_matrix.select(lambda s: s.count('_') >= 2, axis=1)
        adjust_barcodes = lambda s: '-'.join(['TCGA'] + s.split('_')[1:3])
        mutation_matrix = mutation_matrix.rename(columns=adjust_barcodes)
    else:
        mutation_matrix = data_path + '/'.join(['ucsd_processing', 'mutation', 
                                                'mutation_matrix.csv'])
        mutation_matrix = read_csv(mutation_matrix, index_col=0)
        mutation_matrix = mutation_matrix.apply(np.nan_to_num)
        sig_genes = np.array(mutation_matrix[mutation_matrix.sum(1) > 3].index)
    if get_mutsig:
        mutsig = read_table(f + cancer +  '.sig_genes.txt', index_col=1)
        if mutsig.q.dtype != type(.1):
            mutsig.q = mutsig.q.map(lambda s: float(s.replace('<','')))
        sig_genes = np.array(mutsig[mutsig.q < .25].index)
        return mutation_matrix, mutsig
    else:
        sig_genes = np.array(mutation_matrix[mutation_matrix.sum(1) > 3].index)
        return mutation_matrix, sig_genes
 
def get_mutation_matrix(data_path, cancer):
    maf = read_table(data_path + '/'.join(['analyses', cancer, 'MutSigNozzleReport2', 
                                             cancer + '-TP.final_analysis_set.maf']))
    maf = maf.dropna(how='all', axis=[0,1])
    maf.Tumor_Sample_Barcode = maf.Tumor_Sample_Barcode.map(lambda s: s[:12])
    maf = maf.set_index(['Hugo_Symbol','Tumor_Sample_Barcode'])
    non_silent = maf[maf.is_silent == 0]
    non_silent['counter'] = 1
    hit_matrix = non_silent.counter.groupby(level=[0,1]).sum().unstack()
    return hit_matrix   
    
def get_gistic_gene_matrix(data_path, cancer):
    gistic_ext = data_path + '/'.join(['analyses', cancer, 'CopyNumber_Gistic2', ''])
    gistic = read_table(gistic_ext + 'all_thresholded.by_genes.txt', index_col=[2,1,0])
    gistic = gistic.rename(columns=lambda s: s[:12]).astype(float)
    return gistic

def get_cna_matrix(cancer, data_path, cna_type='deletion'):
    gistic_ext = data_path + '/'.join(['broad_analyses', cancer, 'CopyNumber_Gistic2', ''])
    gistic = read_table(gistic_ext + 'all_thresholded.by_genes.txt', index_col=0)
    gistic = gistic.rename(columns=lambda s: s[:12])
    gistic = gistic.ix[:,2:]
    cna_val = {'deletion': -2, 'amplification': 2}[cna_type]
    cnas = (gistic == cna_val).astype(int)
    cnas = cnas[cnas.sum(1) > 0]
    
    gistic_lesions = read_table(gistic_ext + 'all_lesions.conf_99.txt', 
                                index_col=1)
    gistic_lesions = gistic_lesions.rename(columns=lambda s: s[:12])
    call_label = {'amplification': 'Amp', 'deletion': 'Del'}[cna_type]
    calls = gistic_lesions['Unique Name'].apply(lambda s: (s.find('CN') < 0) 
                                                and (s.find(call_label) ==0))
    lesions = gistic_lesions[calls].select(lambda s: 'TCGA' in s, 1) 
    lesions = (lesions == cna_val).astype(int)
    return cnas, lesions

def read_rppa(data_path, cancer):
    ext = '/'.join(['stddata', cancer,  'RPPA_AnnotateWithGene', 
                    cancer + '.rppa.txt'])
    rppa = read_table(data_path +  ext, index_col=0)
    rppa = rppa.select(lambda s: s.split('-')[3].startswith('01'), 1)
    rppa = rppa.sort_index(axis=1).groupby(lambda s: s[:12], 1).first() #get rid of normals
    rppa['protien'] = rppa.index.map(lambda s: s.split('|')[0])
    rppa['antibody'] = rppa.index.map(lambda s: s.split('|')[1])
    rppa = rppa.set_index(['protien','antibody'])
    return rppa

def read_rnaSeq(cancer, data_path, patients=None, average_on_genes=False):
    stddata_path = data_path + 'stddata/' + cancer
    data_types = filter(lambda f: f.startswith('rnaseq'), os.listdir(stddata_path))
    if 'rnaseqv2' in data_types:
        path = [f[0] for f in list(os.walk(stddata_path + '/rnaseqv2')) if 
                'RSEM_genes_normalized/data' in f[0]][0]   
    else:
        path = [f[0] for f in list(os.walk(stddata_path + '/rnaseqv')) if 
                'gene_expression/data' in f[0]][0] 
    rnaSeq = read_table(path + '/data.txt',index_col=0, skiprows=[1])
    rnaSeq = np.log2(rnaSeq).clip_lower(0)
    rnaSeq = rnaSeq.sort_index(axis=1).groupby(lambda s: s[:12]).first() #get rid of normals
    rnaSeq = rnaSeq.select(lambda s: s.split('-')[3].startswith('01'), 1)
    rnaSeq = rnaSeq.rename(columns=lambda s: s[:12])
    if patients is not None:
        rnaSeq  = rnaSeq.ix[:, patients]
    rnaSeq = rnaSeq.dropna(thresh=100)
    if average_on_genes:
        rnaSeq = rnaSeq.groupby(by=lambda n: n.split('|')[0]).mean()
    return rnaSeq

def read_methylation(cancer, data_path, patients=None):
    processed_data_path = data_path + '/'.join(['ucsd_processing', cancer,''])
    data_types = filter(lambda f: f[:11] == 'methylation', os.listdir(processed_data_path))
    data_type = sorted([d for d in data_types if os.path.isfile(processed_data_path + d + 
                                                      '/meta_probes.csv')])[-1]
    meth = read_csv(processed_data_path +  data_type + '/meta_probes.csv', index_col=0)
    meth = meth.sort_index(axis=1).groupby(lambda s: s[:12]).first() #get rid of normals
    meth = meth.select(lambda s: s.split('-')[3].startswith('01'), 1)
    meth = meth.rename(columns=lambda s: s[:12])
    
    if patients is not None:
        meth  = meth.ix[:, patients]
    meth = meth.dropna(thresh=100)
    meth = meth.astype(np.float)
    return meth

def read_mrna(cancer, data_path, patients=None, num_genes='All'):
    stddata_path = data_path + '/'.join(['stddata', cancer,''])
    chips = filter(lambda f: 'transcriptome__agilent' in f, os.listdir(stddata_path))
    data_type = sorted(chips)[-1]
    mrna = read_table(stddata_path +  data_type + 
                      '/unc_lowess_normalization_gene_level_data.txt',
                      index_col=0, skiprows=[1], na_values=['null'])
    mrna = mrna.rename(columns=lambda s: s[:12])
    if patients is not None:
        mrna  = mrna.ix[:, patients]
    mrna = mrna.astype(np.float)
    if num_genes is not 'All':
        variable = np.sort(mrna.mad(axis=1)).tail(num_genes).index
        mrna = mrna.ix[variable]
    return mrna

def read_mirna(cancer, data_path):
    stddata_path = data_path + '/'.join(['stddata', cancer,''])
    data_type = 'mirnaseq__illuminahiseq_mirnaseq__bcgsc_ca__Level_3'
    mirna_ga = read_table(stddata_path +  data_type.replace('hiseq', 'ga') + 
                          '/miR_gene_expression_data.txt', 
                          index_col=0, na_values=['null'])
    try: 
        mirna_hiseq = read_table(stddata_path +  data_type + 
                                 '/miR_gene_expression_data.txt',
                                 index_col=0, na_values=['null'])
        mirna = mirna_hiseq.join(mirna_ga)
    except IOError:
        mirna = mirna_ga
    mirna = mirna.ix[1:,mirna.ix[0] == 'reads_per_million_miRNA_mapped']
    mirna = mirna.rename(columns=lambda s: s[:12])
    #mirna = mirna.ix[:, cancer.patients]
    mirna = np.log2(mirna.astype(np.float))
    #variable = np.sort(mirna.mad(axis=1)).tail(150).index
    #mirna = mirna.ix[variable]
    return mirna

def read_cnmf(cancer, file_name, data_path):
    nmf_file = data_path + '/'.join(['broad_analyses', cancer, file_name, 
                                     'cnmf.normalized.gct'])
    if not os.path.isfile(nmf_file):
        print 'Missing ' + file_name + ' file.'
        return None
    data_type = file_name.split('_')[0]
    data = read_table(nmf_file, index_col=0, skiprows=[0,1], sep='\t')
    data = data.rename(index=lambda s: data_type + '_' + s, 
                       columns=lambda s: s.replace('.','-')[:12])
    data = data.select(lambda s: 'TCGA' in s, axis=1)
    return data

def get_all_nmf(cancer):
    nmf = [f for f in os.listdir(cancer.data_path) if 'CNMF' in f]
    variable_signals = map(read_cnmf, nmf, [cancer.data_path]*len(nmf))
    variable_signals = reduce(DataFrame.append, variable_signals)
    variable_signals = variable_signals.ix[:,cancer.patients]
    return variable_signals

def get_gistic_lesions(cancer_name, data_path):
    gistic_ext = data_path + '/'.join(['analyses', cancer_name, 
                                       'CopyNumber_Gistic2', ''])
    gistic = read_table(gistic_ext + 'all_lesions.conf_99.txt', index_col=[0,1])
    lesions = gistic.select(lambda s: 'TCGA' in s, axis=1) 
    lesions = lesions.select(lambda s: 'values' not in s[0], axis=0) 
    lesions = lesions.rename(columns=lambda s: s[:12])
    from_tuples = pd.MultiIndex.from_tuples
    lesions.index = from_tuples([(s[0].split(' ')[0], s[1].strip(), 'Lesion') 
                                 for s in lesions.index])
    lesions = lesions.groupby(level=[0,1,2]).first()
    lesions.T['Deletion'] = (lesions.T['Deletion']*-1).replace(-0,0)
    return lesions

def get_gistic_genes(cancer_name, data_path, filter_with_rna=True, 
                     collapse_on_bands=True, min_patients=5):
    gistic_ext = data_path + '/'.join(['analyses', cancer_name, 
                                       'CopyNumber_Gistic2', ''])
    gistic = read_table(gistic_ext + 'all_thresholded.by_genes.txt', 
                        index_col=[2,1,0])
    gistic = gistic.rename(columns=lambda s: s[:12]).astype(float)

    deletion = gistic[(gistic == -2).sum(1) > min_patients]
    amp = gistic[(gistic == 2).sum(1) > min_patients]
    ft = pd.MultiIndex.from_tuples #rediculously long pandas names
    deletion.index = ft([('Deletion', s[0], s[2]) for s in deletion.index])
    amp.index = ft([('Amplification', s[0], s[2]) for s in amp.index])
    
    if filter_with_rna:
        rna = read_rnaSeq(cancer_name, data_path, average_on_genes=True)
        def rna_filter(df, val):
            p_vals = Series({g: kruskal_p(df.ix[g] == val, rna.ix[g[-1]]) 
                             for g in df.index if g[-1] in rna.index})
            q_vals = bhCorrection(p_vals)
            return df.ix[q_vals[q_vals < .1].index]
        deletion = rna_filter(deletion, -2)
        amp = rna_filter(amp, 2)
   
    cna_genes = amp.append(deletion)
    if collapse_on_bands == False:
        return cna_genes
    
    cna_genes = DataFrame({(a[0], a[1], tuple(b.index.get_level_values(2))): 
                           b.mean().round() for a,b in 
                           cna_genes.groupby(level=[0,1])}).T
    cna_genes.index = pd.MultiIndex.from_tuples(cna_genes.index)
    return cna_genes

def get_gistic(cancer_name, data_path, filter_with_rna=True, 
               collapse_on_bands=True, min_patients=5):
    lesions = get_gistic_lesions(cancer_name, data_path)
    cna_genes = get_gistic_genes(cancer_name, data_path, filter_with_rna, 
                                 collapse_on_bands, min_patients)
    cna = cna_genes.append(lesions)
    return cna
    