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

def get_mutation_matrix(cancer, data_path, q_cutoff=.25, get_mutsig=False):
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
    data_type = 'RPPA_AnnotateWithGene'
    rppa = read_table(data_path +  'stddata/RPPA_AnnotateWithGene/' + 
                      cancer.cancer + '.rppa.txt', index_col=0)
    rppa = rppa.rename(columns=lambda s: s[:12]) #short barcode format
    rppa = np.log2(rppa.astype(np.float))
    rppa = rppa.ix[:,cancer.patients]
    return rppa

def read_rnaSeq(cancer, data_path, patients=None):
    stddata_path = data_path + '/'.join(['stddata', cancer,''])
    data_types = filter(lambda f: f[:6] == 'rnaseq', os.listdir(stddata_path))
    data_type = sorted(data_types)[-1]
    if data_type.split('__')[0] == 'rnaseqv2':
        rnaSeq = read_table(stddata_path +  data_type + 
                            '/RSEM_genes_normalized_data.txt',
                            index_col=0, skiprows=[1])
    else:
        rnaSeq = read_table(stddata_path +  data_type + 
                            '/gene_expression_data.txt', index_col=0)
        rnaSeq = rnaSeq.ix[1:,rnaSeq.ix[0] == 'RPKM']
    
    rnaSeq = rnaSeq.rename(columns=lambda s: s[:12])
    if patients is not None:
        rnaSeq  = rnaSeq.ix[:, patients]
    rnaSeq = rnaSeq.dropna(thresh=100)
    #rnaSeq = np.log2(rnaSeq.astype(np.float))
    return rnaSeq

def read_methylation(cancer, data_path, patients=None):
    processed_data_path = data_path + '/'.join(['ucsd_processing', cancer,''])
    data_types = filter(lambda f: f[:11] == 'methylation', os.listdir(processed_data_path))
    data_type = sorted([d for d in data_types if os.path.isfile(processed_data_path + d + 
                                                      '/averaged_on_genes.csv')])[-1]
    meth = read_csv(processed_data_path +  data_type + '/averaged_on_genes.csv', index_col=0)
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
    variable = np.sort(mirna.mad(axis=1)).tail(150).index
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
    