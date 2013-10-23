'''
Created on Sep 4, 2012

Set of functions to read in data that has been formated in the 
BROAD GDAC Firehose data processing pipeline. 

Nothing in here should depend on any other modules.

I am relying heavily on pandas for this project, so the main 
goal of these functions is get data from the Firehose's tables
into Pandas data-structures that I can work with. 

There is a little bit of pre-processing that is done to get the 
files from Firehose into the local file-system in a reasonably 
organized hierarchy, that is done at the time of data download 
and in the run's initialization. These dependencies should 
eventually be circumvented.

@author: agross
'''
import os as os
import numpy as np
import pandas as pd

 
def fix_barcode_columns(df, patients=None, tissue_code='All', get_batch=False):
    '''
    Takes TCGA barcode and reformats it into a MultiIndex if all tissue_codes 
    are desired, or just pulls the correct tissue codes and filteres the 
    DataFrame.

    df: pandas DataFrame
    patients: patient list to filter on
    tissue_code: ['01','11','All']  #if all returns MultiIndex

    '''
    if get_batch is False:
        df.columns = pd.MultiIndex.from_tuples([(i[:12], i[13:15]) for i 
                                                in df.columns])
    else:
        df.columns = pd.MultiIndex.from_tuples([(i[:12], i[13:15], i[21:24]) for i 
                                                in df.columns])
    if patients is not None:
        df  = df.ix[:, patients]
    if tissue_code != 'All':
        try:
            df = df.T.xs(tissue_code, level=1).T #pandas bug
            df = df.groupby(axis=1, level=0).first()
        except KeyError: #try different cross-seciton
            new_code = pd.value_counts(df.columns.get_level_values(1)).idxmax()
            df = df.T.xs(new_code, level=1).T #pandas bug
            df = df.groupby(axis=1, level=0).first()
            
    else:
        df = df.groupby(axis=1, level=[0,1]).first()
    return df

def get_dataset_path(data_path, cancer, data_type, ext):
    '''
    This is a helper to get paths to a particular data-set.  
    In processing the data, we develop relatively complicated file hierarchies
    to not have a ton of folders at the top level and make it easier to get 
    to files manually.  This makes it a little tough to track down files 
    automatically so this function fills that role. 
    
    data_type: the top-level data-type of the file (i.e. rnaseqv2, 
               transcriptome,...)
    ext: the file you are looking for (i.e. RSEM_genes_normalized, 
         junction_quantification, ...) 
    '''
    stddata_path = data_path + 'stddata/' + cancer
    if not os.path.isdir(stddata_path):
        return
    data_types = filter(lambda f: f.startswith(data_type), 
                        os.listdir(stddata_path))
    if data_type in data_types: #get the paths
        paths = [f[0] for f in list(os.walk(stddata_path + '/' + data_type)) if 
                (ext + '/data') in f[0]]
    else:
        return
    f = [path + '/'+ f for path in paths 
                       for f in os.listdir(path) 
                       if 'data' in f] #pull the data file
    return f


def get_mutation_matrix(data_path, cancer, tissue_code='01'):
    '''
    Get gene by patient mutation matrix. 
    Here I filter by the is_silent column in the MAF file, 
    so I am returning only non-silent mutations.  
    '''
    path = '{}/analyses/{}/MutSigNozzleReport2/'.format(data_path, cancer)
    maf = pd.read_table(path + cancer + '-TP.final_analysis_set.maf')
    maf = maf.dropna(how='all', axis=[0,1])
    maf = maf.set_index(['Hugo_Symbol','Tumor_Sample_Barcode'])
    non_silent = maf[maf.is_silent == 0]
    non_silent['counter'] = 1
    hit_matrix = non_silent.counter.groupby(level=[0,1]).sum().unstack()
    hit_matrix = fix_barcode_columns(hit_matrix, tissue_code=tissue_code)
    return hit_matrix  

def get_submaf(data_path, cancer, genes='All', fields='basic'):
    '''
    Pull a sub-section of the MAF file for analysis.  

    genes: list of genes for which to return data
    fields: ['basic', 'all']: if basic, returns reduced version of MAF
    '''
    path = '{}/analyses/{}/MutSigNozzleReport2/'.format(data_path, cancer)
    maf = pd.read_table(path + cancer + '-TP.final_analysis_set.maf')
    maf = maf.dropna(how='all', axis=[0,1])
    maf['Tissue_Type'] = maf.Tumor_Sample_Barcode.map(lambda s: s[13:15])
    maf.Tumor_Sample_Barcode = maf.Tumor_Sample_Barcode.map(lambda s: s[:12])
    if genes != 'All':
        maf = maf[maf.Hugo_Symbol.isin(genes)]
    def get_allele(s):
        alleles = [s['Tumor_Seq_Allele1'], s['Tumor_Seq_Allele2']]
        return [a for a in alleles if a != s['Reference_Allele']][0]
                            
    maf['Alt_Allele'] = maf.apply(get_allele, 1)
    if fields == 'basic':
        maf = maf[['Hugo_Symbol','NCBI_Build', 'Chromosome', 'Start_position', 
                   'End_position', 'Strand', 'Reference_Allele', 
                   'Alt_Allele', 'Tumor_Sample_Barcode']]
    maf = maf.set_index('Hugo_Symbol', append=True)
    maf.index = maf.index.swaplevel(0,1)
    return maf 

def get_gistic_gene_matrix(data_path, cancer, tissue_code='01'):
    '''
    Reads in gene by patient copy-number alteration matrix.  
    Index is MultiIndex with ['Cytoband', 'Locus ID', 'Gene Symbol'] 
    on the levels.
    '''
    path = '{}/analyses/{}/CopyNumber_Gistic2/'.format(data_path, cancer)
    gistic = pd.read_table(path + 'all_thresholded.by_genes.txt', 
                           index_col=[2,1,0]) 
    gistic = fix_barcode_columns(gistic, tissue_code=tissue_code)
    return gistic

def get_gistic_arm_values(data_path, cancer, tissue_code='01'):
    '''
    Reads in arm by patient copy-number alteration matrix.  
    '''
    path = '{}/analyses/{}/CopyNumber_Gistic2/'.format(data_path, cancer)
    gistic = pd.read_table(path + 'broad_values_by_arm.txt', index_col=0) 
    gistic = fix_barcode_columns(gistic, tissue_code=tissue_code)
    return gistic

def get_gistic_lesions(data_path, cancer, patients=None, tissue_code='01'):
    '''
    Reads in lesion by patient CNA matrix.
    Returns thresholded calls as made by GISTIC2 in the Firehose pipeline. 
    '''
    path = '{}/analyses/{}/CopyNumber_Gistic2/'.format(data_path, cancer)
    gistic = pd.read_table(path + 'all_lesions.conf_99.txt', index_col=[0,1])
    lesions = gistic.select(lambda s: 'TCGA' in s, axis=1) 
    lesions = lesions.select(lambda s: 'values' not in s[0], axis=0) 
    from_tuples = pd.MultiIndex.from_tuples
    lesions.index = from_tuples([(s[0].split(' ')[0], s[1].strip(), 'Lesion') 
                                 for s in lesions.index])
    lesions = lesions.groupby(level=[0,1,2]).first()
    lesions.T['Deletion'] = (lesions.T['Deletion']*-1).replace(-0,0)
    lesions = fix_barcode_columns(lesions, patients, tissue_code)
    return lesions

def get_gistic(cancer_name, data_path, filter_with_rna=True, 
               collapse_on_bands=True, min_patients=5):
    lesions = get_gistic_lesions(cancer_name, data_path)
    return lesions

def read_rppa(data_path, cancer, patients=None, tissue_code='01'):
    '''
    Reads in antibody by patient reverse-phase protein array matrix. 
    Index is MultiIndex with ['protien','antibody'] on the levels. 
    '''
    path = '{}/stddata/{}/RPPA_AnnotateWithGene/'.format(data_path, cancer)
    rppa = pd.read_table(path +  cancer + '.rppa.txt', index_col=0)
    rppa['protien'] = rppa.index.map(lambda s: s.split('|')[0])
    rppa['antibody'] = rppa.index.map(lambda s: s.split('|')[1])
    rppa = rppa.set_index(['protien','antibody'])
    rppa = fix_barcode_columns(rppa, tissue_code=tissue_code)
    return rppa

def read_rnaSeq(data_path, cancer, patients=None, average_on_genes=True,
                tissue_code='01', get_batch=False):
    '''
    Reads in gene by patient rnaSeq mRNA expression matrix. 
    Data is log-transformed and a lower bound of -3 (1/8 read per million) 
    is set.
    '''
    files = get_dataset_path(data_path, cancer, 'rnaseqv2', 
                             'RSEM_genes_normalized')
    if files is None:
        return 
    
    rnaSeq = pd.concat([pd.read_table(f, index_col=0, skiprows=[1])
                        for f in files])
    rnaSeq = np.log2(rnaSeq).replace(-np.inf, -3.) #close enough to 0
    if average_on_genes:  #Pretty much all duplicates are unknown ('?')
        rnaSeq = rnaSeq.groupby(by=lambda n: n.split('|')[0]).mean()
    rnaSeq = fix_barcode_columns(rnaSeq, patients, tissue_code, get_batch)
    return rnaSeq

def read_rnaSeq_splice_junctions(data_path, cancer, patients=None,
                                 tissue_code='01'):
    '''
    Reads in gene by patient rnaSeq mRNA splice junction matrix. 
    Values are raw counts. 
    '''
    files = get_dataset_path(data_path, cancer, 'rnaseqv2', 
                         'junction_quantification')
    if files is None:
        return 
    rnaSeq = pd.concat([pd.read_table(f, index_col=0, skiprows=[1])
                        for f in files])
    rnaSeq = fix_barcode_columns(rnaSeq, patients, tissue_code)
    return rnaSeq

def read_mrna(data_path, cancer, patients=None, tissue_code='01'):
    '''
    Reads in gene by patient microarray gene expression.
    '''
    files = get_dataset_path(data_path, cancer, 'transcriptome', 
                             'unc_lowess_normalization_gene_level')
    if files is None:
        return 
    mrna = pd.concat([pd.read_table(f,index_col=0, skiprows=[1], 
                                    na_values=['null'])
                      for f in files])
    mrna = fix_barcode_columns(mrna, patients, tissue_code)
    return mrna

def read_miRNASeq(data_path, cancer, patients=None, tissue_code='01'):
    '''
    Reads in miRNA by patient miRNASeq matrix.
    Data is log-transformed and a lower bound of -3 (1/8 read per million) 
    is set.
    They often use both HiSeq and GA2 for this, so I'm merging the two,
    not entirely sure if that's Kosher. 
    '''
    stddata_path = data_path + 'stddata/' + cancer
    paths = [f[0] for f in list(os.walk(stddata_path + '/mirnaseq')) if 
                    'miR_gene_expression/data' in f[0]]
    data = []
    for path in paths:
        f = [f for f in os.listdir(path) if 'data' in f][0]
        mirna = pd.read_table(path + '/' + f, index_col=0, header=None)
        data.append(mirna)
    mirna = pd.concat(data, axis=1)
    mirna = mirna.T.set_index(['miRNA_ID', 'Hybridization REF'])
    mirna = mirna.sortlevel(level=0).ix['reads_per_million_miRNA_mapped']
    mirna = np.log2(mirna.astype(float)).replace(-np.inf, -3.) #close to 0
    mirna = mirna.T
    mirna = mirna.groupby(level=0, axis=1).last() #Use the HiSeq data over GA2
    mirna = fix_barcode_columns(mirna, patients, tissue_code)
    return mirna