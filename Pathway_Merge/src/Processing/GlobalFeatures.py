'''
Created on Dec 17, 2012

@author: agross
'''

from pandas import DataFrame, read_table, concat
from Processing.Helpers import frame_svd
import Data.Firehose as FH
import Data.Intermediate as IM


def get_rates(data_path, cancer):
    mutsig_path = data_path + '/'.join(['analyses', cancer, 
                                    'MutSigNozzleReport2', ''])
    mutsig_path = mutsig_path + cancer + '-TP.'
    try:
        mutation_rates = read_table(mutsig_path + 
                                    'patients.counts_and_rates.txt')
    except:
        return DataFrame([])
    mutation_rates = mutation_rates.set_index('name')
    mutation_rates = mutation_rates.rename(index=lambda s: 
                                           '-'.join(['TCGA'] + s.split('-')[1:3]))
    mutation_rates = mutation_rates[['rate_dbsnp', 'rate_sil', 'rate_non']]
    
    maf = read_table(mutsig_path + 'final_analysis_set.maf')
    maf = maf.dropna(how='all', axis=[0,1])
    maf.Tumor_Sample_Barcode = maf.Tumor_Sample_Barcode.map(lambda s: s[:12])
    maf = maf.set_index(['Hugo_Symbol','Tumor_Sample_Barcode'])
    non_silent = maf[maf.is_silent == 0]
    df = DataFrame({pat: df.categ.value_counts() for pat, df in 
                    non_silent.groupby(level=1)}).T
    
    labels = read_table(mutsig_path + 'mutation_rates.txt')
    df = df.fillna(0)
    df = df.rename(columns=lambda s: labels.category[s-1])
    pct = (df.T / df.sum(1)).T    
    mutation_rates = mutation_rates.join(pct)
    mutation_rates = mutation_rates.rename(columns=lambda s: s.replace('/','_'))
    return mutation_rates


def get_cna_rates(cancer, data_path):
    gistic_ext = data_path + '/'.join(['analyses', cancer, 
                                       'CopyNumber_Gistic2', ''])
    try:
        gistic = read_table(gistic_ext + 'all_thresholded.by_genes.txt', 
                            index_col=0)
    except:
        return DataFrame([])
    gistic = gistic.rename(columns=lambda s: s[:12])
    gistic = gistic.ix[:,2:]
    
    amp_gene_all = (gistic >= 1).astype(int).sum()
    amp_gene_high = (gistic == 2).astype(int).sum()
    del_gene_all = (gistic <= -1).astype(int).sum()
    del_gene_homo = (gistic <= -2).astype(int).sum()
    
    gistic_lesions = read_table(gistic_ext + 'all_lesions.conf_99.txt', 
                                index_col=1)
    gistic_lesions = gistic_lesions.rename(columns=lambda s: s[:12])
    
    calls = gistic_lesions['Unique Name'].apply(lambda s: (s.find('CN') < 0) 
                                                and (s.find('Amp') == 0))
    lesions = gistic_lesions[calls].select(lambda s: 'TCGA' in s, 1) 
    amp_lesion_all = (lesions >= 1).astype(int).sum()
    amp_lesion_high = (lesions == 2).astype(int).sum()
    
    calls = gistic_lesions['Unique Name'].apply(lambda s: (s.find('CN') < 0) 
                                                and (s.find('Del') == 0))
    lesions = gistic_lesions[calls].select(lambda s: 'TCGA' in s, 1)
    del_lesion_all = (lesions >= 1).astype(int).sum()
    del_lesion_homo = (lesions == 2).astype(int).sum()
    
    arm_cn = read_table(gistic_ext + 'broad_values_by_arm.txt', index_col=0)
    arm_cn = arm_cn.rename(columns=lambda s: s[:12])
    chromosomal_instability = arm_cn.abs().mean()
    cna_df = {'gene_amp': amp_gene_all, 'gene_amp_high': amp_gene_high, 
              'gene_del': del_gene_all, 'gene_del_homo': del_gene_homo, 
              'lesion_amp': amp_lesion_all, 'lesion_amp_high': amp_lesion_high, 
              'lesion_del': del_lesion_all, 'lesion_del_homo': del_lesion_homo,
              'chrom_instability': chromosomal_instability}
    cna_df = DataFrame(cna_df)
    return cna_df

def get_global_vars(cancer, data_path, patients=None, filtered_patients=None):
    def patient_filter(df):
        if patients is not None:
            return df[[p for p in df.columns if p in patients]]
        elif filtered_patients is not None:
            return df[[p for p in df.columns if p not in filtered_patients]]
        else:
            return df
    try:
        data_matrix = FH.read_rnaSeq(cancer, data_path, average_on_genes=True)
        data_matrix = patient_filter(data_matrix)
        U, S, vH = frame_svd(data_matrix)
        exp_pc = DataFrame({'pc1': vH[0], 'pc2': vH[1]})
    except:
        exp_pc = DataFrame()
        
    try:
        data_matrix = IM.read_methylation(cancer, data_path)
        data_matrix = patient_filter(data_matrix)
        U, S, vH = frame_svd(data_matrix)
        meth_pc = DataFrame({'pc1': vH[0], 'pc2': vH[1]})
    except:
        meth_pc = DataFrame()
        
    try:
        meth_age, amar = 'FAIL','FAIL'
        #meth_age, amar = get_age_signal(data_path, cancer) 
        meth_pc = meth_pc.join(meth_age).join(amar)
        print 'Should probably check this out'
    except:
        pass
    
    cna_rates = get_cna_rates(cancer, data_path)
    cna_rates = patient_filter(cna_rates.T).T
    mutation_rates = get_rates(data_path, cancer)
    mutation_rates = patient_filter(mutation_rates.T).T
    
    gv = concat([exp_pc, meth_pc, cna_rates, mutation_rates], 
                 keys=['mRNASeq','methylation', 'cna', 'mutation'], axis=1)
    gv = gv.dropna(how='all', axis=1)
    return gv