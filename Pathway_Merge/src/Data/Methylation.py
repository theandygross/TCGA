'''
Created on Oct 29, 2012

@author: agross
'''
import os as os
from subprocess import call, check_output
from pandas import read_table

METH_FOLDER = 'methylation/humanmethylation450/jhu_usc_edu__Level_3/'
METH_FILE = 'within_bioassay_data_set_function_data.txt'

def pull_out_beta_values(folder, probes='All', outfile='beta_values.txt'):
    
    curdir = os.getcwd()
    os.chdir(folder)
    
    if probes != 'All':
        target_file = 'picked.txt'
        query = '\\|'.join(['Hybridization REF'] + list(probes))
        call(['grep', query, METH_FILE], 
             stdout=open(target_file,'wb'))
    else:
        target_file = METH_FILE
    cols = check_output(['awk', '-F\t', '{print NF; exit}', target_file])
    cols = int(cols.strip())
    cols = ','.join(map(str, range(2, int(cols)+1,4)))
    
    call(['cut', '-f' + cols, target_file], stdout=open('tmp','wb'))
    call(['cut', '-f1,3', target_file], stdout=open('idx','wb'))
    call(['paste', 'idx', 'tmp', '-d\t'], stdout=open(outfile,'wb'))
    call(['rm', 'tmp', 'idx'])
    
    os.chdir(curdir)
    
def average_beta_values_on_genes(folder):
    '''
    Memory intensive, I'm reading the ~3-4 GB file into memory to build a 
    pivot table.  Just do this on sauterne.
    ''' 
    table = read_table(folder + 'beta_values.txt', skiprows=[1], index_col=[0])
    table = table.rename(columns=lambda s: s[:12] if s!=table.columns[0] 
                         else 'symbol')
    table = table.groupby(by='symbol').mean()
    table.to_csv(folder + 'averaged_on_genes.csv')
    
     
def run_all_cancers(data_path, probeset='All', recalc=False):     
    for cancer in os.listdir(data_path + 'stddata/'):
        outpath = data_path + '/'.join(['ucsd_processing', cancer, 
                                        'methylation450', ''])
        if not os.path.isdir(outpath):
            os.makedirs(outpath)
        outfile = outpath + ('beta_values.txt' if probeset == 'All' else 
                             'beta_values_picked.txt')
        folder = data_path + '/'.join(['stddata', cancer, METH_FOLDER])
        if os.path.isfile(folder + outfile) and not recalc:
            print 'Using precalculated data for ' + cancer + '.'
            continue
        elif os.path.isdir(folder):
            print cancer
            try:
                pull_out_beta_values(folder, probeset, outfile)
                if probeset == 'All':
                    average_beta_values_on_genes(outpath)
            except:
                print cancer + ' fail.'
        else:
            print 'No data for ' + cancer +'.'
            
            
    