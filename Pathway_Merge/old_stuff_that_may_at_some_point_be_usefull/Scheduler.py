'''
Created on Feb 26, 2013

@author: agross
'''

import os as os
import pickle as pickle
import subprocess
import psutil as psutil #@UnresolvedImport

def parallel_run(process_list):
    processes = set()
    for process in process_list:
        fun, args = process
        processes.add(fun(*args))
        print '\t'.join([fun.func_name] + [str(a) for a in args if a != None])
        if ((len(processes) >= int(psutil.NUM_CPUS*.5)) or 
            (psutil.virtual_memory().percent > 50)):
            os.wait()
            s = [p for p in processes if p.poll() is not None]
            processes.difference_update(s)
            
def run_cancer(script, path, cancer, *args):
    r = subprocess.Popen(['nice', '-n', '5', script, path, cancer] +
                         list(args))
    return r

def initialize_cancer(path, cancer, patients=None, filtered_patients=None):
    if patients is not None:
        p = 'patients=\"' + str(list(patients)) + '\"'
        return run_cancer('InitializeCancer.py', path, cancer, p)
    elif filtered_patients is not None:
        p = 'filtered_patients=\"' + str(list(filtered_patients)) + '\"'
        return run_cancer('InitializeCancer.py', path, cancer, p)
    else:
        return run_cancer('InitializeCancer.py', path, cancer)

def initialize_data(path, cancer, data_type):
    if data_type == 'MAF':
        return run_cancer('InitializeMut.py', path, cancer)
    elif data_type.startswith('CN'):
        return run_cancer('InitializeCN.py', path, cancer, data_type)
    elif data_type == 'RPPA':
        return run_cancer('InitializeRPPA.py', path, cancer)
    else:
        return run_cancer('InitializeReal.py', path, cancer, data_type)

def run_survival(path, cancer, data_type, covariates, survival_test):
    covariates = ','.join(covariates)
    return run_cancer('run_surv.py', path, cancer, data_type, covariates,
                      survival_test)
    
def survival(run, data_type, covariates, survival_test, cancers='All',
             dry_run=False):
    if cancers == 'All':
        cancers = run.cancers
    queue = [(run_survival, (run.report_path, cancer, data_type, covariates,
                             survival_test)) 
             for cancer in cancers] 
    if dry_run:
        return queue
    parallel_run(queue)
    
def initialize(run, data_types, cancers='All', patients=None, 
               filtered_patients=None, init_cancer=False, dry_run=False):
    if cancers == 'All':
        cancers = run.cancers
    if patients is not None:
        patients = list(patients)
    if filtered_patients is not None:
        filtered_patients = list(filtered_patients)
    
    path = run.report_path
    
    if init_cancer:
        queue = [(initialize_cancer, (path, c, patients, filtered_patients))
                  for c in cancers]
        parallel_run(queue)
    
    queue = [(initialize_data, (path, c, data_type)) for c in cancers
              for data_type in data_types]
    if dry_run:
        return queue
    
    parallel_run(queue)
    
    
        
