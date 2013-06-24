import pickle as pickle
from Data.Containers import Cancer, Clinical
from Processing.GlobalFeatures import get_global_vars
from Processing.Helpers import make_path_dump



def initialize_cancer(report_path, cancer_type, patients=None, 
                      filtered_patients=None):
    '''
    Initialize clinical and global variables to instantiate a cancer object.
    
    '''
    run = pickle.load(open('/'.join([report_path, 'RunObject.p']), 'rb'))
    cancer = Cancer(cancer_type, run)
    cancer.patients = patients
    cancer.filtered_patients = filtered_patients
    cancer.path = '/'.join([run.report_path, cancer.name])
    make_path_dump(cancer, cancer.path + '/CancerObject.p')
    
    clinical = Clinical(cancer, run, patients, filtered_patients)
    clinical.artificially_censor(5)
    make_path_dump(clinical, cancer.path + '/Clinical/ClinicalObject.p')
    if type(clinical.drugs) != type(None):
        clinical.drugs.to_csv(cancer.path + '/Clinical/drugs.csv')
    if type(clinical.survival) != type(None):
        clinical.survival.to_csv(cancer.path + '/Clinical/survival.csv')
    clinical.timeline.to_csv(cancer.path + '/Clinical/timeline.csv')
    clinical.clinical.to_csv(cancer.path + '/Clinical/clinical.csv')
    
    global_vars = get_global_vars(cancer.name, run.data_path, patients, 
                                  filtered_patients)
    global_vars.to_csv(cancer.path + '/Global_Vars.csv')