'''
Created on Jan 28, 2013

@author: agross
'''

import os as os
import pickle as pickle

from collections import defaultdict
from pandas import read_table, read_csv
import pandas as pd
from numpy import array, nan

from Data.ProcessClinical import get_clinical
import Data.Firehose as FH
from Processing.Helpers import make_path_dump

def tree(): return defaultdict(tree)

class Run(object):
    '''
    Object for storing meta-data and functions for dealing with Firehose runs.
    '''
    def __init__(self, date, version, data_path, result_path, parameters, description=''):
        self.date = date
        self.data_path = data_path
        self.version = version
        self.result_path = result_path
        self.parameters = parameters
        self.dependency_tree = tree()
        self.description = description
        self.get_meta()
        
        
    def __repr__(self):
        s = 'Run object for TCGA Analysis\n'
        s += 'Firehose run date: ' + self.date + '\n'
        s += 'Code version: ' + self.version + '\n'
        if self.description:
            s += 'Comment: ' + self.description + '\n'
        return s
                 
    def get_meta(self):
        #count_file  = [f for f in os.listdir(self.data_path + 'samples_report') if 
        #               f.startswith('sample_counts')][0]
        self.sample_matrix = read_table(self.data_path + 'meta_data/sample_counts.tsv', 
                                      index_col=0)
        self.data_types = array(self.sample_matrix.columns)
        self.cancers = array(self.sample_matrix.index[:-1])
        
        #self.sample_data = read_csv(self.data_path + 'meta_data/samples.csv', index_col=0)
        self.cancer_codes = read_table(self.data_path + 'meta_data/diseaseStudy.txt', 
                                       index_col=0, squeeze=True)
        self.cancers = list(self.cancer_codes.index)
        
    def load_cancer(self, cancer):
        path = '/'.join([self.report_path, cancer, 'CancerObject.p'])
        obj = pickle.load(open(path, 'rb'))
        return obj
        

class Cancer(object):
    def __init__(self, name, run):
        self.name = name
        try:
            self.full_name = run.cancer_codes.ix[name]
        except:
            self.full_name = name
        counts = run.sample_matrix.ix[name]
        self.samples = counts[counts > 0]
        self.data_types = array(self.samples.index)
        
        #sample_data = run.sample_data[run.sample_data.Disease == name]
        #self.patients = array(sample_data['Participant Number'].drop_duplicates())
    
    def load_clinical(self):
        assert hasattr(self, 'path')
        path = '/'.join([self.path, 'Clinical', 'ClinicalObject.p'])
        obj = pickle.load(open(path, 'rb'))
        return obj
    
    def load_global_vars(self):
        assert hasattr(self, 'path')
        path = '/'.join([self.path, 'Global_Vars.csv'])
        df = read_csv(path, index_col=0)
        df.columns = pd.MultiIndex.from_tuples(map(eval, df.columns))
        return df
    
    def load_data(self, data_type):
        assert hasattr(self, 'path')
        path = '/'.join([self.path, data_type, 'DataObject.p'])
        obj = pickle.load(open(path, 'rb'))
        return obj
        
    def __repr__(self):
        return self.full_name + '(\'' + self.name + '\') cancer object'
    
    
class Clinical(object):
    def __init__(self, cancer, run, patients, filtered_patients):
        self.cancer = cancer.name
        tup = get_clinical(cancer.name, run.data_path, patients, 
                           filtered_patients, **run.clinical_parameters)
        self.clinical, self.drugs, self.followup, self.timeline, self.survival = tup
    
    def __repr__(self):
        return 'Clinical Object for ' + self.cancer
        
    def artificially_censor(self, years):
        for n,s in self.survival.iteritems():
            if n.endswith('y'):
                continue
            df = s.unstack().copy()
            df['event'] = df.event * (df.days < int(365.25*years))
            df['days'] = df.days.clip_upper(int((365.25*years)))
            self.survival[n + '_' + str(years) + 'y'] = df.stack()

def patient_filter(df, can):
    if can.patients is not None:
        return df[[p for p in df.columns if p in can.patients]]
    elif can.filtered_patients is not None:
        return df[[p for p in df.columns if p not in can.filtered_patients]]
    else:
        return df           
    
class Dataset(object):
    def __init__(self, cancer, run, data_type):
        self.data_type = data_type  
        self.path = '/'.join([cancer.path, data_type])
   
        if data_type == 'MAF':
            self.df = FH.get_mutation_matrix(run.data_path, cancer.name)
            self.compressed = False
        elif data_type == 'CN':
            self.df = FH.get_gistic_gene_matrix(run.data_path, cancer.name)
            self.compressed = False
        elif data_type == 'CN_broad':
            self.df = FH.get_gistic(cancer.name, run.data_path,
                                    min_patients=run.parameters['min_patients'])
            self.features = self.df #should probably move
            self.features = patient_filter(self.features, cancer)
            self.compressed = False
        elif data_type == 'mRNA':
            self.df = FH.read_mrna(cancer.name, run.data_path)
            self.compressed = True
        elif data_type == 'mRNASeq':
            self.df = FH.read_rnaSeq(cancer.name, run.data_path, 
                                     average_on_genes=True)
            self.compressed = True
        elif data_type == 'Methylation':
            self.df = FH.read_methylation(cancer.name, run.data_path)
            self.compressed = True
        elif data_type == 'RPPA':
            self.df = FH.read_rppa(run.data_path, cancer.name)
            self.compressed = True
        elif data_type == 'miRNASeq':
            df = FH.read_miRNASeq(cancer.name, run.data_path)
            binary = (df[(df < -1).sum(1) > (df.shape[1]/2)] >= -1)*1.
            binary = binary[binary.sum(1).isin(range(20, df.shape[1]/2))]
            real = df[((df.max(1) - df.min(1)) > 2)]
            real = real.ix[(real == -3).sum(1) < real.shape[1]/2.]
            self.features = pd.concat([real, binary], keys=['real','binary'])
            self.df = df
            self.binary = binary
            self.real = real
            self.compressed = True
        else:
            return
 
        self.df = patient_filter(self.df, cancer)
        self.patients = array(self.df.columns)
        
    def compress(self):
        assert len(self.df.shape) == 2
        self.df = self.df.replace(0, nan).stack()
        if hasattr(self, 'features'):
            self.features = self.features.replace(0, nan).stack()
        self.compressed = True
        
    def uncompress(self):
        assert len(self.df.shape) == 1
        self.df = self.df.unstack().ix[:, self.patients].fillna(0.)
        if hasattr(self, 'features'):
            self.features = self.features.unstack().ix[:, self.patients]
            self.features = self.features.fillna(0.)
        self.compressed = False
        
    def save(self):
        if self.compressed == False:
            self.compress()
        make_path_dump(self, self.path + '/DataObject.p')
            
    def __repr__(self):
        return self.data_type + ' dataset'
            