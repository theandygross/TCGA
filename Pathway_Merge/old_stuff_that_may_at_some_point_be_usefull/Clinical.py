'''
Created on Oct 9, 2012

@author: agross
'''
from numpy import nan, sort, log
from pandas import Series, DataFrame, notnull, read_csv

from Data.Firehose import get_mutation_matrix
from Data.Firehose import read_rnaSeq, read_mrna
from Data.Intermediate import read_methylation
from Data.Pathways import build_meta_matrix
from Processing.Tests import get_cox_ph, anova, fisher_exact_test, pearson_p
from Processing.Helpers import bhCorrection, extract_pc, drop_first_norm_pc
from Processing.Helpers import cluster_down, df_to_binary_vec, run_rate_permutation


MIN_NUM_HITS = 8
GENE_LENGTHS = read_csv('/cellar/users/agross/Data/GeneSets/coding_lengths.csv', 
                        index_col=0, squeeze=True)
    
def get_tests(clinical, survival_tests, real_variables, binary_variables,
              var_type='boolean'):
    '''
    Gets clinical association tests for a boolean or real valued response 
    variable.
    Does some quick checks to see if the data is sufficent for each test. 
    Note that clinical values are now hard coded into each test function.
    ------------------------------------------------------------------------
    clinical: DataFrame of clinical variables
    covariates: names of covariates to be passed to the cox model,
                (must be columns in clinical DataFrame)
    '''
    check_surv = lambda s: ((s['event_var'] in clinical) and 
                 (len(clinical[[s['event_var'], s['time_var']]].dropna()) > 2))
    make_surv = lambda args: lambda vec: get_cox_ph(clinical, vec, **args)
    surv_tests = dict((test, make_surv(args)) for test,args in 
                      survival_tests.iteritems() if check_surv(args))
    
    check_test = lambda t: (t in clinical) and (notnull(clinical[t]).sum() > 10)
    make_anova = lambda test: lambda vec: anova(vec, clinical[test])
    make_anova_r = lambda test: lambda vec: anova(clinical[test], vec)
    make_fisher = lambda test: lambda vec: fisher_exact_test(vec, clinical[test])
    make_pcc = lambda test: lambda vec: pearson_p(vec, clinical[test])
    
    if var_type == 'boolean':
        real_test, bin_test = make_anova, make_fisher
    elif var_type == 'real':
        real_test, bin_test = make_pcc, make_anova_r
            
    real_tests = dict((test, real_test(test)) for test in real_variables
                      if check_test(test))    
    bin_tests = dict((test, bin_test(test)) for test in binary_variables
                      if check_test(test))
    
    return dict(list(surv_tests.items()) + list(real_tests.items()) + 
                list(bin_tests.items()))

def run_tests(tests, data_frame):
    '''
    Runs each test in tests across the rows of the DataFrame.
    Returns the p-values, and the corrected q-values.
    ------------------------------------------------------------
    tests: dictionary mapping test name to test functions
    data_frame: DataFrame of response variables to test against
    '''
    if len(data_frame) == 0:
        return DataFrame(), DataFrame()
    p_values = DataFrame()
    for test, f in tests.iteritems():
        p_vec = Series(dict((p, f(vec)) for p, vec in data_frame.iterrows()), 
                       name=test)
        p_values = p_values.join(p_vec, how='outer')
    q_values = p_values.apply(bhCorrection)
    return p_values, q_values

def run_clinical_bool(cancer, clinical, data_path, gene_sets, 
                      survival_tests, real_variables, binary_variables,
                      data_type='mutation'):
    '''
    Runs clinical tests for boolean type data (mutation, amplification, 
    or deletion).
    '''
    if data_type == 'mutation':
        hit_matrix, _ = get_mutation_matrix(cancer, data_path)
    elif data_type in ['amplification', 'deletion']:
        hit_matrix, lesion_matrix = get_cna_matrix(cancer, data_path, data_type)
    
    single_matrix = lesion_matrix if data_type == 'amplification' else hit_matrix
    single_matrix = single_matrix.clip_upper(1.)
    clinical['rate'] = log(single_matrix.sum(0)).clip_lower(0)
    tests = get_tests(clinical, survival_tests, real_variables, binary_variables,
                      var_type='boolean')
    gene_counts = sort(single_matrix.sum(1))
    good_genes = gene_counts[gene_counts > MIN_NUM_HITS].index[-500:]
    p_genes, q_genes = run_tests(tests, single_matrix.ix[good_genes])
    
    clinical['rate'] = log(hit_matrix.sum(0)).clip_lower(0)
    meta_matrix = build_meta_matrix(gene_sets, hit_matrix, 
                                    set_filter=lambda s: s.clip_upper(1))
    tests = get_tests(clinical, survival_tests, real_variables, binary_variables,
                      var_type='boolean')
    p_pathways, q_pathways = run_tests(tests, meta_matrix.clip_upper(1.))
    
    lengths = GENE_LENGTHS.ix[hit_matrix.index].dropna()
    if 'rate' in tests:
        _res = run_rate_permutation(meta_matrix, hit_matrix, gene_sets, 
                                   lengths, tests['rate'])
        q_pathways['rate']  = _res 
    return locals()

def run_clinical_real(cancer, clinical, data_path, gene_sets,
                      survival_tests, real_variables, binary_variables,
                      data_type='expression', drop_pc=False):
    
    if data_type == 'expression':
        data_matrix = read_rnaSeq(cancer, data_path)
        data_matrix = data_matrix.groupby(by=lambda n: n.split('|')[0]).mean()
    elif data_type == 'expression_array':
        data_matrix = read_mrna(cancer, data_path)
    elif data_type == 'methylation':
        data_matrix = read_methylation(cancer, data_path)
    if drop_pc:
        data_matrix = drop_first_norm_pc(data_matrix)
    pc = dict((p, extract_pc(data_matrix.ix[g])) for p, g in 
              gene_sets.iteritems())
    pc = DataFrame(dict((p, (v - v.mean()) / v.std()) for p,v in pc.iteritems() if 
                   type(v) != type(None))).T
    #clinical['pc'] = extract_pc(data_matrix.dropna(), pc_threshold=0)
    tests  = get_tests(clinical, survival_tests, real_variables, 
                       binary_variables, var_type='real')
    #return locals()
    p_pathways, q_pathways = run_tests(tests, pc)
    return locals()

def create_secondary_clinical_features(clinical):
    '''
    Does a little bit of processing that does not belong in the clinical file.
    Right now we artificially censor at 3 and 5 years.
    '''
    clinical['deceased_5y'] = (clinical.days < (365.25*5)) * clinical.deceased
    clinical['days_5y'] = clinical.days.clip_upper(int(365.25*5))
    clinical['deceased_3y'] = (clinical.days < (365.25*3)) * clinical.deceased
    clinical['days_3y'] = clinical.days.clip_upper(int(365.25*3))

    if hasattr(clinical, 'event'):
        clinical['event'] = clinical[['event','deceased']].sum(1).clip_upper(1.)
        evs = clinical.event_free_survival
        clinical['event_5y'] = (evs < (365.25*5)) * clinical.event
        clinical['event_free_survival_5y'] = evs.clip_upper(int(365.25*5))
        clinical['event_3y'] = (evs < (365.25*3)) * clinical.event
        clinical['event_free_survival_3y'] = evs.clip_upper(int(365.25*3))
    return clinical

class ClinicalObject(object):
    '''
    Wrapper to store the data in a nice object rather than have to unpack it. 
    '''
    def __init__(self, cancer, data_path, gene_sets, data_type='mutation',
                 drop_pc=False, survival_tests={}, real_variables=[],
                 binary_variables=[]):
        clinical = read_csv(data_path + '/'.join(['ucsd_processing', cancer, 
                                      'Clinical','compiled.csv']), index_col=0)
        clinical = create_secondary_clinical_features(clinical)
    
        if data_type in ['mutation', 'amplification','deletion']:
            o_dict = run_clinical_bool(cancer, clinical, data_path, gene_sets, 
                                       survival_tests, real_variables, 
                                       binary_variables, data_type)
        else:
            o_dict = run_clinical_real(cancer, clinical, data_path, 
                                       gene_sets, survival_tests, real_variables, 
                                       binary_variables, data_type, drop_pc)
        self.__dict__ = o_dict
        self.data_path = data_path
        #self.patients = clinical[['days','censored']].dropna().index
        
        
    def filter_bad_pathways(self, gene_lookup):
        for clin_type, p_vals in self.p_genes.iteritems():
            for gene in set(p_vals.index).intersection(set(gene_lookup)):
                for pathway in gene_lookup[gene]:
                    if ((pathway in self.q_pathways[clin_type]) and 
                       (p_vals[gene] < self.p_pathways[clin_type][pathway]) and
                        (self.hit_matrix.ix[gene].sum() > 
                         self.meta_matrix.ix[pathway].sum()*.5)):
                        self.p_pathways[clin_type][pathway] = nan
            self.q_pathways[clin_type] = bhCorrection(self.p_pathways[clin_type])
            
    def cluster_down(self, num_clusters=50, draw_dendrogram=False):
        if self.data_type in ['mutation', 'amplification','deletion']:
            matrix = (self.meta_matrix > 0).astype(int)
            agg_function = df_to_binary_vec
            dist_metric = 'jaccard'
        else:
            matrix = self.pc
            agg_function = extract_pc
            dist_metric = 'euclidean'
        r = cluster_down(matrix, agg_function, dist_metric, num_clusters,
                         draw_dendrogram)
        self.clustered, self.assignments = r[0], r[1]
        if draw_dendrogram:
            return r[2]
        
        