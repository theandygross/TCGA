'''
Created on Dec 5, 2012

@author: agross
'''

from pandas import read_table, read_csv
from pandas import Series, DataFrame
from numpy import triu, nan, array, sqrt, multiply
from scipy.stats import scoreatpercentile

import pandas.rpy.common as com
from rpy2.robjects.packages import importr

from rpy2.rinterface import RRuntimeError
igraph = importr('igraph')
dynamicTreecut = importr('dynamicTreeCut')
r_stats = importr('stats')

triuFrame = lambda df: DataFrame(triu(df, k=1), index=df.index, columns=df.columns)

def get_graph_weight(graph_file, genes):
    '''
    Takes data frame of columns: gene1, gene2, weight
    '''
    graph = read_csv(graph_file, index_col=0)
    graph = graph[(graph.gene1.isin(genes) * graph.gene2.isin(genes))]
    graph = graph.pivot(index='gene1', columns='gene2', values='weight').sort_index()
    graph_weight = graph.ix[sorted(genes), sorted(genes)]
    graph_weight = triuFrame(graph_weight.add(graph_weight.T, fill_value=0))
    graph_weight = graph_weight.replace(0, nan)
    return graph_weight

def get_data_weight(vec):
    '''
    Takes an geometric mean over an outer product for a data-vector.
    '''
    data_weight = sqrt(multiply.outer(array(vec), array(vec)))
    data_weight = DataFrame(data_weight, index=vec.index, columns=vec.index)
    data_weight = data_weight.replace(data_weight.min().min(), 0)
    data_weight = triuFrame(data_weight)
    return data_weight

def sparsify_df(df, sparsity):
    delta = scoreatpercentile(df.as_matrix().flat, (1-sparsity)*100)
    keeper_genes = ((df > delta).sum() * (df > delta).sum(1)) > 0
    df = df.ix[keeper_genes, keeper_genes]
    df = df * (df > delta)
    return df

def nx_to_igraph(nx_graph):
    '''
    Helper to turn networkx graph into r iGraph object.
    '''
    edges_iGraph = [{'source':edge[0],'target':edge[1],
                     'weight':edge[2]['weight']} 
                     for edge in nx_graph.edges(data=True)]
    graph_frame = DataFrame(edges_iGraph)
    graph_frame_r = com.convert_to_r_dataframe(graph_frame)
    r_graph = igraph.graph_data_frame(graph_frame_r, directed=False)
    return r_graph

def pre_cluster(Gr):
    '''
    Random walk-trap clustering. 
    '''     
    g = nx_to_igraph(Gr)
    wt = igraph.walktrap_community(g, modularity=True)
    labels = array(wt.rx('names'), dtype=int)[0]
    membership = array(wt.rx('membership'), dtype=int)[0]
    membership = Series(membership, index=labels)
    return membership

def recluster_genes(geneset, Gr):
    g = nx_to_igraph(Gr.subgraph(geneset))
    wt = igraph.walktrap_community(g, modularity=True)
    labels = array(wt.rx('names'), dtype=int)[0]
    
    dend = r_stats.as_dendrogram(wt, use_modularity=True)
    hc = r_stats.as_hclust(dend)
    cut = dynamicTreecut.cutreeDynamicTree(hc, maxTreeHeight=1000, 
                                           minModuleSize=12)
    membership = array(cut, dtype=int)
    membership = Series(membership, index=labels)
    return membership

def get_candidates(membership, Gr):   
    counts = membership.value_counts()
    recluster = counts[counts > 30].index
    candidates = [membership[membership==i].index for i,count in 
                  counts.iteritems() if (count <= 30) and (count > 2)]
    for i in recluster:
        geneset = membership[membership==i].index
        try:
            membership_sub = recluster_genes(geneset,Gr)
            sub_candidates = get_candidates(membership_sub, Gr)
            candidates.extend(sub_candidates)
        except RRuntimeError: 
            print 'R could not find a module'
            
    return candidates



'''
from Data.Firehose import read_rnaSeq
#graph_weight = get_graph_weight(graph, mutation_scores.index)
stddata_path = firehose_path + 'stddata__' + date_ + '/' + cancer + '/' + date + '/'
data_matrix = read_rnaSeq(stddata_path)
data_matrix = data_matrix.groupby(by=lambda n: n.split('|')[0]).mean()
data_matrix = data_matrix.ix[graph_weight.index]
data_matrix = data_matrix.dropna(thresh=int(data_matrix.shape[1]*.75))
data_weight = data_matrix.T.corr()
''';

