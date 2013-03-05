'''
Created on Dec 5, 2012

@author: agross
'''

class DisplayParams:
    def __init__(self, **entries): 
        self.__dict__.update(entries)


def getDisplayDict(mod, similarityMatrix, geneWeights):
    import networkx as nx
    import numpy as np 
    
    s = similarityMatrix.ix[mod,mod].fillna(0)
    keepers = s.index[(s + s.T).sum() / (s + s.T).sum().sum() > .001]
    #outliers = set(mod).difference(keepers)
    s = s.ix[keepers, keepers]
    
    edges = np.array(np.where(s > 1e-5)).T
    edges = [(mod[edge[0]], mod[edge[1]], s.ix[edge[0]][edge[1]]) 
             for edge in edges]
    graph = nx.Graph()
    graph.add_weighted_edges_from(edges)
    
    bestScore = 1e6
    for i in range(5):
        pos = nx.layout.spring_layout(graph, weight='weight', iterations=1000)
        for node in pos:
            pos[node] = (pos[node] * 2) - 1
            radius = np.sqrt(pow(pos[node][0],2) + pow(pos[node][1],2))
            if radius > 1:
                pos[node][0] = pos[node][0]/(radius*1.3)
                pos[node][1] = pos[node][1]/(radius*1.3)
        #score = sum(numpy.sqrt(numpy.abs(pos.values())))
        score = np.sum(pow(np.sum(pos.values(), axis=0),2))
        #score = sum(pow(np.sum(pos.values(), axis=0),2))
        if score < bestScore:
            bestScore = score
            goodPos = pos

    weights = [g[2]['weight'] for g in graph.edges(data=True)]
    if len(weights)>0:
        weights = 4*(weights/max(weights))
    #sizes = np.nan_to_num(np.abs(alteredGenes).ix[graph.nodes()].sum(axis=1))
    sizes = geneWeights.ix[graph.nodes()]
    if sum(sizes < 0) > 0:
        color = sizes.map(lambda s: 1 if s > 0 else 0)
        sizes = sizes.abs()
    else:
        color = sizes
    sizes = 2000*(sizes/(sum(sizes)+.1))
    cover = geneWeights.ix[list(mod)].sum()
    #cover = alteredGenes.ix[list(mod)].sum().sum()
 
    return dict(graph=graph, pos=goodPos, weights=weights, sizes=sizes,
                cover=cover, color=color)
    
    
def drawMod(genes=None, similarityMatrix=None, geneWeights=None,
            ax=None, displayParams=None, center=(0,0)):
    import networkx as nx
    
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    if displayParams is None:
        if (similarityMatrix is None) or (geneWeights is None):
            print 'Need to compute module display parameters.'
            print 'Function needs similarityMatrix and alteredGene matrix.'
            return
        else:
            dP = DisplayParams(**getDisplayDict(genes, similarityMatrix, 
                                                  geneWeights))
    elif type(displayParams) == dict:
        dP = DisplayParams(**displayParams)
    else:
        dP = displayParams
        
    adjPos = dict([(gene, center + dP.pos[gene]) for gene in dP.pos])
    
    if ax is None: 
        fig,ax = plt.subplots(1,1, figsize=(5,5))
    ax.set_frame_on(False)
    ax.set_xticks([])
    ax.set_yticks([])
    #return dP
    nx.draw_networkx_nodes(dP.graph, adjPos, ax=ax, node_size=dP.sizes, 
                           node_color=dP.color, alpha=.5, linewidths=1)
    nx.draw_networkx_edges(dP.graph, adjPos, ax=ax, node_size=dP.sizes, 
                           node_color='grey', alpha=.5, width=dP.weights);
    
    labelPos = dict([(key, adjPos[key]+(0,0)) for key in adjPos])
    nodesToLabel = list(dP.sizes[dP.sizes > 100].index)
    nx.draw_networkx_labels(dP.graph.subgraph(nodesToLabel), labelPos, ax=ax);
    
    art = mpatches.Circle(center, 1, color='black', alpha=.05, fill=False, lw=5)
    ax.add_patch(art)
    
    ax.set_xbound(-1,1)
    ax.set_ybound(-1,1)
    #ax.annotate(int(dP.cover), (center[0]+.75, center[1]+.75))
    
    return dP