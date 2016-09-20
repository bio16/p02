import igraph 
from copy import deepcopy
import numpy as np
import matplotlib.pyplot as plt
import random as rand
from scipy import linalg

graph = igraph.Graph.Read_Ncol('yeast_AP-MS.txt', directed=False)

def subgraph_centrality(graph_type):

    A = graph_type.get_adjacency()
    A = [A[i] for i in range(A.shape[0])]
    
    eigvalues, eigvectors = linalg.eig(A)

    sc = []
    for i in range(len(graph_type.vs)):
        aux = []
        for j in range(len(eigvectors)):
            aux.append((eigvectors[j][i] ** 2) * np.exp(eigvalues[j]))
        sc.append(np.sum(aux))

    return sc


print subgraph_centrality(graph)
        



