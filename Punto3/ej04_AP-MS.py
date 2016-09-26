import igraph 
from copy import deepcopy
import numpy as np
import matplotlib.pyplot as plt
import random as rand
import networkx as nx

# --- Cargo el grafo con igraph ---- #

graph = igraph.Graph.Read_Ncol('yeast_AP-MS.txt', directed=False)
graph.simplify(multiple = True, loops = False)

# ------- Componente mas grande ----------- #

graph_aux = graph.clusters()
size_of_large_connected_component = max(graph_aux.sizes())

# --------------- Lista de esenciales -------------- #

text = open('Essentials.txt','r').read()
essentials = text.split('\n')

essentials_in_graph = []
for vs in graph.vs:
    vs['essential'] = 0
    for essential in essentials:
        if vs['name'] in essential:
            vs['essential'] = 1
            essentials_in_graph.append(vs['name'])

# ----------- Remuevo esenciales ---------- #

graph_aux = deepcopy(graph)

graph_aux.delete_vertices(essentials_in_graph)

# Calculo y guardo el size de la componente mas grande
graph_aux2 = graph_aux.clusters()
size_max_component =  float(max(graph_aux2.sizes())) \
                      / size_of_large_connected_component

print size_max_component, 
print 'Nodos removidos:', len(essentials_in_graph) 

# ------- Remuevo no esenciales con el mismo grado ------- #

degree_of_essentials = [vs.degree() for vs in graph.vs if vs['essential'] == 1]

# Margen para considerar que dos nodos tienen el mismo grado
delta = 5
same_degree_non_essentials = [[vs['name'] for vs in graph.vs \
                               if vs.degree() <=  (degree + delta) or vs.degree > (degree - delta) \
                               and vs['essential'] == 0] for degree in degree_of_essentials]

data = []
nodes_removed = np.zeros(50)
for i in range(50):

    graph_aux = deepcopy(graph)

    for vertexs in same_degree_non_essentials:

        if vertexs == []:
            continue
        else:
            rand.shuffle(vertexs)
            for vs in vertexs:
                try:
                    graph_aux.delete_vertices(vs)
                    nodes_removed[i] += 1
                    break
                except:
                    pass

    graph_aux2 = graph_aux.clusters()
    size_max_component = float(max(graph_aux2.sizes())) \
                         / size_of_large_connected_component
    data.append(size_max_component)

print np.mean(data), np.std(data)
print 'Nodos no esenciales removidos: ', np.mean(nodes_removed), np.std(nodes_removed)




