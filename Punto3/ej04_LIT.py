import igraph 
from copy import deepcopy
import numpy as np
import matplotlib.pyplot as plt
import random as rand
import networkx as nx

# --- Cargo el grafo con igraph ---- #

graph = igraph.Graph.Read_Ncol('yeast_LIT.txt', directed=False)
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

print len(essentials_in_graph)
print len(set(essentials_in_graph))
exit()
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
print len(degree_of_essentials)
degree_of_non_essentials = [vs.degree() for vs in graph.vs if vs['essential'] == 0]

data = []
nodes_removed = 0
for i in range(1):

    graph_aux = deepcopy(graph)

    rand.shuffle(degree_of_essentials)

    for degree in degree_of_essentials:

        degree_of_non_essentials = [vs.degree() for vs in graph_aux.vs if vs['essential'] == 0]

        vs_name = rand.choice([vs['name'] for vs in graph_aux.vs \
                   if vs.degree() == min(degree_of_non_essentials, key = lambda x: abs(x - degree)) \
                   and vs['essential'] == 0])

        graph_aux.delete_vertices(vs_name)
        nodes_removed += 1
            
        
    graph_aux2 = graph_aux.clusters()
    size_max_component = float(max(graph_aux2.sizes())) \
                         / size_of_large_connected_component
    

print size_max_component, nodes_removed
#    data.append(size_max_component)

#print np.mean(data), np.std(data)
#print 'Nodos no esenciales removidos: ', np.mean(nodes_removed), np.std(nodes_removed)




