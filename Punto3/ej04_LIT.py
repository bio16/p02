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

essentials_in_graph = list(set(essentials_in_graph))

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

# Hago una lista con los grados de los esenciales
degree_of_essentials = [vs.degree() for vs in graph.vs if vs['essential'] == 1]

data = []
nodes_removed = np.zeros(1)
for i in range(1):

    # Creo una copia del grafo
    graph_aux = deepcopy(graph)

    # Mezclo esta lista para no recorrerlos en orden 
    rand.shuffle(degree_of_essentials)

    # Recorro los grados de los esenciales
    for degree in degree_of_essentials:

        # Genero un diccionario de (proteina, grado original) con las 
        # no esenciales presentes en la red auxiliar, 
        # pero tomando el grado que tenian en la red original
        non_essentials = {}
        for vs in graph.vs:
            if vs['essential'] == 0 and vs['name'] in [vs_aux['name'] for vs_aux in graph_aux.vs]:
                non_essentials[vs['name']] = vs.degree()

        non_essentials = non_essentials.items()

        degree_of_non_essentials = [item[1] for item in non_essentials]

        # Genero una lista de nodos posibles a remover, 
        # tomando aquellos con el grado mas cercano a la esencial
        vertexs2remove = [item[0] for item in non_essentials \
                          if item[1] == min(degree_of_non_essentials, key = lambda x: abs(x-degree))]

        vs_name = rand.choice(vertexs2remove)

        graph_aux.delete_vertices(vs_name)
        nodes_removed[i] += 1
            
    graph_aux2 = graph_aux.clusters()
    size_max_component = float(max(graph_aux2.sizes())) \
                         / size_of_large_connected_component
    
    data.append(size_max_component)

print np.mean(data), np.std(data)
print 'Nodos no esenciales removidos: ', np.mean(nodes_removed), np.std(nodes_removed)




