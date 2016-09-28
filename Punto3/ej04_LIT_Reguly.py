import igraph 
from copy import deepcopy
import numpy as np
import matplotlib.pyplot as plt
import random as rand
import networkx as nx

# --- Cargo el grafo con igraph ---- #

graph = igraph.Graph.Read_Ncol('yeast_LIT_Reguly_interactions.txt', directed=False)
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

size_max_component = np.zeros(5)
for conf in range(5):

  # Genero un diccionario de (proteina, grado original) con las 
  # no esenciales presentes en la red auxiliar, 
  # pero tomando el grado que tenian en la red original
  non_essentials = {}
  for vs in graph.vs:
      if vs['essential'] == 0:
          non_essentials[vs['name']] = vs.degree()
  non_essentials = non_essentials.items()

  # Mezclo esta lista para no recorrerlos en orden 
  rand.shuffle(degree_of_essentials)

 # Recorro los grados de los esenciales
  vertexs_2_remove = []
  for degree in degree_of_essentials:

    # Genero un lista con el grado de los no esenciales aun no removidas
    degree_of_non_essentials = [item[1] for item in non_essentials]
 
    # Elijo alguna proteina cuyo grado sea el mas cercano posible al grado de una esencial
    id_to_remove = rand.choice([i for i in range(len(degree_of_non_essentials))\
                                if degree_of_non_essentials[i] == \
                                min(degree_of_non_essentials, key = lambda x: abs(x-degree))])

    # Guardo en una lista las que voy a remover, 
    # pero no las remuevo aun del grafo
    vertex2remove = non_essentials[id_to_remove]
    vertexs_2_remove.append(vertex2remove[0])
    non_essentials.remove(vertex2remove)



  graph_aux = deepcopy(graph)
  # Remuevo las proteinas no esenciales con
  # el grado mas parecido posible a las esenciales
  graph_aux.delete_vertices(vertexs_2_remove)
             
  graph_aux2 = graph_aux.clusters()
  size_max_component[conf] = float(max(graph_aux2.sizes())) \
                         / size_of_large_connected_component
  print conf
    
print np.mean(size_max_component), np.std(size_max_component)

