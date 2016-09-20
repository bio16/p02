import igraph 
from copy import deepcopy
import numpy as np
import matplotlib.pyplot as plt
import random as rand

graph = igraph.Graph.Read_Ncol('yeast_AP-MS.txt', directed=False)

"""

number_of_nodes = len(graph.vs)
nodes_to_remove = int(0.35 * number_of_nodes)

graph_aux = graph.clusters()
size_of_large_connected_component =  max(graph_aux.sizes())

fraction_of_nodes = np.array(range(nodes_to_remove), dtype = float) \
                              / number_of_nodes

# ------------------- Remuevo por betweenness ---------------- #

j = 0
size_max_component = np.zeros(nodes_to_remove, dtype = float)
graph_aux = deepcopy(graph)
  
for j in range(nodes_to_remove):

    # Calculo y guardo el size de la componente mas grande
    graph_aux2 = graph_aux.clusters()
    size_max_component[j] +=  float(max(graph_aux2.sizes())) \
                                  / size_of_large_connected_component

    # Calculo el betweenness
    criteria = graph_aux.betweenness()

    # Elijo el nodo con mayor betweenness
    max_crit = max(criteria)
    vertex_ind = criteria.index(max_crit)

    # Remuevo el vertice
    graph_aux.delete_vertices(vertex_ind)

    j += 1

plt.plot(fraction_of_nodes, size_max_component, label = 'Betweenness', linewidth = 2)

# --------------- Remuevo por eigenvector centrality -------------- #

j = 0
size_max_component = np.zeros(nodes_to_remove, dtype = float)
graph_aux = deepcopy(graph)
  
for j in range(nodes_to_remove):

    # Calculo y guardo el size de la componente mas grande
    graph_aux2 = graph_aux.clusters()
    size_max_component[j] +=  float(max(graph_aux2.sizes())) \
                                  / size_of_large_connected_component

    # Calculo el eigenvector centrality
    criteria = graph_aux.evcent(directed = False)

    # Elijo el nodo con mayor betweenness
    max_crit = max(criteria)
    vertex_ind = criteria.index(max_crit)

    # Remuevo el vertice
    graph_aux.delete_vertices(vertex_ind)

    j += 1

plt.plot(fraction_of_nodes, size_max_component, label = 'Eigenvector', linewidth = 2)


# ----------------------- Remuevo por mayor grado ---------------- #

j = 0
size_max_component = np.zeros(nodes_to_remove, dtype = float)
graph_aux = deepcopy(graph)
  
for j in range(nodes_to_remove):

    # Calculo y guardo el size de la componente mas grande
    graph_aux2 = graph_aux.clusters()
    size_max_component[j] +=  float(max(graph_aux2.sizes())) \
                                  / size_of_large_connected_component

    # Calculo el grado
    criteria = graph_aux.degree()

    # Elijo el nodo con mayor grado
    max_crit = max(criteria)
    vertex_ind = criteria.index(max_crit)

    # Remuevo el vertice
    graph_aux.delete_vertices(vertex_ind)

    j += 1

plt.plot(fraction_of_nodes, size_max_component, label = 'Degree', linewidth = 2)

# ------------------------ Remuevo al azar ----------------------- #

total_conf = 100
size_max_component = np.zeros(nodes_to_remove, dtype = float)

for conf in range(total_conf):

  j = 0
  graph_aux = deepcopy(graph)
  
  for j in range(nodes_to_remove):

    # Calculo y guardo el size de la componente mas grande
    graph_aux2 = graph_aux.clusters()
    size_max_component[j] += float(max(graph_aux2.sizes())) \
                                  / size_of_large_connected_component

    vertex_ind = rand.choice(graph_aux.vs)

    # Remuevo el vertice
    graph_aux.delete_vertices(vertex_ind)

    j += 1

size_max_component = size_max_component / total_conf

plt.plot(fraction_of_nodes, size_max_component, label = 'Random', linewidth = 2)

# -------------- Hago el grafico definitivo ------------- #

plt.axis([0, 0.35, 0, 1.00])
plt.legend()
plt.title('AP-MS network')
plt.grid('on')
plt.xlabel('Fraction of nodes')
plt.ylabel('Largest connected component')
plt.show()
"""


