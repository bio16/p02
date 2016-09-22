import igraph 
from copy import deepcopy
import numpy as np
import matplotlib.pyplot as plt
import random as rand
import networkx as nx

# --- Cargo el grafo con igraph ---- #

graph = igraph.Graph.Read_Ncol('yeast_Y2H.txt', directed=False)

# ---- Cargo con Networkx para calcular algunos observables --- #

graph_nx = nx.read_edgelist('yeast_Y2H.txt', create_using = nx.Graph())



# --------------- Lista de esenciales -------------- #

text = open('Essential_ORFs_paperHe.txt','r').read()
text = text.split('\n')
text2 = [aux.split('\t') for aux in text]

essentials = [row[1] for row in text2[2:] if len(row) == 6]

essentials_in_graph = []
for vs in graph.vs:
    vs['essential'] = 0
    for essential in essentials:
        if vs['name'] in essential:
            vs['essential'] = 1
            essentials_in_graph.append(vs['name'])

fraction_of_essentials = float(len(essentials_in_graph))/len(graph.vs)

number_of_nodes = len(graph.vs)

nodes_to_remove = number_of_nodes / 2

# Tamano del componente mas grande
graph_aux = graph.clusters()
size_of_large_connected_component =  max(graph_aux.sizes())

fraction_of_nodes = np.array(range(nodes_to_remove), dtype = float) \
                              / number_of_nodes

# --------------------- Remuevo los esenciales --------------- #

graph_aux = deepcopy(graph)

graph_aux.delete_vertices(essentials_in_graph)

# Calculo y guardo el size de la componente mas grande
graph_aux2 = graph_aux.clusters()
size_max_component =  float(max(graph_aux2.sizes())) \
                      / size_of_large_connected_component

plt.figure(1)
plt.plot(fraction_of_essentials, size_max_component, '.', markersize = 20, label = 'Esenciales')
plt.figure(2)
plt.plot(fraction_of_essentials, size_max_component, '.', markersize = 20, label = 'Esenciales')

# ------------------- Remuevo por betweenness ---------------- #

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

plt.figure(1)
plt.plot(fraction_of_nodes, size_max_component, label = 'Betweenness (iter)', linewidth = 2)

# --------- Lo hago todo de una ----------- #

size_max_component = np.zeros(nodes_to_remove, dtype = float)
graph_aux = deepcopy(graph)
    
# Calculo el betweenness
criteria = graph_aux.betweenness()
criteria_with_names = [[criteria[i], graph.vs[i]['name']] \
                            for i in  range(len(criteria))]

criteria_ordered = sorted(criteria_with_names, reverse = True)

vertex_name_remove = [criteria_ordered[i][1] for i in range(nodes_to_remove)]
    
j = 0
for name in vertex_name_remove:
    
    # Calculo y guardo el size de la componente mas grande
    graph_aux2 = graph_aux.clusters()
    size_max_component[j] +=  float(max(graph_aux2.sizes())) \
                                  / size_of_large_connected_component

    graph_aux.delete_vertices(name)
    j += 1

plt.figure(2)
plt.plot(fraction_of_nodes, size_max_component, label = 'Betweenness', linewidth = 2)

# ------------------- Remuevo por Bonacich ---------------- #

size_max_component = np.zeros(nodes_to_remove, dtype = float)
graph_aux = deepcopy(graph)
  
for j in range(nodes_to_remove):

    # Calculo y guardo el size de la componente mas grande
    graph_aux2 = graph_aux.clusters()
    size_max_component[j] +=  float(max(graph_aux2.sizes())) \
                                  / size_of_large_connected_component

    # Calculo el betweenness
    criteria = graph_aux.evcent(directed = False)

    # Elijo el nodo con mayor betweenness
    max_crit = max(criteria)
    vertex_ind = criteria.index(max_crit)

    # Remuevo el vertice
    graph_aux.delete_vertices(vertex_ind)

plt.figure(1)
plt.plot(fraction_of_nodes, size_max_component, label = 'Bonacich (iter)', linewidth = 2)

# --------- Lo hago todo de una ----------- #

size_max_component = np.zeros(nodes_to_remove, dtype = float)
graph_aux = deepcopy(graph)
    
criteria = graph_aux.evcent(directed = False)
criteria_with_names = [[criteria[i], graph.vs[i]['name']] \
                            for i in  range(len(criteria))]

criteria_ordered = sorted(criteria_with_names, reverse = True)

vertex_name_remove = [criteria_ordered[i][1] for i in range(nodes_to_remove)]
    
j = 0
for name in vertex_name_remove:
    
    # Calculo y guardo el size de la componente mas grande
    graph_aux2 = graph_aux.clusters()
    size_max_component[j] +=  float(max(graph_aux2.sizes())) \
                                  / size_of_large_connected_component

    graph_aux.delete_vertices(name)
    j += 1

plt.figure(2)
plt.plot(fraction_of_nodes, size_max_component, label = 'Bonacich', linewidth = 2)

# ------------------- Remuevo por grado ---------------- #

size_max_component = np.zeros(nodes_to_remove, dtype = float)
graph_aux = deepcopy(graph)
  
for j in range(nodes_to_remove):

    # Calculo y guardo el size de la componente mas grande
    graph_aux2 = graph_aux.clusters()
    size_max_component[j] +=  float(max(graph_aux2.sizes())) \
                                  / size_of_large_connected_component

    # Calculo el betweenness
    criteria = graph_aux.degree()

    # Elijo el nodo con mayor betweenness
    max_crit = max(criteria)
    vertex_ind = criteria.index(max_crit)

    # Remuevo el vertice
    graph_aux.delete_vertices(vertex_ind)

plt.figure(1)
plt.plot(fraction_of_nodes, size_max_component, label = 'Grado (iter)', linewidth = 2)

# --------- Lo hago todo de una ----------- #

size_max_component = np.zeros(nodes_to_remove, dtype = float)
graph_aux = deepcopy(graph)
    
criteria = graph_aux.degree()
criteria_with_names = [[criteria[i], graph.vs[i]['name']] \
                            for i in  range(len(criteria))]

criteria_ordered = sorted(criteria_with_names, reverse = True)

vertex_name_remove = [criteria_ordered[i][1] for i in range(nodes_to_remove)]
    
j = 0
for name in vertex_name_remove:
    
    # Calculo y guardo el size de la componente mas grande
    graph_aux2 = graph_aux.clusters()
    size_max_component[j] +=  float(max(graph_aux2.sizes())) \
                                  / size_of_large_connected_component

    graph_aux.delete_vertices(name)
    j += 1

plt.figure(2)
plt.plot(fraction_of_nodes, size_max_component, label = 'Grado', linewidth = 2)

# --------------- Remuevo por subgraph centrality -------------- #

size_max_component = np.zeros(nodes_to_remove, dtype = float)
graph_aux = deepcopy(graph)
  
for j in range(nodes_to_remove):

    # Calculo y guardo el size de la componente mas grande
    graph_aux2 = graph_aux.clusters()
    size_max_component[j] +=  float(max(graph_aux2.sizes())) \
                                  / size_of_large_connected_component

    criteria = nx.communicability_centrality(graph_nx).items()
    criteria.sort(reverse = True, key = lambda item: item[1])
    # Tomo el primer elemento, de mayor subgraph centrality
    vertex_ind = criteria[0][0]
    
    # Remuevo el vertice
    graph_aux.delete_vertices(vertex_ind)
    graph_nx.remove_node(vertex_ind)

plt.figure(1)
plt.plot(fraction_of_nodes, size_max_component, label = 'SubGraph (iter)', linewidth = 2)

# --------- Lo hago todo de una ----------- #

size_max_component = np.zeros(nodes_to_remove, dtype = float)

subgraph_centrality = nx.communicability_centrality(graph_nx).items()
subgraph_centrality.sort(reverse = True, key = lambda item: item[1])

graph_aux = deepcopy(graph)

vertex_name_remove = [item[0] for item in subgraph_centrality]
    
j = 0
for name in vertex_name_remove:
    
    # Calculo y guardo el size de la componente mas grande
    graph_aux2 = graph_aux.clusters()
    size_max_component[j] +=  float(max(graph_aux2.sizes())) \
                                  / size_of_large_connected_component

    graph_aux.delete_vertices(name)
    j += 1

plt.figure(2)
plt.plot(fraction_of_nodes, size_max_component, label = 'SubGraph', linewidth = 2)

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

plt.figure(1)
plt.plot(fraction_of_nodes, size_max_component, label = 'Random', linewidth = 2)
plt.figure(2)
plt.plot(fraction_of_nodes, size_max_component, label = 'Random', linewidth = 2)

# -------------- Hago el grafico definitivo ------------- #
plt.figure(1)
plt.axis([0, 0.50, 0, 1.00])
plt.legend()
plt.title('Y2H network')
plt.grid('on')
plt.xlabel('Fraction of nodes')
plt.ylabel('Largest connected component')
plt.savefig('Y2H.eps')

plt.figure(2)
plt.axis([0, 0.50, 0, 1.00])
plt.legend()
plt.title('Y2H network')
plt.grid('on')
plt.xlabel('Fraction of nodes')
plt.ylabel('Largest connected component')
plt.savefig('Y2H_b.eps')
