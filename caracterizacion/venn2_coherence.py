import numpy as np
import argparse as arg
import igraph 
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

argparser = arg.ArgumentParser(description='')
argparser.add_argument('files',help='graph file',nargs=2)

args = argparser.parse_args()

names = []
sets = []       #inicializo una lista de conjuntos
grafos = []     #inicilizo una lista de grafos
for file in args.files: 
    try:
        grafo = igraph.Graph.Read_Ncol(file) 
    except:
        grafo = igraph.Graph.Read_GML(file)
    gset = set(grafo.vs['name'])    # creamos un conjunto de las proteinas del grafo
    sets.append(gset)               # agregamos el conjunto a la lista de conjuntos
    names.append(file.split('/')[-1].split('.')[0])
    grafo = grafo.as_undirected()  
    grafos.append(grafo)            # guardamos el grafo en la lista de grafos
                                    # (ya sabemos que no es dirigido)

inter_all = sets[0] & sets[1]     # Calcula la interseccion de los 2 conjuntos de proiteinas
                                  # en python/conjuntos "&" es el operador de interseccion
print('interseccion: ',len(inter_all))
# Esta seccion es analoga a la creacion de conjuntos anterior pero ahora para los links
# entre proteinas que esten en todos los grafos
link_sets = []
for grafo, name in zip(grafos,names):       
    subgrafo = grafo.subgraph(inter_all)    # crea el subgrafo a partir de la interseccion
    gset = set(subgrafo.get_edgelist())     # crea el conjunto de links del subgrafo
    link_sets.append(gset)                  # guardo el conjunto en link_sets

print('creando diagrama de coherencia...')
venn2(link_sets,names)                      # crea diagrama de coherencia entre links
filename = 'venn2_coherence_'+'-'.join(names)
plt.savefig(filename+'.pdf')
plt.show()

