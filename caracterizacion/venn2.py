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
        grafo = igraph.Graph.Read(file) 
    except:
        grafo = igraph.Graph.Read(file,format='ncol') 

    gset = set(grafo.vs['name'])    # creamos un conjunto de las proteinas del grafo
    sets.append(gset)               # agregamos el conjunto a la lista de conjuntos
    names.append(file.split('/')[-1].split('.')[0])
    grafo = grafo.as_undirected()  
    grafo = grafo.simplify()
    grafos.append(grafo)            # guardamos el grafo en la lista de grafos
                                    # (ya sabemos que no es dirigido)

inter_all = sets[0] & sets[1]     # Calcula la interseccion de los 2 conjuntos de proiteinas
                                  # en python/conjuntos "&" es el operador de interseccion
print('interseccion: ',len(inter_all))
print('creando diagrama de cobertura...')
venn2(sets,names)                      # crea diagrama de cobertura
filename = 'venn2_cobertura_'+'-'.join(names)
plt.savefig(filename+'.pdf')
plt.show()

