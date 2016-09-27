#!/usr/bin/env python 
# -*- coding: utf-8 -*-

# Uso la librería igraph y pandas.
import igraph
import pandas as pd
import random
import matplotlib.pyplot as plt
import numpy as np
from pylab import figure, close, show, find


#--- algunas funciones para des-fragmentar el analisis
def calc_ne(g=None):
    #NOTE: la matriz de adyacencia me da la informacion
    #      si esta o no conectado con cierto nodo
    A = g.get_adjacency()

    ne = 0 # nro de interacc entre nodos esenciales
    for vs in g.vs:
        if vs['essential']: 
            neighbor_index = find(A[vs.index,:]) # indices de todos los vecinos de este esencial
            for nn in neighbor_index: # contemos los vecinos q son esenciales
                ne += 1 if g.vs[nn]['essential'] else 0

    return ne

# ------- Cargo la informacion del problema ------------- #

# Cargo en el objeto graph la red de dolphins.gml.
#graph = igraph.read('dolphins.gml')
graph = igraph.Graph.Read_Ncol('../data/yeast_LIT.txt', directed=False)
graph.simplify(multiple=True, loops=False)

list_raw = open('../data/Essential_ORFs_paperHe.txt','r').readlines()[2:-4]
list_ess = [ nm for nm in 
             [ list_raw[i].split('\t')[1].replace(' ','') for i in range(len(list_raw)) ]
           ]

for vs in graph.vs:
    #print(vs['name'])
    #if vs['name'].startswith('YAL0'): input()
    #vs['essential'] = 1 if vs['name'] in list_ess else 0
    vs['deg'] = vs.degree()
    if vs['name'] in list_ess:
        vs['essential'] = 1
    else:
        vs['essential'] = 0

n_rewire, nbin = 1000, 30
he = np.zeros(nbin, dtype=np.int64)
ne_ = []
for i in range(n_rewire):
    ne = calc_ne(graph)
    print(" i, n_essential: ", i, ne)
    he += np.histogram(ne, bins=nbin, range=[1500.,1900.], normed=False)[0]
    ne_ += [ ne ] # save all values
    graph.rewire()

fig = figure(1, figsize=(6,4))
ax  = fig.add_subplot(111)

ax.plot(he, label='$P_E$')

ax.grid(True)
ax.legend(loc='best')
show()
close()


#EOF
