#!/usr/bin/env ipython 
# -*- coding: utf-8 -*-

import igraph, argparse
import pandas as pd
import random
import matplotlib.pyplot as plt
import numpy as np
from pylab import figure, close, show, find, bar, hist
import funcs as ff


#--- retrieve args
parser = argparse.ArgumentParser(
formatter_class=argparse.ArgumentDefaultsHelpFormatter
)
parser.add_argument(
'-inp', '--fname_inp',
type=str,
default='../data/yeast_LIT.txt',
help='input filename of the network/graph.',
)
parser.add_argument(
'-fig', '--fname_fig',
type=str,
default='./test.png',
help='figure filename',
)
parser.add_argument(
'-bin', '--bins',
type=int,
default=45,
help='number of bins for histogram of IBEPs',
)
parser.add_argument(
'-nr', '--n_rewire',
type=int,
default=100,
help='nro de re-cableados del grafo (manteniendo la distribucion de grado)',
)
parser.add_argument(
'-out', '--fname_out',
type=str,
default='./test.h5',
help='output filename, to save the sequence of the values\
 of the number of essential interactions (IBEPs).',
)
pa = parser.parse_args()


# ------- Cargo la informacion del problema ------------- #

# Cargo en el objeto graph la red de dolphins.gml.
#graph = igraph.read('dolphins.gml')
graph = igraph.Graph.Read_Ncol(pa.fname_inp, directed=False)
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

n_nodes = len(graph.vs) # nro total de nodos del grafo

n_rewire, nbin = pa.n_rewire, pa.bins #1000, 45
he = np.zeros(nbin, dtype=np.int64)
ne_ = []
for i in range(n_rewire):
    ne = ff.calc_ne(graph)
    print(" i, n_essential: ", i, ne)
    #he += np.histogram(ne, bins=nbin, range=[1500.,1900.], normed=False)[0]
    ne_ += [ ne ] # save all values
    graph.rewire()
    #graph.rewire(int(n_nodes/2)) # recablear tantas veces como la mitad
                                 # del nro total de nodos

try:
    # guarda toda la secuencia de nros de interacciones
    # esenciales (IBEPs), para luego hacer histograma en la
    # cantidad de bines q mejor se vea.
    from h5py import File as h5
    fo = h5(pa.fname_out, 'w')
    fo['ne'] = np.array(ne_, dtype=np.int32)
    fo['n_rewire'] = n_rewire
    fo.close()
except ImportError:
    #--- guarda un histograma en ascii
    hc, hx_ = np.histogram(ne_, bins=nbin)
    hx = hx_[1:]-hx_[:-1]
    do = np.array([hx,hc]).T
    np.savetxt(pa.fname_out[:-3]+'.txt', do, fmt='%12.2f')

#--- make fig
fig = figure(1, figsize=(6,4))
ax  = fig.add_subplot(111)

#ax.plot(he, label='$P_E$')
ax.hist(ne_[3:], bins=nbin, label='N:%d'%np.sum(ne_))

ax.grid(True)
ax.legend(loc='best')
#show()
fig.savefig(pa.fname_fig,format='png',dpi=135,bbox_inches='tight')
close()

#EOF
