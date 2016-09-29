#!/usr/bin/env ipython 
# -*- coding: utf-8 -*-

import igraph, argparse
import pandas as pd
import random
import matplotlib.pyplot as plt
import numpy as np
from pylab import figure, close, show, find, bar, hist
from h5py import File as h5
import funcs as ff
from numpy.random import normal as rand_norm # v.a. normal
from numpy import nanmean, nanstd


#--- retrieve args
parser = argparse.ArgumentParser(
formatter_class=argparse.ArgumentDefaultsHelpFormatter
)
parser.add_argument(
'-fig', '--fname_fig',
type=str,
default='./test.png',
help='figure filename',
)
pa = parser.parse_args()


# ------- Cargo la informacion del problema ------------- #
fname_graph = '../data/yeast_LIT.txt'
fname_ess   = '../data/Essential_ORFs_paperHe.txt'


graph = igraph.Graph.Read_Ncol(fname_graph, directed=False)
graph.simplify(multiple=True, loops=False) # remover enlaces repetidos

n_nodes = len(graph.vs) # nro total de nodos del grafo
n_edges = len(graph.es) # nro total de enlaces


#--- leamos los fiteos de la distribucion de IBEPs
with h5('./test.h5', 'r') as f:
    #fit_A  = f['fit/A'].value
    fit_mu = f['fit/mu'].value   # valor medio
    fit_sigma = f['fit/sigma'].value  # sigma de la gaussiana

#--- parametros de la red real
N_e  = ff.count_essential_nodes(fname_graph, fname_ess) # nodos esenciales
N_ie = ff.calc_ne(g=graph, fname_ess=fname_ess)  # interacc esenciales

nbad = 0  # nro de simulaciones "irreales"
N_realiz = 10000 # nro de realizaciones
beta     = ff.nans(N_realiz, dtype=np.float32)
overlap  = ff.nans(N_realiz, dtype=np.float32)
for ir in range(N_realiz):
    # todos los enlaces no-esenciales por defecto
    for id in range(n_edges):
        graph.es[id]['essential'] = 0
    # todos los nodos no-esenciales por defecto
    for node_id in range(len(graph.vs)):
        graph.vs[node_id]['essential'] = 0

    # nro entero aleatorio con distrib normal, dado los 
    # parametros de fiteo de la distribucion del nro de 
    # interacc esenciales (IBEPs).
    n_ie = int(rand_norm(loc=fit_mu,scale=fit_sigma)) #rn()

    for i in range(N_ie-n_ie):
        # asignar enlace esencial `i`
        id = np.random.randint(n_edges) # sorteamos el enlace
        graph.es[id]['essential'] = 1

    # deducir y marcar los nodos esenciales q se generan por
    # culpa de estos enlaces esenciales nuevos
    ff.make_essential_nodes(graph)
    # contemos el nro de nodos ess. agregados!
    n_e = ff.count_essential_nodes(graph=graph) 
    print(" -> falta completar %d nodos esenciales."%(N_e-n_e))
    if N_e < n_e: # muy pocas veces, se supera el caso real
        print(" ----> FUCK, llenamos mas q en el caso real!\n \
                Obviemos este caso.")
        nbad += 1
        continue

    # hacer intentos para atribuir nodos esenciales
    # por causas random. Devuelve el nro total de intentos, y
    # el nro de overlapping con el etiquetado anterior (por 
    # distribuc de IBEPS).
    N_trials, N_overlap = ff.make_random_essential_nodes(graph, N_e)
    # NOTE: notar q debe ser `N_trials >= N_e-n_e`, es decir,
    #       el nro de intentos sera igual o mayor q el nro de 
    #       nodos esenciales a completar, hasta llegar a 
    #       `N_e` (caso real).
    print(" -> se completó con %d intentos."%N_trials)

    beta[ir]    = 1.*N_trials/n_nodes
    # nro de veces q intentamos asignar esencialidad a nodos
    # q ya eran esenciales
    overlap[ir] = 1.*N_overlap/n_e
    print(" -> beta: %3.2f %%"%(100.*beta[ir]))

print(" --> nro de realizaciones irreales: %d/%d"%(nbad,N_realiz))
print(" --> beta: (%3.2f +- %3.2f) %% "%(100.*nanmean(beta), 100.*nanstd(beta)))
print(" --> overlap: (%3.2f +- %3.2f) %% "%(100.*nanmean(overlap), 100.*nanstd(overlap)))

#EOF
