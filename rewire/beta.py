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
from numpy import nanmean, nanstd


#--- retrieve args
parser = argparse.ArgumentParser(
formatter_class=argparse.ArgumentDefaultsHelpFormatter
)
parser.add_argument(
'-inp_txt', '--fname_inp_txt',
type=str,
default='../data/yeast_LIT.txt',
help='ASCII input',
)
parser.add_argument(
'-inp_h5', '--fname_inp_h5',
type=str,
default='./LIT.h5',
help='HDF5 input',
)
pa = parser.parse_args()


# ------- Cargo la informacion del problema ------------- #
fname_graph = pa.fname_inp_txt #'../data/yeast_LIT.txt'
fname_ess   = '../data/Essential_ORFs_paperHe.txt'


graph = igraph.Graph.Read_Ncol(fname_graph, directed=False)
graph.simplify(multiple=True, loops=False) # remover enlaces repetidos

n_nodes = len(graph.vs) # nro total de nodos del grafo
n_edges = len(graph.es) # nro total de enlaces


#--- leamos los fiteos de la distribucion de IBEPs
with h5(pa.fname_inp_h5, 'r') as f:
    #fit_A  = f['fit/A'].value
    fit_mu = f['fit/mu'].value   # valor medio
    fit_sigma = f['fit/sigma'].value  # sigma de la gaussiana

#--- parametros de la red real
N_e  = ff.count_essential_nodes(fname_graph, fname_ess) # nodos esenciales
N_ie = ff.calc_ne(g=graph, fname_ess=fname_ess)  # interacc esenciales

nbad = 0
N_realiz = 10000 # nro de realizaciones
beta     = ff.nans(N_realiz, dtype=np.float32)
overlap  = ff.nans(N_realiz, dtype=np.float32)
for ir in range(N_realiz):
    N_trials, N_overlap, n_e, n_ie = ff.beta_sorting(graph, N_e, N_ie, fit_mu, fit_sigma)
    if N_trials==-1:
        print(" ----> FUCK, llenamos mas q en el caso real!\n \
                Obviemos este caso.")
        nbad += 1

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
