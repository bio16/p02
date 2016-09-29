#!/usr/bin/env ipython
# -*- coding: utf-8 -*-
from pylab import figure, close, show, find
import numpy as np
from h5py import File as h5
import funcs as ff
import igraph, argparse

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
parser.add_argument(
'-fig', '--fname_fig',
type=str,
default='test.png',
help='figure filename.'
)
pa = parser.parse_args()



# ------- Cargo la informacion del problema ------------- #
fname_ess   = '../data/Essential_ORFs_paperHe.txt'
fname_graph = pa.fname_inp_txt #'../data/yeast_LIT.txt'
graph = igraph.Graph.Read_Ncol(fname_graph, directed=False)
graph.simplify(multiple=True, loops=False) # remover enlaces repetidos

n_nodes = len(graph.vs) # nro total de nodos del grafo

#--- leamos los fiteos de la distribucion de IBEPs
with h5(pa.fname_inp_h5, 'r') as f:
    #fit_A  = f['fit/A'].value
    fit_mu = f['fit/mu'].value   # valor medio
    fit_sigma = f['fit/sigma'].value  # sigma de la gaussiana
    ne = f['ne'][...]

#--- histograma
hc, hx_ = np.histogram(ne, bins=50, density=True) 
hx = 0.5*(hx_[:-1] + hx_[1:])  # `hx` son los valores centrados de c/bin

#--- parametros de la red real
#N_e  = ff.count_essential_nodes(fname_graph, fname_ess) # nodos esenciales
N_ie = ff.calc_ne(g=graph, fname_ess=fname_ess)  # interacc esenciales

alpha_mean = (N_ie - fit_mu)/n_nodes
alpha_err  = fit_sigma/n_nodes

fig = figure(1, figsize=(6,4))
ax  = fig.add_subplot(111)

ax.plot(hx, hc, '-ob', label='realizaciones')
ax.plot(hx, ff.gauss(hx, 1., fit_mu, fit_sigma), '-r', lw=3, alpha=0.6, label='ajuste')
ax.axvline(N_ie, ls='--', c='black', lw=3, label='real')

ax.grid(True)
ax.legend(loc='best')

ax.set_xlabel('numero de interacciones entre proteinas esenciales')
ax.set_ylabel('densidad de probabilidad')

fig.savefig(pa.fname_fig, format='png', dpi=200, bbox_inches='tight')
close(fig)

print('alpha: %3.2f +- %3.2f'%(100.*alpha_mean, 100.*alpha_err))


#EOF
