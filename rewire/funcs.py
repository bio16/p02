#!/usr/bin/env ipython
# -*- coding: utf-8 -*-
"""
 algunas funciones para des-fragmentar el analisis
"""
import numpy as np
from pylab import figure, close, show, find
from h5py import File as h5
from numpy import (
exp, sqrt, power, nan, array, ones, zeros
)
import igraph, os

M_PI = np.pi
M_E  = np.e


def calc_ne(g=None, fname_ess=None):
    """
    calcula en nro total de interacciones `ne` (IBEPs) entre nodos
    esenciales, para una realizacion de red `g` dada.
    """
    assert fname_ess is not None or 'essential' in g.vs.attributes(),\
        ' --> falta especificar los nodos esenciales!!'

    if fname_ess is not None:
        list_raw = open(fname_ess,'r').readlines()[2:-4]
        list_ess = [ nm for nm in 
             [ list_raw[i].split('\t')[1].replace(' ','') for i in range(len(list_raw)) ]
           ]
        for vs in g.vs:
            vs['essential'] = 1 if vs['name'] in list_ess else 0

    #NOTE: la matriz de adyacencia me da la informacion
    #      si esta o no conectado con cierto nodo
    #A = g.get_adjacency()
    # list de vecinos `nn_list[i]` para cada nodo `i`
    nn_list = g.get_adjlist()
    ne = 0 # nro de interacc entre nodos esenciales
    for vs in g.vs:
        if vs['essential']: 
            #neighbor_index = find(A[vs.index,:]) # indices de todos los vecinos de este esencial
            neighbor_index = nn_list[vs.index]
            for nn in neighbor_index: # contemos los vecinos q son esenciales
                ne += g.vs[nn]['essential'] # >0 only for essential nodes

    return ne

def gauss_normal(x):		# esto tiene integral =1.
	phi 	= exp(-.5 * x**2.) / sqrt(2.*M_PI)
	return phi

def gauss(x, A, mu, sig):	# esto tmb tiene integral =1.
	return A * (1./sig) * gauss_normal((x - mu) / sig)

def fit_ne_distribution(fname_inp, nbin=50):
    """
    - leer la coleccion de valores del nro de interacc esenciales
    - construir distribucion de dichos valores y fitear con gaussiana
    """
    from scipy.optimize import curve_fit
    # leemos la data de la simulacion anterior (la coleccion de valores
    # del nro de interacciones esenciales)
    with h5(fname_inp,'r') as f:
        ne = f['ne'][3:]  # obviamos sus primeros tres valores para 
                          # olvidarnos un poco de la configuracion original
                          # NOTE: la config original corresponde al `f['ne'][0]`.

    # hallamos su distribucion frecuentista
    # NOTE:
    #   - `hx_` son los bordes del dominio bineado
    #   - density=True es para q devuelva una distribucion con area=1.
    hc, hx_ = np.histogram(ne, bins=nbin, density=True) 
    hx = 0.5*(hx_[:-1] + hx_[1:])  # `hx` son los valores centrados de c/bin
    hx_mean = (hx*hc).sum()/hc.sum()

    # hacemos el ajuste de la distribucion frecuentista para hallar
    # su densidad de probabilidad
    popt, pcov = curve_fit(
        f = gauss, 
        xdata = hx,
        ydata = hc,
        p0 = [1., hx_mean, 30.], # semillas
    )
    return popt, pcov

def count_essential_nodes(fname_graph=None, fname_ess=None, graph=None):
    """
    cuento el nro de nodos esenciales del grafo 
    en `fname_graph`, de acuerdo a la lista
    de esenciales en ´fname_ess´.
    NOTA:
    la informacion de esencialidad lo sacamos o bien de 
    los archivos input o bien del objeto graph.
    """
    assert fname_graph is not None or graph is not None,\
        " --> grafo no especificado!!"

    if fname_graph is not None and fname_ess is not None:
        # grab data
        graph = igraph.Graph.Read_Ncol(fname_graph, directed=False)
        graph.simplify(multiple=True, loops=False) # remover enlaces repetidos

        list_raw = open(fname_ess,'r').readlines()[2:-4]
        list_ess = [ nm for nm in 
             [ list_raw[i].split('\t')[1].replace(' ','') for i in range(len(list_raw)) ]
           ]

    # ahora si contemos
    n = 0
    if fname_ess is not None:
        for vs in graph.vs:
            n += 1 if vs['name'] in list_ess else 0
    elif graph is not None:
        for vs in graph.vs:
            n += vs['essential']
    else:
        raise SystemExit(" --> wrong arguments!")

    return n

def make_essential_nodes(g):
    """
    deducir nodos esenciales a partir de los
    nuevos enlaces esenciales
    """
    for es in g.es:
        if es['essential']:
            for node_id in es.tuple:
                g.vs[node_id]['essential'] = 1

def make_random_essential_nodes(g, N_e):
    """
    calcular nro de intentos para atribuir nodos esenciales por
    causas random (i.e. otras causas ademas de interacciones
    esenciales)
    Ne : nro de nodos esenciales en la red real
    g  : objeto `igraph.Graph`
    """
    n_nodes = len(g.vs)
    __now__N_e = count_essential_nodes(graph=g)
    #assert N_e>=__now__N_e and 'essential' in g.vs.attributes(),\
    #assert 'essential' in g.vs.attributes(),\
    #    " --> falta flaggear nodos esenciales!!"

    if N_e == __now__N_e: # nada q hacer
        return 0

    n_fill = N_e - __now__N_e # nro de nodos efectivos a llenar
    N_trials, nok = 0, 0
    # NOTE: intentamos transformar cada nodo en esencial hasta
    #       q `nok` sea igual al nro de nodos `n_fill` q faltan
    #       transformar.
    while nok < n_fill:
        node_id = np.random.randint(n_nodes)
        if g.vs[node_id]['essential']==0:
            nok += 1 # transformacion efectiva

        g.vs[node_id]['essential'] = 1
        N_trials += 1


    return N_trials


#EOF
