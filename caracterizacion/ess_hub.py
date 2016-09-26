#!/bin/env/python3

from graph_tool.all import * 
import seaborn as sns
import argparse as arg
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import kendalltau, spearmanr

def essential_fraction(df,percent):
    """
    df: DataFrame con columnas 'degree', 'nodes', 'essentials'
    percent: float entre 0 y 1

    retorna 
        hub_frac: fraccion de hubs en la red dado el % de corte
        ess_frac: fraccion de esenciales-hub dado una definicion porcentual de hub

    """
    kmax = df['degrees'].max(axis=0)
    kcut = np.floor(kmax*percent)

    df_hub = df[df['degrees'] >= kcut]

    sum_nodes = df['nodes'].sum(axis=0)
    sum_hubs = df_hub['nodes'].sum(axis=0)
    sum_essentials_hubs = df_hub['essentials'].sum(axis=0)

    hubs_frac = sum_hubs/sum_nodes
    ess_frac = sum_essentials_hubs/sum_hubs

    return hubs_frac, ess_frac

def essential_fraction_array(df,x):
    """
    df: DataFrame con columnas 'degree', 'nodes', 'essentials'
    x : array de porcentajes

    retorna 
        hub_arr: array de fraccion de hubs dado array de porcentajes
        ess_arr: array de fraccion de esenciales dado un array de porcentajes
    """
    ess_arr = np.zeros_like(x)
    hub_arr = np.zeros_like(x)
    for i, percent in enumerate(x):
        hub_arr[i], ess_arr[i] = essential_fraction(df,percent)

    return hub_arr,ess_arr

plt.style.use(['seaborn-talk','seaborn-whitegrid'])
fig,subplot = plt.subplots(ncols=1,nrows=1)

# seccion de argumentos del programa
argparser = arg.ArgumentParser(description='')
argparser.add_argument('essentials',help='essential data')
argparser.add_argument('files',help='graph file',nargs='+')


#argparser.add_argument('--vertex-size','-s',help='set the vertex size (default:5)',default=5,type=float)


args = argparser.parse_args()


essentials = pd.read_csv(args.essentials, usecols = ['ORF_name'], comment='=',sep='\t',skipfooter=4,
        engine='python',squeeze=True)

for i,node in enumerate(essentials.values):
    node_ = ''.join(node.split('-'))                # eliminamos '-' de los nombres de las proteinas
    node_ = node_.upper()                           # setiamos todos en UPPERCASE
    essentials[i] = node_.replace(' ','')           # eliminamos los ' ' inutiles de los nombres de las proteinas




#importa el grafos
data = []
names  = []
for file in args.files:
    try:
        graph = load_graph(file,fmt='gml')
    except OSError:
        graph = load_graph_from_csv(file,string_vals=True,directed=False,
                    csv_options={"delimiter": "\t", "quotechar": "#"})


    name = file.split('/')[-1].split('.')[0]
    names.append(name)

    #primero marcamos los nodos escenciales
    new_prop = graph.new_vertex_property("bool")
    graph.vertex_properties['essential'] = new_prop
    for v in graph.vertices():
        if graph.vp.name[v] in essentials.values:
            graph.vp.essential[v]=True
        else:
            graph.vp.essential[v]=False

    v_degrees = np.array([ v.out_degree() for v in graph.vertices() ])  # lista de grado por id de nodo
    degrees, hist = np.unique(v_degrees,                # degrees: lista de grados existentes en la red
            return_counts=True)                         # hist: cada elemento k es el numero de nodos de grado k


    #creamos histograma de nodos escenciales por grado
    essential_hist = np.zeros_like(hist)
    essential_vertex = np.array([ graph.vp.essential[v] for v in graph.vertices() ])  # lista de escencialidad por id de nodo
    for i,k in enumerate(degrees):
        # (v_degrees == k) es una lista por nodo: 1 si el nodo tiene grado k, 0 si no
        # essential_vertex es una lista de escencialidad de nodo: 1 si el nodo es escencial, 0 si no
        # el producto actua como operador "y"-logico: 1 si es escencial y de grado k
        essential_hist[i] = np.sum(  ( v_degrees == k )*essential_vertex )

    data = pd.DataFrame({'degrees':degrees,'nodes':hist,'essentials':essential_hist}) 

    percent = np.linspace(0,1,100)
    x,y = essential_fraction_array(data,percent)

    subplot.plot(x,y,'-',label=name)

    
    tau, tp_value = kendalltau(x,y) 
    rho, rp_value = spearmanr(x,y)

    print("%25s: %.2f(%g)   %.2f(%.2g)"%(name,tau,tp_value,rho,rp_value))

subplot.set_xlabel('Fraccion de hubs en la red')
subplot.set_ylabel('Fraccion de hubs esenciales')
subplot.legend(loc='best')
plt.savefig('ess_hub.pdf')
plt.show()


