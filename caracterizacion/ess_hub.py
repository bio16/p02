#!/bin/env/python3

from graph_tool.all import * 
import seaborn as sns
import argparse as arg
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import kendalltau, spearmanr

plt.style.use(['seaborn-talk','seaborn-whitegrid'])
fig,subplot = plt.subplots(ncols=1,nrows=1)

# seccion de argumentos del programa
argparser = arg.ArgumentParser(description='')
argparser.add_argument('essentials',help='essential data')
argparser.add_argument('files',help='graph file',nargs='+')
argparser.add_argument('--percent-cut','-p',help='% hubs definition (default:0.2)',default=.2,type=float)
args = argparser.parse_args()


#definicion de funciones utiles
def essential_fraction(df,percent,v_degrees, return_data_cut=False):
    """
    df: DataFrame con columnas 'degree', 'nodes', 'essentials'
    percent: float entre 0 y 1
    v_degrees: lista de grado de cada nodo

    retorna 
        hub_frac: fraccion de hubs en la red dado el % de corte
        ess_frac: fraccion de esenciales-hub dado una definicion porcentual de hub

    """
    sum_nodes = df['nodes'].sum(axis=0)             # numero total de nodos
    ncut = int(np.floor((sum_nodes-1)*percent))     # numero de nodos 'hub' en el corte %
    kcut = np.sort(v_degrees)[::-1][ncut]           # k del nodo de corte
    

    df_hub = df[df['degrees'] >= kcut]              # tomamos solo los nodos con k >= kcut

    if return_data_cut:                             # Retornamos los histogramas y kcut si eso nos interesa
        return df_hub['nodes'], df_hub['essentials'], kcut

    sum_hubs = df_hub['nodes'].sum(axis=0)                  # real numero de hubs dado kcut
    sum_essentials_hubs = df_hub['essentials'].sum(axis=0)  # numero de esenciales dado un kcut

    hubs_frac = sum_hubs/sum_nodes              # fraccion de hubs en la red
    ess_frac = sum_essentials_hubs/sum_hubs     # fraccion de esenciales-hub en la red
    

    return hubs_frac, ess_frac

def essential_fraction_array(df,x,v_degrees):
    """
    df: DataFrame con columnas 'degree', 'nodes', 'essentials'
    x : array de porcentajes
    v_degrees: lista de grado de cada nodo

    retorna 
        hub_arr: array de fraccion de hubs dado array de porcentajes
        ess_arr: array de fraccion de esenciales dado un array de porcentajes
    """
    ess_arr = np.zeros_like(x)
    hub_arr = np.zeros_like(x)
    for i, percent in enumerate(x):
        hub_arr[i], ess_arr[i] = essential_fraction(df,percent,v_degrees)

    return hub_arr,ess_arr



# importacion de nodos esenciales
essentials = pd.read_csv(args.essentials, usecols = ['ORF_name'], comment='=',sep='\t',skipfooter=4,
        engine='python',squeeze=True)

for i,node in enumerate(essentials.values):
    node_ = ''.join(node.split('-'))                # eliminamos '-' de los nombres de las proteinas
    node_ = node_.upper()                           # setiamos todos en UPPERCASE
    essentials[i] = node_.replace(' ','')           # eliminamos los ' ' inutiles de los nombres de las proteinas




print("%25s: tau(pvalue)\trho(pvalue)\tkcut"%' ')
for file in args.files:  #para cada red
    # Crear grafo
    try:
        graph = load_graph(file,fmt='gml')
    except OSError:
        graph = load_graph_from_csv(file,string_vals=True,directed=False,
                    csv_options={"delimiter": "\t", "quotechar": "#"})


    name = file.split('/')[-1].split('.')[0]        # nombre filtrado

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

    data = pd.DataFrame({'degrees':degrees,'nodes':hist,'essentials':essential_hist}) # histogramas como funcion de k


    percent = np.linspace(0,1,100)                          #  
    x,y = essential_fraction_array(data,percent,v_degrees)  #
    subplot.plot(x,y,'-',label=name)                        # 



    # Medicion de correlacion    
    x,y ,kcut= essential_fraction(data,args.percent_cut,v_degrees,return_data_cut=True)

    tau, tp_value = kendalltau(x,y) 
    rho, rp_value = spearmanr(x,y)

    print("%25s: %.2f(%g)\t%.2f(%.2g)\t%3i"%(name,tau,tp_value,rho,rp_value,kcut))

subplot.set_xlabel('Fraccion de hubs en la red',fontsize=20)
subplot.set_ylabel('Fraccion de hubs esenciales',fontsize=20)
subplot.tick_params(labelsize=20)
subplot.legend(loc='best')
plt.savefig('ess_hub.pdf')
plt.show()


