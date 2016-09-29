#!/bin/env python3

import graph_tool.all as gp
import seaborn as sns
import argparse as arg
import numpy as np
import pandas as pd



# seccion de argumentos del programa
############################################################################################################
####        Parseado de argumentos
############################################################################################################
AcceptedFormats = [ 'auto','gml','gt','graphml','xml','dot','csv' ]     # formatos soportados
vcolors = list(sns.xkcd_rgb.keys()) # nombres de colores (http://xkcd.com/color/rgb/)

argparser = arg.ArgumentParser(description='')

argparser.add_argument('essentials',help='essential list')
argparser.add_argument('data',help='graph file')
argparser.add_argument('--format','-f', default='auto',
        help='force format of the file (default:autodetect)',choices=AcceptedFormats)
argparser.add_argument('--is_directed','-d',action='store_true',help='directed graph')
argparser.add_argument('--threshold','-t',default=3,type=int,help='number of undirect neightbors')

argparser.add_argument('--alpha','-a',default=0,type=float,help='alpha probability (default:None)')
argparser.add_argument('--beta','-b',default=0,type=float,help='beta probability (default:None)')

args = argparser.parse_args()


############################################################################################################
####        filtrado de proteinas escenciales
############################################################################################################
essentials = pd.read_csv(args.essentials, usecols = ['ORF_name'], comment='=',sep='\t',skipfooter=4,
        engine='python',squeeze=True)

for i,node in enumerate(essentials.values):
    node_ = ''.join(node.split('-'))                # eliminamos '-' de los nombres de las proteinas
    essentials[i] = node_.replace(' ','').upper()   # eliminamos los ' ' inutiles de los nombres de las proteinas


############################################################################################################
####        Creacion del grafo
############################################################################################################
if args.format == 'csv':   # si el formato es una lista de links hay que tratarlo distinto (csv)
    graph = gp.load_graph_from_csv(args.data,string_vals=True,directed=args.is_directed,
                csv_options={"delimiter": "\t", "quotechar": "#"})
else:  # si no es csv, que lo lea tranqui..
    graph = gp.load_graph(args.data,fmt=args.format)


############################################################################################################
####        Distribucion de grado (todos los nodos)
############################################################################################################
v_degrees = np.array([ v.out_degree() for v in graph.vertices() ])  # lista de grado por id de nodo
degrees, hist = np.unique(v_degrees,                # degrees: lista de grados existentes en la red
        return_counts=True)                         # hist: cada elemento k es el numero de nodos de grado k

############################################################################################################
####        Distribucion de grado (solo nodos escenciales)
############################################################################################################

#primero marcamos los nodos escenciales
new_prop = graph.new_vertex_property("bool")
graph.vertex_properties['essential'] = new_prop
for v in graph.vertices():
    if graph.vp.name[v] in essentials.values:
        graph.vp.essential[v]=True
    else:
        graph.vp.essential[v]=False


essential_vertex = np.array([ graph.vp.essential[v] for v in graph.vertices() ])  # lista de escencialidad por id de nodo

##### Numero de pares no adyacentes que comparten 3+ vecinos:
#     np = .5 * sum_{i,j}  ( 1 - A(i,j)  ) * F(i,j;3) * type(i,j)
#                          |-- no adya --|
#     
#    F(i,j;3) = 1 if A(j) * A(i) = \sum_k A(i,k)A(k,j)  >= 3
#             = 0 else  
#
#   type(i,j)   = 1 if i,j essentials or non-essentials
#               = 0 else
#

def F(x,y,threshold=3):
    common_neighbors = np.sum( x*y )
    if common_neighbors >= threshold:
        return True
    else:
        return False


A = args.alpha
B = args.beta

def P(k,alpha=A,beta=B):
    return 1 - (1-beta)*(1-alpha)**k



adj = gp.adjacency(graph).toarray()
N = graph.num_vertices()

one = np.ones_like(adj)


disjoin = one - adj

pairs = 0
sameType_pairs = 0
expected_pairs = 0
for i in range(N):
    for j in range(i+1,N):
        if disjoin[i,j] and F(adj[i,:],adj[j,:],args.threshold):
            pairs += 1
            if graph.vp.essential[i] == graph.vp.essential[j]:
                sameType_pairs += 1
            if A and B:
                ki = v_degrees[i]
                kj = v_degrees[j]
                expected_pairs += P(kj)*P(ki) + (1-P(kj))*(1-P(ki))
#for i,arri in enumerate(adj):
#    for j, arrj in enumerate(adj):
#        if i != j and graph.vp.essential[i] == graph.vp.essential[j] and F(arri,arrj) and (1-arri[j]):
#            pairs += 1

print('pairs',pairs)
print('same type pairs',sameType_pairs)
if A and B:
    print('expected pairs',expected_pairs)
