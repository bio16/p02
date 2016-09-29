#!/bin/env python3

from graph_tool.all import *
import seaborn as sns
import argparse as arg
import lmfit as lmf
import matplotlib.pyplot as plt
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
argparser.add_argument('--data-color','-c',help='set the data color on the plot (xkcd palette http://xkcd.com/color/rgb/)',
        default='red')

argparser.add_argument('--kcut','-k',default=None,help='Max degree to take account (default:kmax)')
argparser.add_argument('--bins-number','-b',default=None,help='Number of bins (default: sqrt(kmax))')

args = argparser.parse_args()


############################################################################################################
####        filtrado de proteinas escenciales
############################################################################################################
essentials = pd.read_csv(args.essentials, usecols = ['ORF_name'], comment='=',sep='\t',skipfooter=4,
        engine='python',squeeze=True)

for i,node in enumerate(essentials.values):
    node_ = ''.join(node.split('-'))                # eliminamos '-' de los nombres de las proteinas
    essentials[i] = node_.replace(' ','')           # eliminamos los ' ' inutiles de los nombres de las proteinas


############################################################################################################
####        Creacion del grafo
############################################################################################################
if args.format == 'csv':   # si el formato es una lista de links hay que tratarlo distinto (csv)
    graph = load_graph_from_csv(args.data,string_vals=True,directed=args.is_directed,
                csv_options={"delimiter": "\t", "quotechar": "#"})
else:  # si no es csv, que lo lea tranqui..
    graph = load_graph(args.data,fmt=args.format)


############################################################################################################
####        Distribucion de grado (todos los nodos)
############################################################################################################
v_degrees = np.array([ v.out_degree() for v in graph.vertices() ])  # lista de grado por id de nodo
degrees, hist = np.unique(v_degrees,                # degrees: lista de grados existentes en la red
        return_counts=True)                         # hist: cada elemento k es el numero de nodos de grado k

max_degree = np.max(degrees)                # Grado maximo en la red
if args.kcut != None:
    max_degree = int(args.kcut)
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

#creamos histograma de nodos escenciales por grado
essential_hist = np.zeros_like(hist)
essential_vertex = np.array([ graph.vp.essential[v] for v in graph.vertices() ])  # lista de escencialidad por id de nodo
for i,k in enumerate(degrees):
    # (v_degrees == k) es una lista por nodo: 1 si el nodo tiene grado k, 0 si no
    # essential_vertex es una lista de escencialidad de nodo: 1 si el nodo es escencial, 0 si no
    # el producto actua como operador "y"-logico: 1 si es escencial y de grado k
    essential_hist[i] = np.sum(  ( v_degrees == k )*essential_vertex )


############################################################################################################
####        Logaritmic binning
############################################################################################################

#Extension de los histogramas a todos los grados entre 0 y max_degree
all_degrees = np.arange(0,max_degree+1)
all_degrees_hist = np.zeros_like(all_degrees,dtype=float)
all_degrees_essential_hist = np.zeros_like(all_degrees,dtype=float)
for k,tnodes,enodes in zip(degrees,hist,essential_hist):
    try:
        all_degrees_hist[k] = tnodes
        all_degrees_essential_hist[k] = enodes
    except IndexError:
        pass


nbins = args.bins_number                    # Numero de bins pasado por argumentos
if nbins == None:                           # si no
    nbins = np.sqrt(graph.num_vertices())             # tomar nbins = sqrt(kmax)
    
# creamos una escala logaritmica para los bins
log_bins = np.unique(np.ceil(np.logspace(0,np.log10(max_degree), nbins)))
log_bin_centers = np.array([ .5*(log_bins[i]+log_bins[i+1]) for i in range(len(log_bins)-1) ])

# rellenamos los log_bins segun numero de nodos que caen en el
log_hist = np.zeros_like(log_bin_centers)
log_essential_hist = np.zeros_like(log_bin_centers)
for i in range(len(log_bins)-1):
    kmin = int(log_bins[i])     # limite inferior del bin
    kmax = int(log_bins[i+1])   # limite superior del bin
    if kmax == max_degree :     # correccion si es el ultimo bin
        kmax += 1               # ^
    
    for k in range(kmin,kmax):
        log_hist[i] += all_degrees_hist[k]
        log_essential_hist[i] += all_degrees_essential_hist[k]
    

non_essential_probability_by_degree = 1 - log_essential_hist/log_hist

############################################################################################################
####        Fiteo
############################################################################################################
#fiteo linear del log de los datos
lin_model = lmf.models.LinearModel()
params = lin_model.guess(np.log(non_essential_probability_by_degree),x=log_bin_centers)
result = lin_model.fit(np.log(non_essential_probability_by_degree),params,x=log_bin_centers)
print(result.fit_report())
print()
slope = result.params['slope']
intercept = result.params['intercept']
print('alpha =', 1 - np.exp( slope.value) ,'+-',np.abs(np.exp( slope.value )*slope.stderr ))
print('beta  =', 1 - np.exp( intercept.value), '+-',np.abs(np.exp( intercept.value )*intercept.stderr ))
print()


def lin(x,n,m):
    return n*x + m



############################################################################################################
####        Ploteo
############################################################################################################

filename = args.data.split('/')[-1].split('.')[0]
plt.style.use(['seaborn-talk','seaborn-whitegrid'])
fig,subplot = plt.subplots(ncols=1,nrows=1)



color = args.data_color
subplot.plot(log_bin_centers,np.log(non_essential_probability_by_degree),'o',
        color=color, label = 'data')

#fit plot
xmin,xmax = subplot.get_xlim()   # para dibujar una linea de extremo a extremo
x = np.linspace(xmin,xmax,10)
y = lin(x,result.params['slope'].value,result.params['intercept'].value)
subplot.plot(x,y,'k--',label='fit')

subplot.set_xlabel('$k$',fontsize=20)
subplot.set_ylabel('$\log(1-P_E)$',fontsize=20)
subplot.tick_params(labelsize=20)
plt.savefig('./schemes/'+filename+'.pdf')
plt.show()

