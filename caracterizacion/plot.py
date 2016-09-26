#!/bin/env/python3

from graph_tool.all import * 
import seaborn as sns
import argparse as arg
import sys

# seccion de argumentos del programa
AcceptedFormats = [ 'auto','gml','gt','graphml','xml','dot','csv' ]     # formatos soportados
vcolors = list(sns.xkcd_rgb.keys()) # nombres de colores (http://xkcd.com/color/rgb/)

argparser = arg.ArgumentParser(description='')
argparser.add_argument('file',help='graph file')
argparser.add_argument('--format','-f', default='auto',
        help='force format of the file (default:autodetect)',choices=AcceptedFormats)
argparser.add_argument('--is_directed','-d',action='store_true',help='directed graph')
argparser.add_argument('--vertex-color','-v',help='set the vertex color (xkcd palette http://xkcd.com/color/rgb/)',
        default='red')

argparser.add_argument('--vertex-size','-s',help='set the vertex size (default:5)',default=5,type=float)
argparser.add_argument('--edge-width','-lw',help='set the edge size (default:1.2)',default=1.2,type=float)
argparser.add_argument('--gamma',help='set the Strength of the attractive force between connected components (default:1)',default=1,type=float)


args = argparser.parse_args()
if args.vertex_color not in vcolors: #si no conocemos el color tirar error !!
    print('No conocemos el color seleccionado, por favor mira la lista en http://xkcd.com/color/rgb/')
    sys.exit(1)  

#importa el grafo
if args.format == 'csv':   # si el formato es una lista de links hay que tratarlo distinto (csv)
    graph = load_graph_from_csv(args.file,string_vals=True,directed=args.is_directed,
                csv_options={"delimiter": "\t", "quotechar": "#"})
else:  # si no es csv, que lo lea tranqui..
    graph = load_graph(args.file,fmt=args.format)

#chequea que todo sea como "debe ser"
print('chequeando...',args.file)
print('vertices',len(list(graph.vertices()))) #numero de vertices
print('edges',len(list(graph.edges())))       #numero de links


# seteo de algunas argumentos de la creacion del grafo
filename = '-'.join(args.file.split('/')[-1].split('.'))    #extraccion de simbolos raros en el nombre del archivo
pos = sfdp_layout(graph , gamma=args.gamma)                                    # layout que nos gusta mas
cmap = sns.cubehelix_palette(8, start=.5, rot=-.75, as_cmap=True)  # estaba en el ejemplo de guia...
vcolor = sns.xkcd_rgb[args.vertex_color]  # color que elegimos para nuestro grafo (por defecto rojo)


#dibuja el grafico en el archivo filename.pdf
graph_draw(graph, pos, output_size=(1000, 1000), vertex_color=[0,0,0,.6],
                   vertex_fill_color= vcolor,
                   vertex_pen_width=0.8,
                   vertex_size=args.vertex_size, edge_pen_width=args.edge_width,
                   vcmap=cmap, 
                   output=filename+'.pdf')




