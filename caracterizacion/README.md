
# ploteo de grafos
crea un plot del grafo a partir de un set de datos
ej:
```bash
 plot.py ../data/yeast_LIT_Reguly.gml -v 'bright blue' -s 5 -lw 0.5 
```

# observables
retorna observables estadisticos de una red, en el caso en que la red no este en formato GML, usar opcion ```-f ncol```

```bash
 observables.py ../data/yeast_LIT_Reguly.txt -f ncol
 ```
 
# cobertura y coherencia
## cobertura
 entrega diagrama de venn de las proteinas de los sets de proteinas cubiertas en las redes
 ```bash
 venn2.py red1 red2
 ```
## coherencia
 entrega el diagrama de venn del set de enlaces cubiertos por las proteinas comunes entre 2 redes
 ```bash
 venn2_coherence.py red1 red2
 ```

# regla centralidad-letalidad
retorna grafico de relacion entre esencialidad y hubs, ademas retorna coeficientes de correlacion de Rendall y Spearman dado un umbral de definicion de hub (opcion -p)
```bash
ess_hub.py datos_esenciales red1 red2 red3 ... [-p porcentaje]
```
