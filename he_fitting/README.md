

# Fiteo
Retorna probabilidades y estadisticas del fiteo para el modelo de He
```bash
fitting.py [-h]  [-f formato] [-k grado_de_corte] data_esencial red
```

#Conteo de pares
Retorna numero de pares totales de proteinas no-vecinas que comparten un numero mayor o igual a threshold de vecinos.
Ademas retorna el numero de pares de igual tipo (Esencial-esencial o no-esencial no-esencial)
Si se dan probabilidades alpha y beta retorna tambien la estimacion a partir del modelo de He para el numero de pares de igual tipo
```bash
count_pairs.py [-h] [-f format] [-t threshold] [-a alpha] [-b beta] data_escencial red
```
