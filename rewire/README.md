# run

---
Para el alpha, simulamos muchas realizaciones de la red, haciendo
un recableado de enlaces, manteniendo su distribucion de grados Pk.
Esto genera una secuencia de los valores del nro de interacciones esenciales,
y los guarda en un archivo `./test.h5`.
```bash
./alpha -- --bins 50 --n_rewire 6000 -fig ./LIT.png -inp ../data/yeast_LIT.txt -out ./LIT.h5
```
Notar q si se usa `n_rewire` cerca a 10000, el histograma se hace raro. Tal vez
se deba a q la generacion de nro aleatorio entra en un periodo de ciclo?


---
Hallar los parametros de fiteo de la distribucion del nro de interacc esenciales.
```bash
./fitear.py -- -inp ./LIT.h5  ## anexa la informacion en el archivo `./test.h5`
```


---
Simulacion de redes, asumiendo q `N_ie-n_ie` interacciones se debe a causas 
random, y que `n_ie` se debe a IBEPs; siendo `n_ie` un nro entero aleatorio q 
se deduce de la distribucion mencionada antes, y `N_ie` la cantidad de interacciones
esenciales en la red real.
A partir de estas realizaciones de red, se deduce un beta promedio (y su error), asi
como tambien la fraccion de overlapping (por intentos de asignar esencialidad a nodos
q ya eran esenciales).
```bash
./beta.py -- -inp_txt ../data/yeast_LIT.txt -inp_h5 ./LIT.h5
```

---
Para hacer las figuras de distribuciones y de sus ajustes:
```bash
./histos.py -- -inp_txt ../data/yeast_LIT.txt -inp_h5 LIT.h5 -fig hist_LIT.png
```
