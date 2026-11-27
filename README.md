# Vertex Cover: voraz, programación dinámica, aleatorizado y el híbrido HyVC (C++)

Proyecto de *Topics on Advanced Algorithms* (UNAL) — Samuel Felipe Palacios Olaya y Juan Diego Díaz Ávila.

| Archivo | Algoritmo | Artículo |
|---|---|---|
| `src/greedy_wh.hpp` | Voraz de grado máximo, O(n+m) | Weigt & Hartmann, PRL 84 (2000); PRE 63 (2001) |
| `src/dp_td.hpp` | Programación dinámica sobre descomposición en árbol, O(2^w w n) | Bannach & Berndt, Algorithms 12(8) (2019); Alber, Dorn & Niedermeier, DAM 145 (2005) |
| `src/random_mis.hpp` | MIS en orden aleatorio, 2-aprox. en esperanza, O(n+m) | Veldt, SOSA 2024 |
| `src/divide_conquer.hpp` | Divide y vencerás exacto (punto de comparación exponencial) | Weigt & Hartmann, PRE 63 (2001) |
| `src/reductions.hpp`, `src/bounds.hpp`, `src/hybrid.hpp` | **HyVC**: reducciones + componentes + PD + cotas + ramificación | este trabajo |

## Compilar y probar

```bash
make                 # bin/vcsolve, bin/gen, bin/test_bruteforce
make test            # 2000 grafos aleatorios contra fuerza bruta
```

## Ejecutar

```bash
./bin/vcsolve hybrid data/pace2019/.../vc-exact_001.gr --time 10
./bin/vcsolve greedy  grafo.gr
./bin/vcsolve rmis    grafo.gr --runs 64 --seed 1
./bin/vcsolve dp      grafo.gr --width 22
./bin/vcsolve dc      grafo.gr --time 10
```
Salida (CSV): `alg,archivo,n,m,tamaño,tiempo_s,estado,extra`, con estado `ok` (óptimo si el algoritmo es exacto), `heur`, `anytime`, `timeout` o `toowide`.
Formatos aceptados: PACE (`p td n m`), DIMACS (`p edge n m` / `e u v`) y lista de aristas.

## Datos (públicos)

- PACE 2019 Vertex Cover (200 instancias): https://zenodo.org/records/3368306 — descomprimir en `data/pace2019/`.
- PACE 2017 treewidth exact (200 instancias, CC0): https://github.com/PACE-challenge/Treewidth-PACE-2017-instances — `.gr` en `data/pace2017/`.
- Sintéticas: `bash bench/make_synthetic.sh`.

## Reproducir los experimentos

```bash
bash bench/run_all.sh                       # todos los algoritmos, 10 s por instancia
python3 bench/cpsat_ref.py results/<ds>_need_ref.txt results/<ds>_cpsat.csv 60 2   # óptimos faltantes
python3 bench/make_tables.py                # tablas y figuras -> paper/
cd paper && latexmk -pdf main.tex
```

## Versión anterior

La primera versión del proyecto (Python, con greedy/PD/divide y vencerás/aleatorizado de Gao et al.) está en `python_v1/`.
