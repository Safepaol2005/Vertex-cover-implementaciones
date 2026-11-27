# HyVC — Vertex Cover híbrido (proyecto Topics on Advanced Algorithms, UNAL)

Samuel Felipe Palacios Olaya · Juan Diego Díaz Ávila

Implementaciones en Python de cuatro paradigmas para *Minimum Vertex Cover*, tomadas de
artículos publicados, y de **HyVC**, un híbrido construido comparando sus secciones.

| Archivo | Algoritmo | Fuente |
|---|---|---|
| `vc/greedy.py` | `max_degree_greedy`, `edge_greedy_fastvc`, `matching_2approx` | Cai (IJCAI 2015); Cai, Lin & Luo (JAIR 2017); Weigt & Hartmann (PRE 2001) |
| `vc/dp_treedec.py` | `dp_vertex_cover` — PD sobre descomposición en árbol, O*(2^tw) | Bannach & Berndt (ESA 2018 / Algorithms 2019); Alber, Dorn & Niedermeier (DAM 2005) |
| `vc/divide_conquer.py` | `divide_and_conquer` — ramificación + componentes conexas | Weigt & Hartmann (PRE 2001) |
| `vc/randomized.py` | `randomized_greedy` (deg^γ, Monte Carlo), `random_order_mis` | Gao, Friedrich, Neumann & Hercher (GECCO 2018); Veldt (arXiv 2022) |
| `vc/reductions.py` | Reglas R0–R4 (grado 0/1/2, folding, dominancia) | Ásgeirsson & Stein (SIDMA 2009); Akiba & Iwata (TCS 2016) |
| `vc/hybrid.py` | **HyVC** | este trabajo |
| `vc/exact_ref.py` | CP-SAT (referencia) y fuerza bruta (tests) | OR-Tools |

## Uso

```bash
pip install numpy matplotlib ortools        # ortools solo para la referencia y los tests
python tests/test_correctness.py            # vs. fuerza bruta (180 grafos)
python tests/test_vs_cpsat.py               # exactos vs. CP-SAT (33 grafos)
python bench/run_experiments.py --out results/main.csv           # benchmark principal (~45 min)
python bench/run_experiments.py --ablation --out results/ablation.csv
python bench/run_experiments.py --gamma --out results/gamma.csv
python bench/make_figures.py                # tablas y figuras -> paper/tables, paper/figs
cd paper && latexmk -pdf main.tex
```

```python
from vc import hybrid_vc
from vc.graph import erdos_renyi
g = erdos_renyi(200, 4.0, seed=1)      # dict[int, set[int]]
cover, info = hybrid_vc(g, time_limit=10)
print(len(cover), info["exact"])       # exact=True => óptimo certificado
```
