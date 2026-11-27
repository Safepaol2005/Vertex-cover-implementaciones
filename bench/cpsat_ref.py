"""Óptimo de referencia con OR-Tools CP-SAT para instancias que ningún algoritmo exacto nuestro certificó.

Uso: python3 bench/cpsat_ref.py lista.txt salida.csv [segundos] [hilos]
lista.txt: una ruta de grafo por línea. Salida: instancia,cpsat_size,cpsat_optimal,cpsat_bound,time_s
"""
import csv
import sys
import time

from ortools.sat.python import cp_model


def read_graph(path):
    edges, n = [], None
    with open(path) as f:
        for line in f:
            if not line.strip() or line[0] in "c%#":
                continue
            p = line.split()
            if p[0] == "p":
                n = int(p[2])
                continue
            if p[0] == "e":
                p = p[1:]
            u, v = int(p[0]) - 1, int(p[1]) - 1
            if u != v:
                edges.append((min(u, v), max(u, v)))
    return n, sorted(set(edges))


def solve(path, limit, workers):
    n, E = read_graph(path)
    m = cp_model.CpModel()
    x = [m.NewBoolVar("") for _ in range(n)]
    for u, v in E:
        m.AddBoolOr([x[u], x[v]])
    m.Minimize(sum(x))
    s = cp_model.CpSolver()
    s.parameters.max_time_in_seconds = limit
    s.parameters.num_workers = workers
    t = time.time()
    st = s.Solve(m)
    ok = st in (cp_model.OPTIMAL, cp_model.FEASIBLE)
    return (int(s.ObjectiveValue()) if ok else "", st == cp_model.OPTIMAL,
            int(s.BestObjectiveBound()) if ok else "", round(time.time() - t, 2))


if __name__ == "__main__":
    paths = [l.strip() for l in open(sys.argv[1]) if l.strip()]
    limit = float(sys.argv[3]) if len(sys.argv) > 3 else 60
    workers = int(sys.argv[4]) if len(sys.argv) > 4 else 2
    with open(sys.argv[2], "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["file", "cpsat_size", "cpsat_optimal", "cpsat_bound", "time_s"])
        for p in paths:
            w.writerow([p, *solve(p, limit, workers)])
            f.flush()
            print(p, flush=True)
