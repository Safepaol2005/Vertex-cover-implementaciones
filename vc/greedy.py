"""Algoritmos voraces (greedy) para Vertex Cover.

- max_degree_greedy: regla clásica "toma el vértice de mayor grado residual"
  (heurística usada por Weigt & Hartmann 2001 como cota superior; Θ(log n)-aprox. en el peor caso).
- edge_greedy_fastvc: construcción de FastVC (Cai 2015; Cai, Lin & Luo 2017, JAIR):
  recorre las aristas, si una está descubierta toma el extremo de mayor grado, y después
  elimina vértices redundantes. O(|E|).
- matching_2approx: 2-aproximación por matching maximal (referencia teórica).
"""
from __future__ import annotations

import heapq
from typing import Set

from .graph import Adj, greedy_maximal_matching, remove_redundant


def max_degree_greedy(adj: Adj) -> Set[int]:
    deg = {v: len(n) for v, n in adj.items()}
    heap = [(-d, v) for v, d in deg.items() if d > 0]
    heapq.heapify(heap)
    alive = {v: set(n) for v, n in adj.items()}
    C: Set[int] = set()
    while heap:
        d, v = heapq.heappop(heap)
        if v in C or -d != deg[v] or deg[v] == 0:
            continue  # entrada obsoleta (lazy deletion)
        C.add(v)
        for u in alive[v]:
            alive[u].discard(v)
            deg[u] -= 1
            if deg[u] > 0:
                heapq.heappush(heap, (-deg[u], u))
        alive[v].clear()
        deg[v] = 0
    return C


def edge_greedy_fastvc(adj: Adj) -> Set[int]:
    C: Set[int] = set()
    deg = {v: len(n) for v, n in adj.items()}
    for u in adj:
        for v in adj[u]:
            if u < v and u not in C and v not in C:
                C.add(u if deg[u] >= deg[v] else v)
    # "read-one": FastVC elimina los vértices con pérdida 0 (todos sus vecinos cubiertos)
    return remove_redundant(adj, C, order=sorted(C, key=lambda x: deg[x]))


def matching_2approx(adj: Adj) -> Set[int]:
    C: Set[int] = set()
    for u, v in greedy_maximal_matching(adj):
        C.add(u)
        C.add(v)
    return C
