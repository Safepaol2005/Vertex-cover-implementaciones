"""Algoritmos probabilísticos para Vertex Cover.

- randomized_greedy (Gao, Friedrich, Neumann & Hercher, GECCO 2018): en cada paso elige
  un vértice u del grafo residual con probabilidad deg(u)^γ / Σ deg(v)^γ.
  Con γ = 1 la razón de aproximación esperada es 2 (Teorema 3.1); repetir el algoritmo
  R veces y quedarse con el mejor es un esquema Monte Carlo: la cubierta siempre es válida,
  pero alcanzar una α-aproximación (α < 2) solo ocurre con cierta probabilidad.
- random_order_mis (Veldt, arXiv:2209.04673): recorre los vértices en orden aleatorio
  uniforme y construye un conjunto independiente maximal I; V \\ I es una cubierta
  2-aproximada en esperanza. O(|E|).
"""
from __future__ import annotations

import random
from typing import Set

import numpy as np

from .graph import Adj, remove_redundant


def randomized_greedy_once(adj: Adj, gamma: float, rng: np.random.Generator) -> Set[int]:
    nodes = list(adj)
    idx = {v: i for i, v in enumerate(nodes)}
    nbrs = [[idx[u] for u in adj[v]] for v in nodes]
    deg = np.array([len(x) for x in nbrs], dtype=float)
    alive = [set(x) for x in nbrs]
    C: Set[int] = set()
    remaining = int(deg.sum()) // 2
    while remaining > 0:
        w = np.power(deg, gamma)
        w[deg == 0] = 0.0
        i = int(rng.choice(len(nodes), p=w / w.sum()))
        C.add(nodes[i])
        for j in alive[i]:
            alive[j].discard(i)
            deg[j] -= 1
        remaining -= len(alive[i])
        alive[i].clear()
        deg[i] = 0
    return C


def randomized_greedy(adj: Adj, gamma: float = 1.0, runs: int = 1, seed: int = 0,
                      prune: bool = False) -> Set[int]:
    """Monte Carlo: R ejecuciones independientes, se devuelve la cubierta más pequeña."""
    rng = np.random.default_rng(seed)
    best: Set[int] | None = None
    for _ in range(runs):
        C = randomized_greedy_once(adj, gamma, rng)
        if prune:
            C = remove_redundant(adj, C)
        if best is None or len(C) < len(best):
            best = C
    return best if best is not None else set()


def random_order_mis(adj: Adj, seed: int = 0) -> Set[int]:
    rng = random.Random(seed)
    order = list(adj)
    rng.shuffle(order)
    I: Set[int] = set()
    blocked: Set[int] = set()
    for v in order:
        if v not in blocked:
            I.add(v)
            blocked.add(v)
            blocked |= adj[v]
    return set(adj) - I
