"""Reglas de reducción exactas (preservan el óptimo).

Tomadas de la fase de "reducciones óptimas" de Ásgeirsson & Stein (SIAM J. Discrete Math.
23(3), 2009) y del conjunto de reglas de Akiba & Iwata (TCS 609, 2016):
  R0 grado 0 : se borra el vértice.
  R1 grado 1 : su único vecino entra en la cubierta.
  R2 grado 2 con vecinos adyacentes (triángulo): ambos vecinos entran.
  R3 grado 2 con vecinos no adyacentes (folding): v, u, w se contraen en un vértice x;
               OPT(G) = OPT(G') + 1. Al deshacer: x ∈ C' ⇒ {u, w}; x ∉ C' ⇒ {v}.
  R4 dominancia: si u ∈ N(v) y N[u] ⊆ N[v], existe una cubierta óptima que contiene a v.
"""
from __future__ import annotations

import itertools
from typing import List, Set, Tuple

from .graph import Adj

Fold = Tuple[int, int, int, int]  # (x, v, u, w)

_new_ids = itertools.count(10 ** 9)


def reduce_graph(adj: Adj, use_fold: bool = True, use_dominance: bool = True
                 ) -> Tuple[Adj, Set[int], List[Fold]]:
    """Trabaja sobre `adj` (se modifica). Devuelve (grafo reducido, forzados, folds)."""
    forced: Set[int] = set()
    folds: List[Fold] = []

    def take(v: int) -> None:
        forced.add(v)
        for a in adj.pop(v):
            adj[a].discard(v)
            queue.add(a)

    def drop(v: int) -> None:
        for a in adj.pop(v):
            adj[a].discard(v)
            queue.add(a)

    queue: Set[int] = set(adj)
    changed = True
    while changed:
        changed = False
        while queue:
            v = queue.pop()
            if v not in adj:
                continue
            d = len(adj[v])
            if d == 0:
                del adj[v]
            elif d == 1:
                (u,) = adj[v]
                take(u)
                changed = True
            elif d == 2:
                u, w = adj[v]
                if w in adj[u]:
                    take(u)
                    take(w)
                    changed = True
                elif use_fold:
                    x = next(_new_ids)
                    nx = (adj[u] | adj[w]) - {v, u, w}
                    for z in (v, u, w):
                        for a in adj.pop(z):
                            adj.get(a, set()).discard(z)
                    adj[x] = set(nx)
                    for a in nx:
                        adj[a].add(x)
                        queue.add(a)
                    queue.add(x)
                    folds.append((x, v, u, w))
                    changed = True
        if use_dominance:
            for v in sorted(adj, key=lambda z: -len(adj[z])):
                if v not in adj:
                    continue
                Nv = adj[v]
                dv = len(Nv)
                for u in Nv:
                    if len(adj[u]) <= dv and adj[u] - {v} <= Nv:
                        take(v)
                        changed = True
                        break
    return adj, forced, folds


def unfold(cover: Set[int], forced: Set[int], folds: List[Fold]) -> Set[int]:
    C = set(cover) | forced
    for x, v, u, w in reversed(folds):
        if x in C:
            C.discard(x)
            C.add(u)
            C.add(w)
        else:
            C.add(v)
    return C
