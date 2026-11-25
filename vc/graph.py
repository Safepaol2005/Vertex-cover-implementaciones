"""Utilidades de grafos para Vertex Cover.

Representación: lista de adyacencia como dict[int, set[int]] (grafos simples, no dirigidos).
"""
from __future__ import annotations

import random
import time
from typing import Dict, Iterable, List, Set, Tuple

Adj = Dict[int, Set[int]]


class Timeout(Exception):
    """Se lanza cuando un algoritmo excede su presupuesto de tiempo."""


class Deadline:
    """Reloj compartido: los algoritmos exponenciales llaman a check() periódicamente."""

    def __init__(self, seconds: float | None):
        self.t0 = time.perf_counter()
        self.limit = None if seconds is None else self.t0 + seconds
        self._calls = 0

    def check(self) -> None:
        self._calls += 1
        if self.limit is not None and (self._calls & 63) == 0 and time.perf_counter() > self.limit:
            raise Timeout()

    def expired(self) -> bool:
        return self.limit is not None and time.perf_counter() > self.limit

    def remaining(self) -> float:
        return float("inf") if self.limit is None else self.limit - time.perf_counter()


# --------------------------------------------------------------------------- básicos

def copy_adj(adj: Adj) -> Adj:
    return {v: set(n) for v, n in adj.items()}


def from_edges(edges: Iterable[Tuple[int, int]], nodes: Iterable[int] = ()) -> Adj:
    adj: Adj = {v: set() for v in nodes}
    for u, v in edges:
        if u == v:
            continue
        adj.setdefault(u, set()).add(v)
        adj.setdefault(v, set()).add(u)
    return adj


def edges_of(adj: Adj) -> List[Tuple[int, int]]:
    return [(u, v) for u in adj for v in adj[u] if u < v]


def num_edges(adj: Adj) -> int:
    return sum(len(n) for n in adj.values()) // 2


def is_vertex_cover(adj: Adj, cover: Set[int]) -> bool:
    return all(u in cover or v in cover for u, v in edges_of(adj))


def remove_vertex(adj: Adj, v: int) -> None:
    for u in adj.pop(v, ()):
        adj[u].discard(v)


def remove_vertices(adj: Adj, vs: Iterable[int]) -> None:
    for v in list(vs):
        if v in adj:
            remove_vertex(adj, v)


def induced(adj: Adj, vs: Iterable[int]) -> Adj:
    s = set(vs)
    return {v: adj[v] & s for v in s}


def connected_components(adj: Adj) -> List[Set[int]]:
    seen: Set[int] = set()
    comps: List[Set[int]] = []
    for s in adj:
        if s in seen:
            continue
        comp = {s}
        stack = [s]
        seen.add(s)
        while stack:
            x = stack.pop()
            for y in adj[x]:
                if y not in seen:
                    seen.add(y)
                    comp.add(y)
                    stack.append(y)
        comps.append(comp)
    return comps


def greedy_maximal_matching(adj: Adj) -> List[Tuple[int, int]]:
    matched: Set[int] = set()
    M = []
    # visitar primero vértices de grado bajo da matchings más grandes en la práctica
    for u in sorted(adj, key=lambda x: len(adj[x])):
        if u in matched:
            continue
        best = None
        for w in adj[u]:
            if w not in matched and (best is None or len(adj[w]) < len(adj[best])):
                best = w
        if best is not None:
            matched.add(u)
            matched.add(best)
            M.append((u, best))
    return M


def lower_bound(adj: Adj) -> int:
    """Cota inferior barata: max(|matching maximal|, ceil(m / Δ))."""
    m = num_edges(adj)
    if m == 0:
        return 0
    delta = max(len(n) for n in adj.values())
    return max(len(greedy_maximal_matching(adj)), -(-m // delta))


def remove_redundant(adj: Adj, cover: Set[int], order: Iterable[int] | None = None) -> Set[int]:
    """Post-procesamiento tipo FastVC: quita vértices cuyos vecinos ya están todos en la cubierta."""
    C = set(cover)
    seq = list(order) if order is not None else sorted(C, key=lambda x: len(adj[x]))
    for v in seq:
        if v in C and all(u in C for u in adj[v]):
            C.remove(v)
    return C


# --------------------------------------------------------------------------- generadores

def erdos_renyi(n: int, c: float, seed: int) -> Adj:
    """G(n, p) con p = c/(n-1): grado medio c (régimen de Weigt-Hartmann)."""
    rng = random.Random(seed)
    p = c / (n - 1)
    adj: Adj = {v: set() for v in range(n)}
    for u in range(n):
        for v in range(u + 1, n):
            if rng.random() < p:
                adj[u].add(v)
                adj[v].add(u)
    return adj


def barabasi_albert(n: int, m: int, seed: int) -> Adj:
    """Ley de potencias (preferential attachment), similar a redes reales."""
    rng = random.Random(seed)
    adj: Adj = {v: set() for v in range(n)}
    targets = list(range(m))
    repeated: List[int] = []
    for v in range(m, n):
        for t in set(targets):
            adj[v].add(t)
            adj[t].add(v)
        repeated.extend(targets)
        repeated.extend([v] * m)
        targets = set()
        while len(targets) < m:
            targets.add(rng.choice(repeated))
        targets = list(targets)
    return adj


def partial_k_tree(n: int, k: int, keep: float, seed: int) -> Adj:
    """k-árbol aleatorio con aristas conservadas con probabilidad `keep`: treewidth <= k."""
    rng = random.Random(seed)
    adj: Adj = {v: set() for v in range(n)}
    cliques = [tuple(range(k + 1))]
    full_edges = [(u, v) for u in range(k + 1) for v in range(u + 1, k + 1)]
    for v in range(k + 1, n):
        base = rng.choice(cliques)
        drop = rng.randrange(k + 1)
        sep = base[:drop] + base[drop + 1:]
        for u in sep:
            full_edges.append((u, v))
        cliques.append(sep + (v,))
    for u, v in full_edges:
        if rng.random() < keep:
            adj[u].add(v)
            adj[v].add(u)
    return adj


def banded(n: int, b: int, p: float, seed: int) -> Adj:
    """Grafo "banda": arista i-j (|i-j| <= b) con probabilidad p. Pathwidth <= b, no cordal."""
    rng = random.Random(seed)
    adj: Adj = {v: set() for v in range(n)}
    for i in range(n):
        for j in range(i + 1, min(n, i + b + 1)):
            if rng.random() < p:
                adj[i].add(j)
                adj[j].add(i)
    return adj


def greedy_trap(k: int) -> Adj:
    """Instancia bipartita clásica (Johnson 1974) donde el greedy de grado máximo
    obtiene Θ(log k) veces el óptimo. L = {0..k-1} es la cubierta óptima.
    Para cada i = 2..k, el bloque R_i tiene floor(k/i) vértices, cada uno unido a i
    vértices distintos de L. Los vértices de R reciben los ids más pequeños, de modo que
    el desempate por id (el natural) favorece a R: el caso adverso para el greedy."""
    nR = sum(k // i for i in range(2, k + 1))
    L0 = nR  # L = {nR, ..., nR+k-1}
    adj: Adj = {v: set() for v in range(nR + k)}
    r = 0
    for i in range(k, 1, -1):
        for j in range(k // i):
            for t in range(j * i, (j + 1) * i):
                adj[r].add(L0 + t)
                adj[L0 + t].add(r)
            r += 1
    return adj
