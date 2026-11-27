"""Programación dinámica sobre descomposición en árbol (exacto, O*(2^tw)).

Esquema "folklore" descrito por Bannach & Berndt (ESA 2018 / Algorithms 2019) y evaluado
experimentalmente para VC planar por Alber, Dorn & Niedermeier (DAM 2005):
  - estado = subconjunto S de la bolsa (bitmask) que pertenece a la cubierta,
  - un estado es válido si cubre todas las aristas ORIGINALES dentro de la bolsa,
  - "forget" = mínimo sobre el bit del vértice olvidado,
  - "join"   = suma de las contribuciones de los hijos sobre la intersección.

La descomposición se obtiene con la heurística de eliminación por grado mínimo:
bag(v) = {v} ∪ N(v) en el momento de eliminar v; el padre de bag(v) es el vecino de v que
se elimina primero después de v. Con esta construcción bag(hijo) \\ {hijo} ⊆ bag(padre),
así que cada hijo aporta una tabla ya "olvidada" indexada por bits de la bolsa del padre.
"""
from __future__ import annotations

import heapq
from typing import Dict, List, Optional, Set, Tuple

import numpy as np

from .graph import Adj, Deadline

INF = np.int64(1 << 40)


def min_degree_elimination(adj: Adj, max_width: int | None = None
                           ) -> Optional[Tuple[List[int], Dict[int, List[int]], int]]:
    """Devuelve (orden, bolsas, ancho). bolsas[v] = [vecinos posteriores..., v] (v al final).
    Si el ancho supera max_width aborta y devuelve None."""
    g = {v: set(n) for v, n in adj.items()}
    heap = [(len(n), v) for v, n in g.items()]
    heapq.heapify(heap)
    order: List[int] = []
    bags: Dict[int, List[int]] = {}
    width = 0
    done: Set[int] = set()
    while heap:
        d, v = heapq.heappop(heap)
        if v in done or d != len(g[v]):
            continue
        N = list(g[v])
        width = max(width, len(N))
        if max_width is not None and width > max_width:
            return None
        bags[v] = N + [v]
        order.append(v)
        done.add(v)
        for a in N:
            g[a].discard(v)
        for i, a in enumerate(N):  # fill-in: N(v) se vuelve clique
            for b in N[i + 1:]:
                if b not in g[a]:
                    g[a].add(b)
                    g[b].add(a)
        for a in N:
            heapq.heappush(heap, (len(g[a]), a))
        del g[v]
    return order, bags, width


def dp_cost(adj: Adj, max_width: int = 22) -> Optional[int]:
    """Costo estimado de la PD = suma de 2^|bolsa| (celdas de tabla). None si el ancho > max_width."""
    r = min_degree_elimination(adj, max_width)
    if r is None:
        return None
    return sum(1 << len(B) for B in r[1].values())


def treewidth_upper_bound(adj: Adj, cap: int | None = None) -> int:
    r = min_degree_elimination(adj, cap)
    return (cap + 1) if r is None else r[2]


_POP = None


def _popcount_table(bits: int) -> np.ndarray:
    global _POP
    if _POP is None or len(_POP) < (1 << bits):
        n = 1 << max(bits, 1)
        p = np.zeros(n, dtype=np.int64)
        for b in range(max(bits, 1)):
            p[1 << b:1 << (b + 1)] = p[:1 << b] + 1
        _POP = p
    return _POP[: 1 << bits]


def dp_vertex_cover(adj: Adj, max_width: int = 16, deadline: Deadline | None = None
                    ) -> Optional[Set[int]]:
    """Cubierta mínima exacta si la heurística encuentra ancho <= max_width; si no, None."""
    if not adj:
        return set()
    r = min_degree_elimination(adj, max_width)
    if r is None:
        return None
    order, bags, _ = r
    pos = {v: i for i, v in enumerate(order)}
    parent: Dict[int, Optional[int]] = {}
    children: Dict[int, List[int]] = {v: [] for v in order}
    for v in order:
        later = bags[v][:-1]
        if later:
            p = min(later, key=lambda x: pos[x])
            parent[v] = p
            children[p].append(v)
        else:
            parent[v] = None

    forgotten: Dict[int, np.ndarray] = {}   # g_c[m]: tabla del hijo con su vértice olvidado
    argbit: Dict[int, np.ndarray] = {}      # bit elegido del vértice olvidado (reconstrucción)
    proj_pos: Dict[int, np.ndarray] = {}    # posiciones de bag(c)\{c} dentro de bag(padre)

    for v in order:  # los hijos se eliminan antes que el padre
        if deadline is not None:
            deadline.check()
        B = bags[v]
        b = len(B)
        local = {x: i for i, x in enumerate(B)}
        masks = np.arange(1 << b, dtype=np.int64)
        pop = _popcount_table(b)
        # validez: si x no está en S, todos sus vecinos originales en la bolsa deben estar en S
        valid = np.ones(1 << b, dtype=bool)
        for x in B:
            i = local[x]
            nb = 0
            for y in adj[x]:
                j = local.get(y)
                if j is not None:
                    nb |= 1 << j
            if nb:
                valid &= ((masks >> i) & 1).astype(bool) | ((masks & nb) == nb)
        T = np.where(valid, pop.copy(), INF)
        for c in children[v]:
            Bc = bags[c]
            pp = np.array([local[x] for x in Bc[:-1]], dtype=np.int64)
            proj_pos[c] = pp
            idx = np.zeros(1 << b, dtype=np.int64)
            for j, p in enumerate(pp):
                idx |= ((masks >> p) & 1) << j
            gc = forgotten[c]
            contrib = gc[idx] - pop[idx]
            T = np.minimum(T + contrib, INF)
        half = 1 << (b - 1)  # v ocupa el bit más alto: olvidarlo = min de las dos mitades
        lo, hi = T[:half], T[half:]
        argbit[v] = hi < lo
        forgotten[v] = np.minimum(lo, hi)
        # Mantener T completo solo si v es raíz (se reconstruye desde ahí)
    # reconstrucción descendente
    cover: Set[int] = set()
    chosen_mask: Dict[int, int] = {}
    for v in reversed(order):
        B = bags[v]
        b = len(B)
        if parent[v] is None:
            m = 0  # la bolsa de una raíz es {v}: la parte "olvidada" es vacía
        else:
            pm = chosen_mask[parent[v]]
            m = 0
            for j, p in enumerate(proj_pos[v]):
                m |= ((pm >> int(p)) & 1) << j
        top = int(argbit[v][m])
        full = m | (top << (b - 1))
        chosen_mask[v] = full
        if top:
            cover.add(v)
    return cover
