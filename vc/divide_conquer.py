"""Divide y vencerás exacto (Weigt & Hartmann, Phys. Rev. E 63, 056127, 2001).

En cada nodo:
  1. toma el vértice libre i de mayor grado actual;
  2. rama izquierda: i entra en la cubierta (se borra i);
     rama derecha: i queda fuera, así que TODOS sus vecinos entran (se borra N[i]);
  3. tras cada rama, el grafo libre se parte en componentes conexas que se resuelven
     de forma independiente y cuyos tamaños se suman (el "divide").
  4. se devuelve la rama más pequeña.
Debajo del umbral de percolación (grado medio c < 1) las componentes son O(log N) y el
algoritmo es rápido; por encima explota.
"""
from __future__ import annotations

from typing import Set

from .graph import Adj, Deadline, connected_components, induced


def _solve(adj: Adj, dl: Deadline) -> Set[int]:
    dl.check()
    live = {v: n for v, n in adj.items() if n}
    if not live:
        return set()
    comps = connected_components(live)
    if len(comps) > 1:
        out: Set[int] = set()
        for comp in comps:
            out |= _solve(induced(live, comp), dl)
        return out
    i = max(live, key=lambda x: (len(live[x]), -x))
    # rama 1: i en la cubierta
    g1 = {v: n - {i} for v, n in live.items() if v != i}
    c1 = _solve(g1, dl) | {i}
    # rama 2: i fuera => N(i) dentro
    Ni = live[i]
    drop = Ni | {i}
    g2 = {v: n - drop for v, n in live.items() if v not in drop}
    c2 = _solve(g2, dl) | Ni
    return c1 if len(c1) <= len(c2) else c2


def divide_and_conquer(adj: Adj, deadline: Deadline | None = None) -> Set[int]:
    return _solve(adj, deadline or Deadline(None))
