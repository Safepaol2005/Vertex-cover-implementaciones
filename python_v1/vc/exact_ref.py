"""Óptimo de referencia con OR-Tools CP-SAT (solo para medir la calidad, no compite)."""
from __future__ import annotations

from typing import Optional, Set, Tuple

from .graph import Adj, edges_of


def cpsat_vertex_cover(adj: Adj, time_limit: float = 60.0, workers: int = 1
                       ) -> Tuple[Optional[Set[int]], bool]:
    from ortools.sat.python import cp_model

    m = cp_model.CpModel()
    x = {v: m.NewBoolVar(f"x{v}") for v in adj}
    for u, v in edges_of(adj):
        m.AddBoolOr([x[u], x[v]])
    m.Minimize(sum(x.values()))
    s = cp_model.CpSolver()
    s.parameters.max_time_in_seconds = time_limit
    s.parameters.num_workers = workers
    st = s.Solve(m)
    if st not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        return None, False
    return {v for v in adj if s.Value(x[v])}, st == cp_model.OPTIMAL


def brute_force(adj: Adj) -> Set[int]:
    """Solo para tests (n <= 18)."""
    import itertools

    nodes = list(adj)
    E = edges_of(adj)
    for k in range(len(nodes) + 1):
        for S in itertools.combinations(nodes, k):
            s = set(S)
            if all(u in s or v in s for u, v in E):
                return s
    return set(nodes)
