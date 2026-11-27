"""HyVC: algoritmo híbrido para Minimum Vertex Cover.

Cada fase toma la "sección" que mejor funcionó en la comparación de los artículos:

  Fase 1  REDUCIR      reglas exactas R0-R4 (Ásgeirsson & Stein 2009; Akiba & Iwata 2016).
  Fase 2  DIVIDIR      componentes conexas resueltas de forma independiente
                       (Weigt & Hartmann 2001), también dentro de cada nodo de la búsqueda.
  Fase 3  DESPACHAR    si el costo estimado de la PD (suma de 2^|bolsa| en la descomposición
                       por grado mínimo) es <= cost_cap, la componente se resuelve EXACTA por
                       programación dinámica sobre la descomposición
                       (Bannach & Berndt 2018/19; Alber, Dorn & Niedermeier 2005).
  Fase 4  ACOTAR       cota superior inicial Monte Carlo: R corridas del greedy aleatorizado
                       deg^γ (Gao et al. 2018) + la construcción de FastVC (Cai et al.),
                       ambas con eliminación de redundantes. Cota inferior: matching maximal.
  Fase 5  RAMIFICAR    branch & reduce sobre el vértice de grado máximo (la regla que mejor
                       funciona según Akiba & Iwata) con poda UB/LB; en cada nodo se vuelve a
                       reducir, dividir y despachar a DP.
  Fase 6  ANYTIME      si se agota el tiempo en una componente, se devuelve la mejor cubierta
                       encontrada y se reporta la brecha certificada |C| - LB.

Garantías: la salida SIEMPRE es una cubierta válida; si `info["exact"]` es True es óptima.
"""
from __future__ import annotations

from typing import Dict, Optional, Set

from .dp_treedec import dp_cost, dp_vertex_cover
from .graph import (Adj, Deadline, Timeout, connected_components, copy_adj, induced,
                    lower_bound, remove_redundant)
from .greedy import edge_greedy_fastvc, max_degree_greedy
from .randomized import randomized_greedy
from .reductions import reduce_graph, unfold


class _Ctx:
    def __init__(self, deadline: Deadline, cost_cap: int, use_dp: bool, use_split: bool,
                 use_reduce: bool):
        self.dl = deadline
        self.cost_cap = cost_cap
        self.use_dp = use_dp
        self.use_split = use_split
        self.use_reduce = use_reduce
        self.nodes = 0
        self.dp_calls = 0


MAX_W = 22
INNER_CAP = 1 << 16
BIG_COMP = 100


def _dp_ok(adj: Adj, cost_cap: int) -> bool:
    c = dp_cost(adj, MAX_W)
    return c is not None and c <= cost_cap


def _cheap_ub(adj: Adj) -> Set[int]:
    a = remove_redundant(adj, max_degree_greedy(adj))
    b = edge_greedy_fastvc(adj)
    return a if len(a) <= len(b) else b


def _solve(adj: Adj, ub: int, ctx: _Ctx) -> Optional[Set[int]]:
    """Cubierta mínima de `adj` si su tamaño es < ub; None si no existe ninguna < ub."""
    ctx.dl.check()
    ctx.nodes += 1
    g = copy_adj(adj)
    if ctx.use_reduce:
        g, forced, folds = reduce_graph(g)
    else:
        g, forced, folds = {v: n for v, n in g.items() if n}, set(), []
    off = len(forced) + len(folds)
    if off >= ub:
        return None
    if not any(g.values()):
        return unfold(set(), forced, folds)
    budget = ub - off
    comps = [c for c in connected_components(g) if len(c) > 1] if ctx.use_split else [set(g)]
    subs = [induced(g, c) for c in comps]
    lbs = [lower_bound(s) for s in subs]
    if sum(lbs) >= budget:
        return None
    order = sorted(range(len(subs)), key=lambda i: len(subs[i]))
    total: Set[int] = set()
    used = 0
    rest_lb = sum(lbs)
    for i in order:
        rest_lb -= lbs[i]
        cb = budget - used - rest_lb  # presupuesto exclusivo para esta componente
        sol = _solve_component(subs[i], cb, lbs[i], ctx)
        if sol is None:
            return None
        total |= sol
        used += len(sol)
    return unfold(total, forced, folds)


def _solve_component(adj: Adj, ub: int, lb: int, ctx: _Ctx) -> Optional[Set[int]]:
    if lb >= ub:
        return None
    # Despacho a programación dinámica si el treewidth estimado es pequeño
    # Dentro de la búsqueda la PD puede invocarse miles de veces. En componentes pequeñas la
    # ramificación es barata, así que solo se acepta una PD cara (> INNER_CAP celdas) si la
    # componente es grande (> BIG_COMP vértices), donde ramificar no tiene esperanza.
    cap = ctx.cost_cap if len(adj) > BIG_COMP else min(ctx.cost_cap, INNER_CAP)
    if ctx.use_dp and _dp_ok(adj, cap):
        ctx.dp_calls += 1
        sol = dp_vertex_cover(adj, max_width=MAX_W, deadline=ctx.dl)
        return sol if sol is not None and len(sol) < ub else None
    best: Optional[Set[int]] = None
    h = _cheap_ub(adj)
    if len(h) < ub:
        best, ub = h, len(h)
    if lb >= ub:
        return best
    v = max(adj, key=lambda x: (len(adj[x]), -x))
    # rama 1: v en la cubierta
    g1 = {x: n - {v} for x, n in adj.items() if x != v}
    s1 = _solve(g1, ub - 1, ctx)
    if s1 is not None:
        best, ub = s1 | {v}, len(s1) + 1
    # rama 2: v fuera => N(v) en la cubierta
    Nv = adj[v]
    if len(Nv) < ub:
        drop = Nv | {v}
        g2 = {x: n - drop for x, n in adj.items() if x not in drop}
        s2 = _solve(g2, ub - len(Nv), ctx)
        if s2 is not None:
            best = s2 | Nv
    return best


def hybrid_vc(adj: Adj, time_limit: float | None = 10.0, cost_cap: int = 1 << 23, runs: int = 16,
              gamma: float = 2.0, seed: int = 0, use_dp: bool = True,
              use_split: bool = True, use_reduce: bool = True, use_mc: bool = True
              ) -> tuple[Set[int], Dict]:
    dl = Deadline(time_limit)
    ctx = _Ctx(dl, cost_cap, use_dp, use_split, use_reduce)
    g = copy_adj(adj)
    if use_reduce:
        g, forced, folds = reduce_graph(g)
    else:
        g, forced, folds = {v: n for v, n in g.items() if n}, set(), []
    comps = [c for c in connected_components(g) if len(c) > 1] if use_split else ([set(g)] if g else [])
    cover: Set[int] = set()
    exact = True
    lb_total = len(forced) + len(folds)
    stats = {"components": len(comps), "dp_components": 0, "bb_components": 0,
             "timeouts": 0, "reduced_n": len(g), "n": len(adj)}
    for comp in sorted(comps, key=len):
        sub = induced(g, comp)
        lb = lower_bound(sub)
        lb_total += lb
        if use_dp and _dp_ok(sub, cost_cap):
            sol = dp_vertex_cover(sub, max_width=MAX_W)
            stats["dp_components"] += 1
            cover |= sol
            continue
        # cota superior Monte Carlo + FastVC
        ub_sol = _cheap_ub(sub)
        if use_mc:
            mc = randomized_greedy(sub, gamma=gamma, runs=runs, seed=seed, prune=True)
            if len(mc) < len(ub_sol):
                ub_sol = mc
        stats["bb_components"] += 1
        if dl.expired():
            exact = False
            stats["timeouts"] += 1
            cover |= ub_sol
            continue
        try:
            better = _solve_component(sub, len(ub_sol), lb, ctx)
            cover |= better if better is not None else ub_sol  # None => ub_sol es óptima
        except Timeout:
            exact = False
            stats["timeouts"] += 1
            cover |= ub_sol
    cover = unfold(cover, forced, folds)
    if not exact:
        cover = remove_redundant(adj, cover)
    stats.update(exact=exact, lower_bound=lb_total, nodes=ctx.nodes, inner_dp=ctx.dp_calls)
    return cover, stats
