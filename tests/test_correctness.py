"""Pruebas de corrección contra fuerza bruta en grafos pequeños.

Ejecutar desde la raíz del proyecto:  python -m pytest -q   (o  python tests/test_correctness.py)
"""
import os
import random
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from vc import (divide_and_conquer, dp_vertex_cover, edge_greedy_fastvc, hybrid_vc,
                is_vertex_cover, matching_2approx, max_degree_greedy, random_order_mis,
                randomized_greedy)
from vc.exact_ref import brute_force
from vc.graph import erdos_renyi, greedy_trap, partial_k_tree
from vc.reductions import reduce_graph, unfold
from vc.graph import copy_adj


def small_graphs(count=150):
    rng = random.Random(7)
    for i in range(count):
        n = rng.randint(1, 14)
        c = rng.uniform(0.5, 6.0)
        yield erdos_renyi(n, c if n > 1 else 0, seed=i) if n > 1 else {0: set()}
    for i in range(30):
        yield partial_k_tree(rng.randint(5, 16), rng.randint(1, 4), 0.8, seed=i)


def test_all():
    for g in small_graphs():
        opt = len(brute_force(g))
        for name, f in [("maxdeg", max_degree_greedy), ("fastvc", edge_greedy_fastvc),
                        ("match", matching_2approx), ("mis", lambda a: random_order_mis(a, 1)),
                        ("rand", lambda a: randomized_greedy(a, 1.0, 4, 1))]:
            C = f(g)
            assert is_vertex_cover(g, C), name
            assert len(C) >= opt, name
        assert len(matching_2approx(g)) <= 2 * opt
        C = dp_vertex_cover(g, max_width=20)
        assert C is not None and is_vertex_cover(g, C) and len(C) == opt, "dp"
        C = divide_and_conquer(g)
        assert is_vertex_cover(g, C) and len(C) == opt, "dc"
        C, info = hybrid_vc(g, time_limit=None)
        assert info["exact"] and is_vertex_cover(g, C) and len(C) == opt, ("hybrid", len(C), opt)
        # híbrido sin DP (obliga a usar branch & reduce + folding)
        C, info = hybrid_vc(g, time_limit=None, use_dp=False)
        assert is_vertex_cover(g, C) and len(C) == opt, ("hybrid-noDP", len(C), opt)
        # reducciones: OPT(G) = OPT(G') + |forced| + |folds|
        gr, forced, folds = reduce_graph(copy_adj(g))
        sub = brute_force(gr)
        assert len(sub) + len(forced) + len(folds) == opt
        assert is_vertex_cover(g, unfold(sub, forced, folds))


def test_trap():
    g = greedy_trap(12)
    assert len(max_degree_greedy(g)) > 12  # el greedy cae en la trampa
    C, info = hybrid_vc(g, time_limit=None)
    assert len(C) == 12


if __name__ == "__main__":
    test_all()
    test_trap()
    print("OK")
