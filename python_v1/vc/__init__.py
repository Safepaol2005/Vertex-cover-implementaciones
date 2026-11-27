"""Vertex Cover: implementaciones base (greedy, DP, divide y vencerás, aleatorizado) e híbrido HyVC."""
from .divide_conquer import divide_and_conquer
from .dp_treedec import dp_vertex_cover, treewidth_upper_bound
from .graph import Deadline, Timeout, is_vertex_cover
from .greedy import edge_greedy_fastvc, matching_2approx, max_degree_greedy
from .hybrid import hybrid_vc
from .randomized import random_order_mis, randomized_greedy

__all__ = [
    "divide_and_conquer", "dp_vertex_cover", "treewidth_upper_bound", "Deadline", "Timeout",
    "is_vertex_cover", "edge_greedy_fastvc", "matching_2approx", "max_degree_greedy",
    "hybrid_vc", "random_order_mis", "randomized_greedy",
]
