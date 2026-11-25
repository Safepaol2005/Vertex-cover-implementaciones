"""Benchmark: compara las implementaciones base y HyVC.

Uso (desde la raíz del proyecto):
    python bench/run_experiments.py --limit 20 --out results/main.csv
    python bench/run_experiments.py --ablation --out results/ablation.csv
    python bench/run_experiments.py --gamma --out results/gamma.csv
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from vc import (Deadline, Timeout, divide_and_conquer, dp_vertex_cover, edge_greedy_fastvc,
                hybrid_vc, is_vertex_cover, matching_2approx, max_degree_greedy,
                random_order_mis, randomized_greedy, treewidth_upper_bound)
from vc.exact_ref import cpsat_vertex_cover
from vc.graph import (banded, barabasi_albert, erdos_renyi, greedy_trap, num_edges,
                      partial_k_tree)


def instances():
    for c in (1.0, 2.0, 3.0, 4.0, 6.0):
        for s in range(3):
            yield "ER-n200", f"er200_c{c:g}_s{s}", erdos_renyi(200, c, 100 + s), c
    for c in (2.0, 4.0):
        for s in range(2):
            yield "ER-n500", f"er500_c{c:g}_s{s}", erdos_renyi(500, c, 200 + s), c
    for k in (4, 8, 12):
        for s in range(3):
            yield "PartialKTree", f"pkt300_k{k}_s{s}", partial_k_tree(300, k, 0.6, 300 + s), k
    for b in (8, 12, 16):
        for s in range(2):
            yield "Banded", f"band400_b{b}_s{s}", banded(400, b, 0.4, 600 + s), b
    for m in (2, 3, 5):
        for s in range(2):
            yield "BA-n1000", f"ba1000_m{m}_s{s}", barabasi_albert(1000, m, 400 + s), m
    for k in (20, 40, 80):
        yield "GreedyTrap", f"trap_k{k}", greedy_trap(k), k
    for s in range(3):
        yield "ER-dense", f"er60_c15_s{s}", erdos_renyi(60, 15.0, 500 + s), 15.0


def run_alg(name, g, limit):
    t0 = time.perf_counter()
    status = "ok"
    info = {}
    try:
        if name == "G-MaxDeg":
            C = max_degree_greedy(g)
        elif name == "G-FastVC":
            C = edge_greedy_fastvc(g)
        elif name == "2-Approx":
            C = matching_2approx(g)
        elif name == "R-MC(g=1,R=16)":
            C = randomized_greedy(g, gamma=1.0, runs=16, seed=1)
        elif name == "R-MIS":
            C = random_order_mis(g, seed=1)
        elif name == "DP-TD":
            C = dp_vertex_cover(g, max_width=16, deadline=Deadline(limit))
            if C is None:
                status, C = "tw>16", None
        elif name == "DC-WH":
            C = divide_and_conquer(g, Deadline(limit))
        elif name == "HyVC":
            C, info = hybrid_vc(g, time_limit=limit, runs=16, gamma=2.0, seed=1)
            status = "ok" if info["exact"] else "anytime"
        else:
            raise ValueError(name)
    except Timeout:
        status, C = "timeout", None
    dt = time.perf_counter() - t0
    if C is not None:
        assert is_vertex_cover(g, C), name
    return (len(C) if C is not None else ""), dt, status, info


ALGS = ["G-MaxDeg", "G-FastVC", "2-Approx", "R-MC(g=1,R=16)", "R-MIS", "DP-TD", "DC-WH", "HyVC"]


def main_bench(limit, out, ref_limit):
    with open(out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["family", "instance", "param", "n", "m", "tw_ub", "alg", "size", "time_s",
                    "status", "opt", "opt_proven", "hy_dp_comp", "hy_bb_comp", "hy_nodes",
                    "hy_reduced_n"])
        for fam, name, g, param in instances():
            n, m = len(g), num_edges(g)
            tw = treewidth_upper_bound(g)
            ref, proven = cpsat_vertex_cover(g, ref_limit, workers=2)
            rows = []
            best_exact = None
            for a in ALGS:
                size, dt, st, info = run_alg(a, g, limit)
                if a in ("DP-TD", "DC-WH") and st == "ok" or (a == "HyVC" and st == "ok"):
                    best_exact = size
                rows.append((a, size, dt, st, info))
            opt = len(ref) if ref is not None else ""
            if not proven and best_exact is not None:
                opt, proven = best_exact, True  # un algoritmo exacto nuestro lo certificó
            for a, size, dt, st, info in rows:
                w.writerow([fam, name, param, n, m, tw, a, size, f"{dt:.4f}", st, opt, proven,
                            info.get("dp_components", ""), info.get("bb_components", ""),
                            info.get("nodes", ""), info.get("reduced_n", "")])
            f.flush()
            print(f"{name:22s} n={n:5d} m={m:6d} tw<={tw:3d} opt={opt}{'' if proven else '?'} | " +
                  " ".join(f"{a}={s}({st[0]},{dt:.2f}s)" for a, s, dt, st, _ in rows), flush=True)


def ablation(limit, out):
    variants = {
        "HyVC completo": {},
        "sin reducciones": {"use_reduce": False},
        "sin DP": {"use_dp": False},
        "sin division": {"use_split": False},
        "sin Monte Carlo": {"use_mc": False},
    }
    with open(out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["family", "instance", "variant", "size", "time_s", "exact", "nodes"])
        for fam, name, g, param in instances():
            if fam not in ("ER-n200", "Banded", "BA-n1000") or name.endswith("s2"):
                continue
            for vname, kw in variants.items():
                t0 = time.perf_counter()
                C, info = hybrid_vc(g, time_limit=limit, runs=16, gamma=2.0, seed=1, **kw)
                dt = time.perf_counter() - t0
                w.writerow([fam, name, vname, len(C), f"{dt:.4f}", info["exact"], info["nodes"]])
                f.flush()
                print(name, vname, len(C), f"{dt:.2f}", info["exact"], info["nodes"], flush=True)


def gamma_sweep(out):
    graphs = [("GreedyTrap", "trap_k80", greedy_trap(80), 80)] + \
             [("ER-n500", f"er500_c4_s{s}", erdos_renyi(500, 4.0, 200 + s), None) for s in range(2)] + \
             [("BA-n1000", "ba1000_m3_s0", barabasi_albert(1000, 3, 400), None)]
    with open(out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["instance", "opt", "gamma", "runs", "size", "time_s"])
        for fam, name, g, opt in graphs:
            if opt is None:
                C, info = hybrid_vc(g, time_limit=120)
                opt = len(C) if info["exact"] else ""
            for gamma in (0.0, 1.0, 2.0, 4.0, 8.0, 32.0):
                for runs in (1, 16):
                    t0 = time.perf_counter()
                    C = randomized_greedy(g, gamma=gamma, runs=runs, seed=3)
                    w.writerow([name, opt, gamma, runs, len(C), f"{time.perf_counter() - t0:.4f}"])
                    print(name, opt, gamma, runs, len(C), flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=float, default=20.0)
    ap.add_argument("--ref-limit", type=float, default=30.0)
    ap.add_argument("--out", default="results/main.csv")
    ap.add_argument("--ablation", action="store_true")
    ap.add_argument("--gamma", action="store_true")
    a = ap.parse_args()
    if a.ablation:
        ablation(a.limit, a.out)
    elif a.gamma:
        gamma_sweep(a.out)
    else:
        main_bench(a.limit, a.out, a.ref_limit)
