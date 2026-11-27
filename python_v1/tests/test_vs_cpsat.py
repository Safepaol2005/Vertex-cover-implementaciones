"""Verificación en grafos medianos contra CP-SAT (requiere ortools)."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from vc import hybrid_vc, dp_vertex_cover, divide_and_conquer, is_vertex_cover, Deadline, Timeout
from vc.exact_ref import cpsat_vertex_cover
from vc.graph import erdos_renyi, barabasi_albert, partial_k_tree

def main():
    cases = [(f"er{n}_{c}_{s}", erdos_renyi(n, c, s)) for n in (40, 70) for c in (2, 3, 5) for s in range(4)]
    cases += [(f"ba{n}_{s}", barabasi_albert(n, 3, s)) for n in (80, 150) for s in range(3)]
    cases += [(f"pkt{s}", partial_k_tree(120, 6, 0.7, s)) for s in range(3)]
    bad = 0
    for name, g in cases:
        opt, ok = cpsat_vertex_cover(g, 30)
        if not ok:
            C, info = hybrid_vc(g, time_limit=60)
            print(name, 'CP-SAT sin probar optimo:', len(opt), 'HyVC:', len(C), info['exact'])
            assert len(C) <= len(opt)
            continue
        print(name, len(g), len(opt), flush=True)
        C, info = hybrid_vc(g, time_limit=60, cost_cap=1 << 12)
        assert is_vertex_cover(g, C)
        if info["exact"] and len(C) != len(opt):
            bad += 1; print("MISMATCH hybrid", name, len(C), len(opt))
        C2, _ = hybrid_vc(g, time_limit=60, use_dp=False)
        if len(C2) != len(opt):
            bad += 1; print("MISMATCH hybrid-noDP", name, len(C2), len(opt))
        D = dp_vertex_cover(g, max_width=14)
        if D is not None and len(D) != len(opt):
            bad += 1; print("MISMATCH dp", name, len(D), len(opt))
        try:
            E = divide_and_conquer(g, Deadline(5))
            if len(E) != len(opt):
                bad += 1; print("MISMATCH dc", name, len(E), len(opt))
        except Timeout:
            pass
    print("cases", len(cases), "mismatches", bad)
    assert bad == 0

if __name__ == "__main__":
    main()
