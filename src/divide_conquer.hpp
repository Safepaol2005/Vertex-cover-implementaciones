// divide_conquer.hpp — Divide y vencerás exacto de Weigt & Hartmann (PRE 63, 056127, 2001).
//
// Toma el vértice libre de mayor grado i y abre dos ramas: i ∈ C (borrar i) o i ∉ C
// (entonces N(i) ⊆ C; borrar N[i]). Tras cada rama el grafo residual se parte en
// componentes conexas que se resuelven por separado y se suman. Sin cotas: es el punto de
// comparación exponencial de este trabajo.
#pragma once
#include "graph.hpp"

namespace vc {

namespace detail {
inline std::vector<int> dc_solve(const Graph& g, Deadline& dl) {
    dl.check();
    if (g.m() == 0) return {};
    auto comps = components(g);
    std::vector<int> local;
    if (comps.size() > 1) {
        std::vector<int> out;
        for (auto& c : comps) {
            if (c.size() < 2) continue;
            Graph h = induced(g, c, local);
            for (int x : dc_solve(h, dl)) out.push_back(c[x]);
        }
        return out;
    }
    int i = 0;
    for (int v = 1; v < g.n; ++v)
        if (g.adj[v].size() > g.adj[i].size()) i = v;
    // rama 1: i en la cubierta
    std::vector<int> rest1;
    for (int v = 0; v < g.n; ++v) if (v != i) rest1.push_back(v);
    Graph g1 = induced(g, rest1, local);
    std::vector<int> c1;
    for (int x : dc_solve(g1, dl)) c1.push_back(rest1[x]);
    c1.push_back(i);
    // rama 2: i fuera => N(i) dentro
    std::vector<char> drop(g.n, 0);
    drop[i] = 1;
    for (int w : g.adj[i]) drop[w] = 1;
    std::vector<int> rest2;
    for (int v = 0; v < g.n; ++v) if (!drop[v]) rest2.push_back(v);
    Graph g2 = induced(g, rest2, local);
    std::vector<int> c2(g.adj[i].begin(), g.adj[i].end());
    for (int x : dc_solve(g2, dl)) c2.push_back(rest2[x]);
    return c1.size() <= c2.size() ? c1 : c2;
}
}  // namespace detail

inline std::vector<int> divide_and_conquer(const Graph& g, Deadline& dl) {
    return detail::dc_solve(g, dl);
}

}  // namespace vc
