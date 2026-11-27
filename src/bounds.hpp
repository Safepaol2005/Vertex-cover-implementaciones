// bounds.hpp — Cotas inferiores para Vertex Cover usadas por HyVC.
//
//  * LP (Nemhauser & Trotter 1975): el óptimo de la relajación lineal
//      min Σ x_v  s.a.  x_u + x_v >= 1,  0 <= x <= 1
//    es semientero y vale |M*|/2, donde M* es un matching máximo del grafo bipartito
//    "doble" B(G) = (V_L ∪ V_R, {u_L v_R : uv ∈ E}). Se calcula con Hopcroft–Karp
//    en O(m √n).  Cota: ceil(|M*| / 2).
//  * Cobertura por cliques: si V se parte en cliques Q_1..Q_k, toda cubierta toma al menos
//    |Q_i| - 1 vértices de cada Q_i; cota n - k. Partición voraz en O(m).
//  * Matching maximal (2-aprox. clásica, Gavril): cota |M|.
#pragma once
#include <functional>
#include <queue>

#include "graph.hpp"

namespace vc {

// Hopcroft–Karp sobre el doble bipartito: ambos lados son copias de V.
inline int lp_bound(const Graph& g) {
    const int n = g.n;
    if (g.m() == 0) return 0;
    const int NIL = -1;
    std::vector<int> mL(n, NIL), mR(n, NIL), dist(n);
    std::vector<int> q(n);
    auto bfs = [&]() {
        int h = 0, t = 0;
        bool found = false;
        for (int u = 0; u < n; ++u) {
            if (mL[u] == NIL) { dist[u] = 0; q[t++] = u; }
            else dist[u] = -1;
        }
        while (h < t) {
            int u = q[h++];
            for (int v : g.adj[u]) {
                int w = mR[v];
                if (w == NIL) found = true;
                else if (dist[w] < 0) { dist[w] = dist[u] + 1; q[t++] = w; }
            }
        }
        return found;
    };
    std::vector<size_t> it(n);
    std::function<bool(int)> dfs = [&](int u) -> bool {
        for (size_t& i = it[u]; i < g.adj[u].size(); ++i) {
            int v = g.adj[u][i];
            int w = mR[v];
            if (w == NIL || (dist[w] == dist[u] + 1 && dfs(w))) {
                mL[u] = v; mR[v] = u; ++i;
                return true;
            }
        }
        dist[u] = -1;
        return false;
    };
    int M = 0;
    while (bfs()) {
        std::fill(it.begin(), it.end(), 0);
        for (int u = 0; u < n; ++u)
            if (mL[u] == NIL && dfs(u)) ++M;
    }
    return (M + 1) / 2;
}

inline int clique_cover_bound(const Graph& g) {
    const int n = g.n;
    std::vector<int> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int a, int b) { return g.adj[a].size() < g.adj[b].size(); });
    std::vector<int> cl(n, -1), mark(n, -1);
    std::vector<std::vector<int>> members;
    for (int v : order) {
        // marcar cuántos miembros de cada clique son vecinos de v
        std::vector<int> cand;
        for (int w : g.adj[v]) {
            int c = cl[w];
            if (c < 0) continue;
            if (mark[c] != v) { mark[c] = v; cand.push_back(c); }
        }
        int best = -1;
        for (int c : cand) {
            int cnt = 0;
            for (int w : g.adj[v]) if (cl[w] == c) ++cnt;
            if (cnt == (int)members[c].size() && (best < 0 || members[c].size() > members[best].size())) best = c;
        }
        if (best < 0) { best = (int)members.size(); members.emplace_back(); }
        members[best].push_back(v);
        cl[v] = best;
    }
    return n - (int)members.size();
}

inline int maximal_matching_bound(const Graph& g) { return lower_bound_vc(g); }

inline int best_lower_bound(const Graph& g) {
    if (g.m() == 0) return 0;
    return std::max({maximal_matching_bound(g), lp_bound(g), clique_cover_bound(g)});
}

}  // namespace vc
