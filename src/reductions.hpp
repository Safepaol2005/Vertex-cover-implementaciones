// reductions.hpp — Reglas de reducción exactas (preservan el óptimo).
//   R0 grado 0     : borrar v.
//   R1 grado 1     : el único vecino entra en la cubierta.
//   R2 triángulo   : grado 2 con vecinos adyacentes u, w -> u, w entran.
//   R3 plegado     : grado 2 con vecinos u, w no adyacentes -> contraer {u, v, w} en x;
//                    OPT(G) = OPT(G') + 1. Al deshacer: x ∈ C' ⇒ {u, w}; si no ⇒ {v}.
//   R4 dominancia  : u ∈ N(v) y N[u] ⊆ N[v] -> v entra.
//   R4' no confinado (Xiao & Nagamochi; usada por Akiba & Iwata 2016): generaliza la
//                    dominancia. Si v es "no confinado", existe una cubierta mínima que lo contiene.
//   R5 LP (Nemhauser & Trotter 1975): con un óptimo semientero x de la relajación LP,
//                    existe una cubierta mínima que contiene {x=1} y evita {x=0}.
//                    x se obtiene de un matching máximo (Hopcroft–Karp) en el doble bipartito
//                    y de la cubierta de Kőnig correspondiente. O(m √n) por aplicación.
// Las reglas son clásicas (Cygan et al. 2015, cap. 2-3; Akiba & Iwata 2016).
#pragma once
#include <functional>

#include "graph.hpp"

namespace vc {

struct Fold { int x, v, u, w; };

struct Kernel {
    Graph k;                  // núcleo con ids locales
    std::vector<int> work;    // id de trabajo de cada vértice del núcleo
    std::vector<int> forced;  // ids de trabajo forzados en la cubierta
    std::vector<Fold> folds;  // en ids de trabajo
    int n_input = 0;          // ids de trabajo < n_input son vértices de la entrada
    int offset() const { return (int)forced.size() + (int)folds.size(); }
};

// Óptimo semientero de la relajación LP: devuelve 2·x_v ∈ {0,1,2} para cada vértice de g.
inline std::vector<int> lp_half_integral(const Graph& g) {
    const int n = g.n, NIL = -1;
    std::vector<int> mL(n, NIL), mR(n, NIL), dist(n), q(n);
    auto bfs = [&]() {
        int h = 0, t = 0; bool found = false;
        for (int u = 0; u < n; ++u) { if (mL[u] == NIL) { dist[u] = 0; q[t++] = u; } else dist[u] = -1; }
        while (h < t) {
            int u = q[h++];
            for (int v : g.adj[u]) { int w = mR[v]; if (w == NIL) found = true; else if (dist[w] < 0) { dist[w] = dist[u] + 1; q[t++] = w; } }
        }
        return found;
    };
    std::vector<size_t> it(n);
    std::function<bool(int)> dfs = [&](int u) -> bool {
        for (size_t& i = it[u]; i < g.adj[u].size(); ++i) {
            int v = g.adj[u][i], w = mR[v];
            if (w == NIL || (dist[w] == dist[u] + 1 && dfs(w))) { mL[u] = v; mR[v] = u; ++i; return true; }
        }
        dist[u] = -1; return false;
    };
    while (bfs()) { std::fill(it.begin(), it.end(), 0); for (int u = 0; u < n; ++u) if (mL[u] == NIL) dfs(u); }
    // Kőnig: Z = vértices alcanzables desde los L libres por caminos alternantes.
    std::vector<char> zL(n, 0), zR(n, 0);
    std::vector<int> st;
    for (int u = 0; u < n; ++u) if (mL[u] == NIL) { zL[u] = 1; st.push_back(u); }
    while (!st.empty()) {
        int u = st.back(); st.pop_back();
        for (int v : g.adj[u]) if (!zR[v] && mL[u] != v) { zR[v] = 1; int w = mR[v]; if (w != NIL && !zL[w]) { zL[w] = 1; st.push_back(w); } }
    }
    // cubierta mínima del bipartito: (L \ Z) ∪ (R ∩ Z);  2x_v = [v_L en ella] + [v_R en ella]
    std::vector<int> x2(n);
    for (int v = 0; v < n; ++v) x2[v] = (!zL[v]) + (zR[v] ? 1 : 0);
    return x2;
}

inline Kernel reduce(const Graph& g, bool use_fold = true, bool use_dominance = true, bool use_lp = true,
                     bool use_unconfined = true) {
    Kernel K;
    K.n_input = g.n;
    std::vector<std::vector<int>> adj = g.adj;
    std::vector<char> alive(g.n, 1);
    std::vector<int> deg(g.n), mark(g.n, 0);
    int stamp = 0;
    for (int v = 0; v < g.n; ++v) deg[v] = (int)adj[v].size();
    std::vector<int> queue;
    std::vector<char> inq(g.n, 1);
    for (int v = g.n - 1; v >= 0; --v) queue.push_back(v);
    auto push = [&](int a) { if (alive[a] && !inq[a]) { inq[a] = 1; queue.push_back(a); } };
    auto kill_take = [&](int x) {
        alive[x] = 0;
        K.forced.push_back(x);
        for (int a : adj[x]) if (alive[a]) { --deg[a]; push(a); }
    };
    auto new_vertex = [&]() {
        adj.emplace_back(); alive.push_back(1); deg.push_back(0); mark.push_back(0); inq.push_back(0);
        return (int)adj.size() - 1;
    };
    auto alive_nbrs = [&](int v, std::vector<int>& out) {
        out.clear();
        for (int a : adj[v]) if (alive[a]) out.push_back(a);
    };
    std::vector<int> N1, N2, NU;
    bool changed = true;
    while (changed) {
        changed = false;
        while (!queue.empty()) {
            int v = queue.back();
            queue.pop_back();
            inq[v] = 0;
            if (!alive[v]) continue;
            if (deg[v] == 0) { alive[v] = 0; continue; }
            if (deg[v] == 1) {
                alive_nbrs(v, N1);
                kill_take(N1[0]);
                changed = true;
                continue;
            }
            if (deg[v] == 2) {
                alive_nbrs(v, N1);
                int u = N1[0], w = N1[1];
                bool adjuw = false;
                for (int a : adj[u]) if (a == w && alive[a]) { adjuw = true; break; }
                if (adjuw) {
                    kill_take(u);
                    if (alive[w]) kill_take(w);
                    changed = true;
                } else if (use_fold) {
                    ++stamp;
                    alive_nbrs(u, N1);
                    alive_nbrs(w, N2);
                    NU.clear();
                    for (auto* L : {&N1, &N2})
                        for (int a : *L)
                            if (a != v && mark[a] != stamp) { mark[a] = stamp; NU.push_back(a); }
                    alive[v] = alive[u] = alive[w] = 0;
                    int x = new_vertex();
                    if ((int)mark.size() <= x) mark.resize(x + 1, 0);
                    adj[x] = NU;
                    deg[x] = (int)NU.size();
                    for (int a : NU) {
                        adj[a].push_back(x);
                        int d = 0;
                        for (int b : adj[a]) d += alive[b];
                        deg[a] = d;
                        push(a);
                    }
                    push(x);
                    K.folds.push_back({x, v, u, w});
                    changed = true;
                }
            }
        }
        if (use_dominance) {
            std::vector<int> order;
            for (int v = 0; v < (int)adj.size(); ++v) if (alive[v] && deg[v] > 0) order.push_back(v);
            std::sort(order.begin(), order.end(), [&](int a, int b) { return deg[a] > deg[b]; });
            for (int v : order) {
                if (!alive[v] || deg[v] == 0) continue;
                ++stamp;
                if ((int)mark.size() < (int)adj.size()) mark.resize(adj.size(), 0);
                for (int a : adj[v]) if (alive[a]) mark[a] = stamp;
                mark[v] = stamp;
                bool dom = false;
                for (int u : adj[v]) {
                    if (!alive[u] || deg[u] > deg[v]) continue;
                    bool sub = true;
                    for (int b : adj[u]) if (alive[b] && mark[b] != stamp) { sub = false; break; }
                    if (sub) { dom = true; break; }
                }
                if (dom) { kill_take(v); changed = true; }
            }
        }
        if (!changed && use_unconfined) {
            // S = {v}; mientras exista u ∈ N(S) con |N(u) ∩ S| = 1: si N(u) \ N[S] = ∅, v es
            // no confinado; si N(u) \ N[S] = {w} (para el u que minimiza), S ← S ∪ {w}; si no, parar.
            std::vector<int> S, inS_mark(adj.size(), 0), inNS_mark(adj.size(), 0);
            int st2 = 0;
            for (int v = 0; v < (int)adj.size(); ++v) {
                if (!alive[v] || deg[v] == 0) continue;
                ++st2;
                S.assign(1, v);
                inS_mark[v] = st2;
                inNS_mark[v] = st2;
                for (int a2 : adj[v]) if (alive[a2]) inNS_mark[a2] = st2;
                bool unconf = false;
                for (int iter = 0; iter < 64; ++iter) {
                    int best_cnt = 1 << 30, best_w = -1;
                    for (int s0 : S) {
                        for (int u : adj[s0]) {
                            if (!alive[u] || inS_mark[u] == st2) continue;
                            int inS = 0, out = 0, w = -1;
                            for (int b2 : adj[u]) {
                                if (!alive[b2]) continue;
                                if (inS_mark[b2] == st2) ++inS;
                                else if (inNS_mark[b2] != st2) { ++out; w = b2; }
                            }
                            if (inS == 1 && out < best_cnt) { best_cnt = out; best_w = w; }
                            if (best_cnt == 0) break;
                        }
                        if (best_cnt == 0) break;
                    }
                    if (best_cnt == 0) { unconf = true; break; }
                    if (best_cnt != 1) break;
                    S.push_back(best_w);
                    inS_mark[best_w] = st2;
                    inNS_mark[best_w] = st2;
                    for (int a2 : adj[best_w]) if (alive[a2]) inNS_mark[a2] = st2;
                }
                if (unconf) { kill_take(v); changed = true; }
            }
        }
        if (!changed && use_lp) {
            // R5: construir el grafo vivo, resolver la relajación LP y aplicar Nemhauser–Trotter
            std::vector<int> ids, loc(adj.size(), -1);
            for (int v = 0; v < (int)adj.size(); ++v) if (alive[v] && deg[v] > 0) { loc[v] = (int)ids.size(); ids.push_back(v); }
            if (!ids.empty()) {
                Graph h((int)ids.size());
                for (int i = 0; i < (int)ids.size(); ++i)
                    for (int a2 : adj[ids[i]]) if (alive[a2] && loc[a2] >= 0) h.adj[i].push_back(loc[a2]);
                h.finalize();
                auto x2 = lp_half_integral(h);
                for (int i = 0; i < h.n; ++i) if (x2[i] == 2 && alive[ids[i]]) { kill_take(ids[i]); changed = true; }
                for (int i = 0; i < h.n; ++i)
                    if (x2[i] == 0 && alive[ids[i]]) {
                        int v = ids[i];  // todos sus vecinos tienen x = 1 y ya entraron
                        alive[v] = 0;
                        for (int a2 : adj[v]) if (alive[a2]) { --deg[a2]; push(a2); }
                        changed = true;
                    }
            }
        }
    }
    // construir el núcleo
    std::vector<int> loc(adj.size(), -1);
    for (int v = 0; v < (int)adj.size(); ++v)
        if (alive[v] && deg[v] > 0) { loc[v] = (int)K.work.size(); K.work.push_back(v); }
    K.k = Graph((int)K.work.size());
    for (int i = 0; i < (int)K.work.size(); ++i)
        for (int a : adj[K.work[i]])
            if (alive[a] && loc[a] >= 0) K.k.adj[i].push_back(loc[a]);
    K.k.finalize();
    return K;
}

// Deshace las reglas: recibe la cubierta del núcleo (ids locales) y devuelve ids de la entrada.
inline std::vector<int> unfold(const Kernel& K, const std::vector<int>& kcover) {
    int total = K.n_input + (int)K.folds.size() + 1;
    for (auto& f : K.folds) total = std::max(total, f.x + 1);
    std::vector<char> in(total, 0);
    for (int x : kcover) in[K.work[x]] = 1;
    for (int x : K.forced) in[x] = 1;
    for (int i = (int)K.folds.size() - 1; i >= 0; --i) {
        const Fold& f = K.folds[i];
        if (in[f.x]) { in[f.x] = 0; in[f.u] = in[f.w] = 1; }
        else in[f.v] = 1;
    }
    std::vector<int> out;
    for (int v = 0; v < K.n_input; ++v) if (in[v]) out.push_back(v);
    return out;
}

}  // namespace vc
