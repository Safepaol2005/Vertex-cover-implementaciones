// dp_td.hpp — Programación dinámica sobre una descomposición en árbol (exacta, O*(2^tw)).
//
// Esquema de Bannach & Berndt (ESA 2018; Algorithms 12(8), 2019) y Alber, Dorn &
// Niedermeier (DAM 145, 2005): una tabla por bolsa B_t indexada por el subconjunto S ⊆ B_t
// que está en la cubierta; S es válido si cubre las aristas originales dentro de B_t.
//
// Descomposición: eliminación por grado mínimo. bag(v) = N(v) ∪ {v} al eliminar v; el padre
// de bag(v) es el vecino de v eliminado primero después de v. Con esta construcción
// bag(hijo) \ {hijo} ⊆ bag(padre), de modo que cada hijo entrega una tabla ya "olvidada".
#pragma once
#include <queue>
#include <unordered_set>

#include "graph.hpp"

namespace vc {

struct TreeDec {
    std::vector<int> order;                 // orden de eliminación
    std::vector<std::vector<int>> bag;      // bag[v] = vecinos posteriores..., v (v al final)
    int width = -1;
    unsigned long long cost = 0;            // Σ 2^{|bag|}  (celdas de tabla)
};

// Devuelve false si el ancho supera max_width (aborta temprano).
inline bool min_degree_td(const Graph& g, int max_width, TreeDec& td) {
    const int n = g.n;
    std::vector<std::unordered_set<int>> h(n);
    for (int v = 0; v < n; ++v) h[v].insert(g.adj[v].begin(), g.adj[v].end());
    using P = std::pair<int, int>;
    std::priority_queue<P, std::vector<P>, std::greater<P>> pq;
    for (int v = 0; v < n; ++v) pq.push({(int)h[v].size(), v});
    std::vector<char> done(n, 0);
    td.order.clear();
    td.bag.assign(n, {});
    td.width = 0;
    td.cost = 0;
    while (!pq.empty()) {
        auto [d, v] = pq.top();
        pq.pop();
        if (done[v] || d != (int)h[v].size()) continue;
        std::vector<int> N(h[v].begin(), h[v].end());
        td.width = std::max(td.width, (int)N.size());
        if (td.width > max_width) return false;
        td.bag[v] = N;
        td.bag[v].push_back(v);
        td.cost += 1ULL << td.bag[v].size();
        td.order.push_back(v);
        done[v] = 1;
        for (int a : N) h[a].erase(v);
        for (size_t i = 0; i < N.size(); ++i)
            for (size_t j = i + 1; j < N.size(); ++j)
                if (h[N[i]].insert(N[j]).second) h[N[j]].insert(N[i]);
        for (int a : N) pq.push({(int)h[a].size(), a});
        h[v].clear();
    }
    return true;
}

// Celdas de tabla que costaría la PD, o UINT64_MAX si el ancho supera max_width.
inline unsigned long long dp_cost(const Graph& g, int max_width) {
    TreeDec td;
    if (!min_degree_td(g, max_width, td)) return ~0ULL;
    return td.cost;
}

// Cubierta mínima exacta; devuelve false si el ancho supera max_width.
// También devuelve false si la tabla total (Σ 2^|bolsa|) supera max_cost (límite de memoria).
inline bool dp_vertex_cover(const Graph& g, int max_width, std::vector<int>& cover,
                            Deadline* dl = nullptr, unsigned long long max_cost = 1ULL << 27) {
    cover.clear();
    if (g.n == 0) return true;
    TreeDec td;
    if (!min_degree_td(g, max_width, td)) return false;
    if (td.cost > max_cost) return false;
    const int n = g.n;
    std::vector<int> pos(n);
    for (int i = 0; i < n; ++i) pos[td.order[i]] = i;
    std::vector<int> parent(n, -1);
    std::vector<std::vector<int>> children(n);
    for (int v : td.order) {
        const auto& B = td.bag[v];
        int p = -1;
        for (size_t i = 0; i + 1 < B.size(); ++i)
            if (p < 0 || pos[B[i]] < pos[p]) p = B[i];
        parent[v] = p;
        if (p >= 0) children[p].push_back(v);
    }
    const uint32_t INF = 0x3fffffff;
    std::vector<std::vector<uint32_t>> forgotten(n);  // g_c[m]: tabla del hijo, vértice olvidado
    std::vector<std::vector<uint8_t>> argbit(n);      // decisión del vértice olvidado
    std::vector<std::vector<int>> proj(n);            // posiciones de bag(c)\{c} en bag(padre)
    std::vector<int> local(n, -1);
    std::vector<uint32_t> T;
    for (int v : td.order) {
        if (dl) dl->check();
        const auto& B = td.bag[v];
        const int b = (int)B.size();
        for (int i = 0; i < b; ++i) local[B[i]] = i;
        // máscara de vecinos ORIGINALES dentro de la bolsa
        std::vector<uint32_t> nb(b, 0);
        for (int i = 0; i < b; ++i)
            for (int w : g.adj[B[i]])
                if (local[w] >= 0) nb[i] |= 1u << local[w];
        const uint32_t full = 1u << b;
        T.assign(full, 0);
        for (uint32_t S = 0; S < full; ++S) {
            bool ok = true;
            for (int i = 0; i < b && ok; ++i)
                if (!(S >> i & 1) && (S & nb[i]) != nb[i]) ok = false;
            T[S] = ok ? (uint32_t)__builtin_popcount(S) : INF;
        }
        for (int c : children[v]) {
            const auto& Bc = td.bag[c];
            auto& pp = proj[c];
            pp.resize(Bc.size() - 1);
            for (size_t j = 0; j + 1 < Bc.size(); ++j) pp[j] = local[Bc[j]];
            const auto& gc = forgotten[c];
            for (uint32_t S = 0; S < full; ++S) {
                if (T[S] >= INF) continue;
                uint32_t idx = 0;
                for (size_t j = 0; j < pp.size(); ++j) idx |= ((S >> pp[j]) & 1u) << j;
                uint32_t val = gc[idx];
                if (val >= INF) { T[S] = INF; continue; }
                T[S] += val - (uint32_t)__builtin_popcount(idx);
            }
            std::vector<uint32_t>().swap(forgotten[c]);
        }
        for (int i = 0; i < b; ++i) local[B[i]] = -1;
        // v es el bit más alto: olvidarlo = mínimo de las dos mitades
        const uint32_t half = full >> 1;
        forgotten[v].resize(half);
        argbit[v].resize(half);
        for (uint32_t S = 0; S < half; ++S) {
            uint32_t lo = T[S], hi = T[S | half];
            argbit[v][S] = hi < lo;
            forgotten[v][S] = std::min(lo, hi);
        }
    }
    // reconstrucción descendente
    std::vector<uint32_t> chosen(n, 0);
    for (int i = n - 1; i >= 0; --i) {
        int v = td.order[i];
        const int b = (int)td.bag[v].size();
        uint32_t m = 0;
        if (parent[v] >= 0) {
            uint32_t pm = chosen[parent[v]];
            for (size_t j = 0; j < proj[v].size(); ++j) m |= ((pm >> proj[v][j]) & 1u) << j;
        }
        uint32_t top = argbit[v][m];
        chosen[v] = m | (top << (b - 1));
        if (top) cover.push_back(v);
    }
    return true;
}

}  // namespace vc
