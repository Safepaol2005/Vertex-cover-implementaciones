// greedy_wh.hpp — Heurística voraz de grado máximo (Weigt & Hartmann, PRE 63, 056127, 2001).
//
// En cada paso se elige el vértice de mayor grado residual, se pone en la cubierta y se
// borran sus aristas. Implementación con cubetas por grado: O(n + m).
#pragma once
#include "graph.hpp"

namespace vc {

inline std::vector<int> greedy_wh(const Graph& g) {
    const int n = g.n;
    int D = g.max_degree();
    std::vector<int> deg(n), pos(n);
    std::vector<std::vector<int>> bucket(D + 1);
    for (int v = 0; v < n; ++v) {
        deg[v] = (int)g.adj[v].size();
        pos[v] = (int)bucket[deg[v]].size();
        bucket[deg[v]].push_back(v);
    }
    std::vector<char> removed(n, 0);
    std::vector<int> cover;
    auto bucket_remove = [&](int v) {  // O(1): intercambio con el último
        auto& b = bucket[deg[v]];
        int last = b.back();
        b[pos[v]] = last;
        pos[last] = pos[v];
        b.pop_back();
    };
    int cur = D;
    while (cur > 0) {
        if (bucket[cur].empty()) { --cur; continue; }
        // desempate determinista: el último de la cubeta
        int v = bucket[cur].back();
        bucket_remove(v);
        removed[v] = 1;
        cover.push_back(v);
        for (int w : g.adj[v]) {
            if (removed[w]) continue;
            bucket_remove(w);
            --deg[w];
            pos[w] = (int)bucket[deg[w]].size();
            bucket[deg[w]].push_back(w);
        }
        deg[v] = 0;
    }
    return cover;
}

}  // namespace vc
