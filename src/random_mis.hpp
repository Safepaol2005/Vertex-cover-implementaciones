// random_mis.hpp — Algoritmo aleatorizado de Veldt (SOSA 2024):
// "Growing a random maximal independent set produces a 2-approximate vertex cover".
//
// Se recorre una permutación uniforme de V; cada vértice todavía libre entra al
// conjunto independiente I y sus vecinos libres pasan a la cubierta C = V \ I.
// Tiempo O(n + m) por corrida; E|C| <= 2 OPT.  Repetir R veces y quedarse con la mejor
// cubierta es un esquema Monte Carlo (la cubierta siempre es válida).
#pragma once
#include "graph.hpp"

namespace vc {

inline std::vector<int> random_mis_once(const Graph& g, std::mt19937_64& rng, std::vector<int>& perm) {
    const int n = g.n;
    if ((int)perm.size() != n) { perm.resize(n); std::iota(perm.begin(), perm.end(), 0); }
    std::shuffle(perm.begin(), perm.end(), rng);
    std::vector<char> state(n, 0);  // 0 libre, 1 en I, 2 en C
    std::vector<int> cover;
    for (int v : perm) {
        if (state[v]) continue;
        state[v] = 1;
        for (int w : g.adj[v])
            if (!state[w]) { state[w] = 2; cover.push_back(w); }
    }
    return cover;
}

inline std::vector<int> random_mis(const Graph& g, int runs = 1, unsigned long long seed = 1,
                                   bool prune = false) {
    std::mt19937_64 rng(seed);
    std::vector<int> perm, best;
    bool have = false;
    for (int r = 0; r < runs; ++r) {
        auto c = random_mis_once(g, rng, perm);
        if (prune) c = remove_redundant(g, std::move(c));
        if (!have || c.size() < best.size()) { best = std::move(c); have = true; }
    }
    return best;
}

}  // namespace vc
