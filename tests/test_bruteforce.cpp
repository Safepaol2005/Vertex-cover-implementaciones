// test_bruteforce.cpp — compara los algoritmos con fuerza bruta en grafos pequeños aleatorios.
#include <cstdio>

#include "../src/divide_conquer.hpp"
#include "../src/dp_td.hpp"
#include "../src/greedy_wh.hpp"
#include "../src/hybrid.hpp"
#include "../src/random_mis.hpp"
#include "../src/reductions.hpp"

using namespace vc;

static int brute(const Graph& g) {
    int n = g.n, best = n;
    for (uint32_t S = 0; S < (1u << n); ++S) {
        int k = __builtin_popcount(S);
        if (k >= best) continue;
        bool ok = true;
        for (int u = 0; u < n && ok; ++u)
            if (!(S >> u & 1))
                for (int w : g.adj[u]) if (!(S >> w & 1)) { ok = false; break; }
        if (ok) best = k;
    }
    return best;
}

int main() {
    std::mt19937_64 rng(7);
    int cases = 0, fails = 0;
    for (int it = 0; it < 2000; ++it) {
        int n = 1 + rng() % 16;
        double p = std::uniform_real_distribution<double>(0.05, 0.8)(rng);
        Graph g(n);
        for (int u = 0; u < n; ++u)
            for (int v = u + 1; v < n; ++v)
                if (std::uniform_real_distribution<double>(0, 1)(rng) < p) g.add_edge(u, v);
        g.finalize();
        int opt = brute(g);
        ++cases;
        auto chk = [&](const char* name, const std::vector<int>& C, bool exact) {
            if (!is_vertex_cover(g, C) || (exact && (int)C.size() != opt) || (int)C.size() < opt) {
                std::printf("FALLO %s it=%d n=%d opt=%d got=%zu\n", name, it, n, opt, C.size());
                ++fails;
            }
        };
        chk("greedy", greedy_wh(g), false);
        chk("rmis", random_mis(g, 8, it), false);
        auto r = random_mis(g, 1, it);
        if ((int)r.size() > n) ++fails;
        std::vector<int> D;
        dp_vertex_cover(g, 20, D);
        chk("dp", D, true);
        Deadline dl;
        chk("dc", divide_and_conquer(g, dl), true);
        for (int variant = 0; variant < 4; ++variant) {
            HybridParams hp;
            hp.time_limit = 1e9;
            hp.use_dp = variant != 1;
            hp.use_reduce = variant != 2;
            hp.use_split = variant != 3;
            if (variant == 1) hp.inner_cap = 0;
            HybridStats st;
            auto H = hybrid_vc(g, hp, st);
            chk("hybrid", H, true);
            if (!st.exact) ++fails;
        }
        // reducciones: OPT(G) = OPT(núcleo) + |forzados| + |plegados|
        Kernel K = reduce(g);
        std::vector<int> kc;
        dp_vertex_cover(K.k, 24, kc);
        if (brute(K.k) + K.offset() != opt) { std::printf("FALLO reducciones it=%d\n", it); ++fails; }
        chk("unfold", unfold(K, kc), true);
    }
    std::printf("%d casos, %d fallos\n", cases, fails);
    return fails ? 1 : 0;
}
