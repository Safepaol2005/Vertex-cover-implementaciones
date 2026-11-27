// hybrid.hpp — HyVC: algoritmo híbrido construido a partir de los tres algoritmos estudiados.
//
//   1. REDUCIR    reglas exactas R0-R5 (reductions.hpp), incluida la LP de Nemhauser–Trotter.
//   2. DIVIDIR    componentes conexas, resueltas por separado (idea de Weigt & Hartmann).
//   3. DESPACHAR  si la tabla de la PD es barata (Σ 2^|bolsa| <= tope) -> PD exacta.
//   4. ACOTAR     cota superior: voraz de grado máximo y R corridas del MIS aleatorio de
//                 Veldt, sin redundantes. Cota inferior: max(LP de Nemhauser–Trotter,
//                 cobertura por cliques, matching maximal).
//   5. RAMIFICAR  sobre el vértice de grado máximo (v ∈ C  |  N(v) ⊆ C), con presupuesto;
//                 en cada nodo se vuelve a reducir, dividir y despachar.
//   6. ANYTIME    si se agota el tiempo se devuelve la mejor cubierta conocida.
// La salida siempre es una cubierta válida; si `exact` es verdadero, es mínima.
#pragma once
#include <optional>

#include "bounds.hpp"
#include "dp_td.hpp"
#include "greedy_wh.hpp"
#include "random_mis.hpp"
#include "reductions.hpp"

namespace vc {

struct HybridParams {
    double time_limit = 30.0;
    unsigned long long cost_cap = 1ULL << 26;   // tope de celdas de la PD en la raíz
    unsigned long long inner_cap = 1ULL << 16;  // tope en nodos internos con componentes pequeñas
    int big_comp = 100;                         // componentes > big_comp usan cost_cap
    int max_width = 22;
    int mc_runs = 64;
    unsigned long long seed = 1;
    bool use_reduce = true, use_split = true, use_dp = true, use_mc = true;
    bool use_strong = true;  // reglas R4' (no confinado) y R5 (LP de Nemhauser–Trotter)
};

struct HybridStats {
    bool exact = true;
    long long nodes = 0, dp_calls = 0;
    int kernel_n = 0, components = 0, lower_bound = 0;
};

class Hybrid {
public:
    Hybrid(const HybridParams& p) : P(p), dl(p.time_limit) {}

    std::vector<int> run(const Graph& g, HybridStats& st) {
        stats = &st;
        Kernel K = P.use_reduce ? reduce(g, true, true, P.use_strong, P.use_strong) : trivial_kernel(g);
        st.kernel_n = K.k.n;
        st.lower_bound = K.offset();
        std::vector<int> kc;
        auto comps = P.use_split ? components(K.k) : std::vector<std::vector<int>>{all(K.k.n)};
        std::sort(comps.begin(), comps.end(), [](auto& a, auto& b) { return a.size() < b.size(); });
        std::vector<int> local;
        for (auto& c : comps) {
            if (c.size() < 2) continue;
            ++st.components;
            Graph h = induced(K.k, c, local);
            int lb = best_lower_bound(h);
            st.lower_bound += lb;
            std::vector<int> sol;
            if (P.use_dp && dp_ok(h, P.cost_cap)) {
                ++st.dp_calls;
                dp_vertex_cover(h, P.max_width, sol);
            } else {
                std::vector<int> ub = upper_bound(h, true);
                if (dl.expired()) { st.exact = false; sol = ub; }
                else {
                    try {
                        auto better = solve_component(h, (int)ub.size(), lb);
                        sol = better ? *better : ub;  // sin solución < |ub| => ub es óptima
                    } catch (Timeout&) { st.exact = false; sol = ub; }
                }
            }
            for (int x : sol) kc.push_back(c[x]);
        }
        auto cover = unfold(K, kc);
        if (!st.exact) cover = remove_redundant(g, cover);
        return cover;
    }

private:
    HybridParams P;
    Deadline dl;
    HybridStats* stats = nullptr;

    static std::vector<int> all(int n) { std::vector<int> v(n); std::iota(v.begin(), v.end(), 0); return v; }

    static Kernel trivial_kernel(const Graph& g) {
        Kernel K;
        K.n_input = g.n;
        for (int v = 0; v < g.n; ++v) if (!g.adj[v].empty()) K.work.push_back(v);
        std::vector<int> local;
        K.k = induced(g, K.work, local);
        return K;
    }

    bool dp_ok(const Graph& h, unsigned long long cap) const {
        if (h.m() > (long long)h.n * P.max_width) return false;  // degeneración > ancho máx.
        return dp_cost(h, P.max_width) <= cap;
    }

    std::vector<int> upper_bound(const Graph& h, bool root) {
        auto best = remove_redundant(h, greedy_wh(h));
        if (P.use_mc) {
            auto r = random_mis(h, root ? P.mc_runs : 4, P.seed + (unsigned long long)stats->nodes, true);
            if (r.size() < best.size()) best = std::move(r);
        }
        return best;
    }

    // Cubierta mínima de h si su tamaño es < ub; nullopt si no existe ninguna < ub.
    std::optional<std::vector<int>> solve(const Graph& g, int ub) {
        dl.check();
        ++stats->nodes;
        Kernel K = P.use_reduce ? reduce(g, true, true, P.use_strong, P.use_strong) : trivial_kernel(g);
        int off = K.offset();
        if (off >= ub) return std::nullopt;
        if (K.k.m() == 0) return unfold(K, {});
        int budget = ub - off;
        auto comps = P.use_split ? components(K.k) : std::vector<std::vector<int>>{all(K.k.n)};
        std::sort(comps.begin(), comps.end(), [](auto& a, auto& b) { return a.size() < b.size(); });
        std::vector<Graph> subs;
        std::vector<std::vector<int>> keep;
        std::vector<int> lbs, local;
        long long sum_lb = 0;
        for (auto& c : comps) {
            if (c.size() < 2) continue;
            subs.push_back(induced(K.k, c, local));
            keep.push_back(c);
            lbs.push_back(best_lower_bound(subs.back()));
            sum_lb += lbs.back();
        }
        if (sum_lb >= budget) return std::nullopt;
        std::vector<int> kc;
        long long used = 0, rest = sum_lb;
        for (size_t i = 0; i < subs.size(); ++i) {
            rest -= lbs[i];
            int cb = (int)(budget - used - rest);
            auto s = solve_component(subs[i], cb, lbs[i]);
            if (!s) return std::nullopt;
            for (int x : *s) kc.push_back(keep[i][x]);
            used += (long long)s->size();
        }
        return unfold(K, kc);
    }

    std::optional<std::vector<int>> solve_component(const Graph& h, int ub, int lb) {
        if (lb >= ub) return std::nullopt;
        unsigned long long cap = h.n > P.big_comp ? P.cost_cap : std::min(P.cost_cap, P.inner_cap);
        if (P.use_dp && dp_ok(h, cap)) {
            ++stats->dp_calls;
            std::vector<int> s;
            dp_vertex_cover(h, P.max_width, s, &dl);
            if ((int)s.size() < ub) return s;
            return std::nullopt;
        }
        std::optional<std::vector<int>> best;
        auto u = upper_bound(h, false);
        if ((int)u.size() < ub) { ub = (int)u.size(); best = u; }
        if (lb >= ub) return best;
        int v = 0;
        for (int x = 1; x < h.n; ++x) if (h.adj[x].size() > h.adj[v].size()) v = x;
        std::vector<int> local;
        // rama 1: v en la cubierta
        {
            std::vector<int> rest;
            for (int x = 0; x < h.n; ++x) if (x != v) rest.push_back(x);
            Graph g1 = induced(h, rest, local);
            auto s = solve(g1, ub - 1);
            if (s) {
                std::vector<int> c;
                for (int x : *s) c.push_back(rest[x]);
                c.push_back(v);
                best = c;
                ub = (int)c.size();
            }
        }
        // rama 2: v fuera => N(v) en la cubierta
        int dv = (int)h.adj[v].size();
        if (dv < ub) {
            std::vector<char> drop(h.n, 0);
            drop[v] = 1;
            for (int w : h.adj[v]) drop[w] = 1;
            std::vector<int> rest;
            for (int x = 0; x < h.n; ++x) if (!drop[x]) rest.push_back(x);
            Graph g2 = induced(h, rest, local);
            auto s = solve(g2, ub - dv);
            if (s) {
                std::vector<int> c(h.adj[v].begin(), h.adj[v].end());
                for (int x : *s) c.push_back(rest[x]);
                best = c;
            }
        }
        return best;
    }
};

inline std::vector<int> hybrid_vc(const Graph& g, const HybridParams& p, HybridStats& st) {
    Hybrid H(p);
    return H.run(g, st);
}

}  // namespace vc
