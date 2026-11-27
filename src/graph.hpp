// graph.hpp — grafo simple no dirigido, lectura de formatos públicos y utilidades.
#pragma once
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace vc {

struct Graph {
    int n = 0;
    std::vector<std::vector<int>> adj;  // listas ordenadas y sin repetidos

    explicit Graph(int n_ = 0) : n(n_), adj(n_) {}

    void add_edge(int u, int v) {
        if (u == v) return;
        adj[u].push_back(v);
        adj[v].push_back(u);
    }
    // Ordena y elimina aristas repetidas.
    void finalize() {
        for (auto& l : adj) {
            std::sort(l.begin(), l.end());
            l.erase(std::unique(l.begin(), l.end()), l.end());
        }
    }
    long long m() const {
        long long s = 0;
        for (auto& l : adj) s += (long long)l.size();
        return s / 2;
    }
    int max_degree() const {
        int d = 0;
        for (auto& l : adj) d = std::max(d, (int)l.size());
        return d;
    }
};

// Subgrafo inducido por `vs` (ids locales 0..|vs|-1, orig[i] = vs[i]).
inline Graph induced(const Graph& g, const std::vector<int>& vs, std::vector<int>& local_of) {
    Graph h((int)vs.size());
    if ((int)local_of.size() < g.n) local_of.assign(g.n, -1);
    for (int i = 0; i < (int)vs.size(); ++i) local_of[vs[i]] = i;
    for (int i = 0; i < (int)vs.size(); ++i)
        for (int w : g.adj[vs[i]])
            if (local_of[w] >= 0) h.adj[i].push_back(local_of[w]);
    for (int v : vs) local_of[v] = -1;
    for (auto& l : h.adj) std::sort(l.begin(), l.end());
    return h;
}

inline std::vector<std::vector<int>> components(const Graph& g) {
    std::vector<int> comp(g.n, -1);
    std::vector<std::vector<int>> out;
    std::vector<int> st;
    for (int s = 0; s < g.n; ++s) {
        if (comp[s] >= 0) continue;
        int c = (int)out.size();
        out.emplace_back();
        comp[s] = c;
        st.push_back(s);
        while (!st.empty()) {
            int x = st.back();
            st.pop_back();
            out[c].push_back(x);
            for (int y : g.adj[x])
                if (comp[y] < 0) { comp[y] = c; st.push_back(y); }
        }
    }
    return out;
}

inline bool is_vertex_cover(const Graph& g, const std::vector<int>& cover) {
    std::vector<char> in(g.n, 0);
    for (int v : cover) {
        if (v < 0 || v >= g.n) return false;
        in[v] = 1;
    }
    for (int u = 0; u < g.n; ++u)
        if (!in[u])
            for (int w : g.adj[u])
                if (!in[w]) return false;
    return true;
}

// Post-proceso común: quita vértices cuyos vecinos ya están todos en la cubierta.
inline std::vector<int> remove_redundant(const Graph& g, std::vector<int> cover) {
    std::vector<char> in(g.n, 0);
    for (int v : cover) in[v] = 1;
    std::sort(cover.begin(), cover.end(), [&](int a, int b) { return g.adj[a].size() < g.adj[b].size(); });
    std::vector<int> out;
    for (int v : cover) {
        bool red = true;
        for (int w : g.adj[v]) if (!in[w]) { red = false; break; }
        if (red) in[v] = 0;
    }
    for (int v = 0; v < g.n; ++v) if (in[v]) out.push_back(v);
    return out;
}

// Cota inferior barata: max(|matching maximal|, ceil(m / Δ)).
inline int lower_bound_vc(const Graph& g) {
    long long m = g.m();
    if (m == 0) return 0;
    std::vector<int> order(g.n);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int a, int b) { return g.adj[a].size() < g.adj[b].size(); });
    std::vector<char> mt(g.n, 0);
    int M = 0;
    for (int u : order) {
        if (mt[u]) continue;
        int best = -1;
        for (int w : g.adj[u])
            if (!mt[w] && (best < 0 || g.adj[w].size() < g.adj[best].size())) best = w;
        if (best >= 0) { mt[u] = mt[best] = 1; ++M; }
    }
    int D = g.max_degree();
    return std::max<long long>(M, (m + D - 1) / D);
}

// ------------------------------------------------------------------ lectura de archivos
// Acepta: PACE ("p td n m" / "p tw n m" + "u v"), DIMACS ("p edge n m" + "e u v"),
// y lista de aristas simple "u v" (0- o 1-indexada, se detecta).
inline Graph read_graph(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("no se pudo abrir " + path);
    std::string line;
    int n = -1;
    std::vector<std::pair<long long, long long>> E;
    long long maxid = -1, minid = (1LL << 60);
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == 'c' || line[0] == '%' || line[0] == '#') continue;
        std::istringstream ss(line);
        if (line[0] == 'p') {
            std::string p, fmt;
            long long nn, mm;
            ss >> p >> fmt >> nn >> mm;
            n = (int)nn;
            continue;
        }
        if (line[0] == 'e') { char e; ss >> e; }
        long long u, v;
        if (ss >> u >> v) {
            E.emplace_back(u, v);
            maxid = std::max({maxid, u, v});
            minid = std::min({minid, u, v});
        }
    }
    bool one_based = (n >= 0) ? true : (minid >= 1);
    if (n < 0) n = (int)(one_based ? maxid : maxid + 1);
    Graph g(n);
    for (auto [u, v] : E) g.add_edge((int)(u - one_based), (int)(v - one_based));
    g.finalize();
    return g;
}

struct Timer {
    std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();
    double seconds() const {
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    }
};

struct Timeout : std::exception {};

struct Deadline {
    Timer t;
    double limit;
    unsigned long long calls = 0;
    explicit Deadline(double lim = 1e18) : limit(lim) {}
    // Se consulta el reloj en cada llamada (cuesta ~20 ns): con nodos de búsqueda caros, revisar
    // solo cada 256 llamadas dejaba pasar varios segundos de más.
    void check() {
        ++calls;
        if (t.seconds() > limit) throw Timeout();
    }
    bool expired() const { return t.seconds() > limit; }
};

}  // namespace vc
