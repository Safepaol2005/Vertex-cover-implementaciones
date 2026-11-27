// gen.cpp — generadores de instancias sintéticas reproducibles (formato PACE .gr).
//   gen er   <n> <c> <seed>        Erdős–Rényi G(n, c/(n-1))
//   gen ba   <n> <m> <seed>        Barabási–Albert
//   gen band <n> <b> <p> <seed>    banda: arista i-j (|i-j|<=b) con prob. p  (pathwidth <= b)
//   gen trap <k>                   trampa bipartita de Johnson (1974) para el voraz
//   gen bip  <a> <b> <c> <seed>    bipartito aleatorio a+b vértices, grado medio c
#include <cstdio>
#include <random>
#include <set>
#include <string>
#include <vector>

int main(int argc, char** argv) {
    if (argc < 3) { std::fprintf(stderr, "uso: gen <er|ba|band|trap> ...\n"); return 2; }
    std::string t = argv[1];
    std::vector<std::pair<int, int>> E;
    int n = 0;
    if (t == "er") {
        n = std::stoi(argv[2]);
        double c = std::stod(argv[3]);
        std::mt19937_64 rng(std::stoull(argv[4]));
        std::uniform_real_distribution<double> U(0, 1);
        double p = c / (n - 1);
        for (int u = 0; u < n; ++u)
            for (int v = u + 1; v < n; ++v)
                if (U(rng) < p) E.push_back({u, v});
    } else if (t == "ba") {
        n = std::stoi(argv[2]);
        int m = std::stoi(argv[3]);
        std::mt19937_64 rng(std::stoull(argv[4]));
        std::vector<int> rep;
        std::vector<int> targets;
        for (int i = 0; i < m; ++i) targets.push_back(i);
        for (int v = m; v < n; ++v) {
            for (int x : targets) { E.push_back({x, v}); rep.push_back(x); rep.push_back(v); }
            std::set<int> s;
            while ((int)s.size() < m) s.insert(rep[rng() % rep.size()]);
            targets.assign(s.begin(), s.end());
        }
    } else if (t == "band") {
        n = std::stoi(argv[2]);
        int b = std::stoi(argv[3]);
        double p = std::stod(argv[4]);
        std::mt19937_64 rng(std::stoull(argv[5]));
        std::uniform_real_distribution<double> U(0, 1);
        for (int i = 0; i < n; ++i)
            for (int j = i + 1; j <= std::min(n - 1, i + b); ++j)
                if (U(rng) < p) E.push_back({i, j});
    } else if (t == "bip") {
        int a = std::stoi(argv[2]), b = std::stoi(argv[3]);
        double c = std::stod(argv[4]);
        std::mt19937_64 rng(std::stoull(argv[5]));
        std::uniform_real_distribution<double> U(0, 1);
        n = a + b;
        double p = c * n / (2.0 * a * b);
        for (int u = 0; u < a; ++u)
            for (int v = 0; v < b; ++v)
                if (U(rng) < p) E.push_back({u, a + v});
    } else if (t == "trap") {
        int k = std::stoi(argv[2]);
        int nR = 0;
        for (int i = 2; i <= k; ++i) nR += k / i;
        n = nR + k;  // R recibe los ids pequeños: el desempate favorece a R (caso adverso)
        int r = 0;
        for (int i = k; i >= 2; --i)
            for (int j = 0; j < k / i; ++j, ++r)
                for (int t2 = j * i; t2 < (j + 1) * i; ++t2) E.push_back({r, nR + t2});
    }
    std::printf("p td %d %zu\n", n, E.size());
    for (auto [u, v] : E) std::printf("%d %d\n", u + 1, v + 1);
    return 0;
}
