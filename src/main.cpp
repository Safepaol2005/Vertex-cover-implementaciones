// main.cpp — línea de comandos para correr un algoritmo sobre un grafo.
//
//   vcsolve <algoritmo> <archivo> [--time T] [--seed S] [--runs R] [--width W]
//           [--no-reduce] [--no-split] [--no-dp] [--no-mc] [--cover salida.txt]
//
// algoritmos: greedy | dp | rmis | dc | hybrid
// Imprime una línea CSV:  alg,file,n,m,size,time_s,status,extra
//   status = ok (óptimo si el algoritmo es exacto) | timeout | toowide | anytime | heur
#include <cstdio>
#include <cstring>
#include <iostream>

#include "divide_conquer.hpp"
#include "dp_td.hpp"
#include "graph.hpp"
#include "greedy_wh.hpp"
#include "hybrid.hpp"
#include "random_mis.hpp"

using namespace vc;

int main(int argc, char** argv) {
    if (argc < 3) {
        std::fprintf(stderr, "uso: %s <greedy|dp|rmis|dc|hybrid> <grafo> [opciones]\n", argv[0]);
        return 2;
    }
    std::string alg = argv[1], file = argv[2], cover_out;
    double T = 30;
    unsigned long long seed = 1;
    int runs = 1, width = 22;
    HybridParams hp;
    for (int i = 3; i < argc; ++i) {
        std::string a = argv[i];
        auto next = [&]() { return std::string(argv[++i]); };
        if (a == "--time") T = std::stod(next());
        else if (a == "--seed") seed = std::stoull(next());
        else if (a == "--runs") runs = std::stoi(next());
        else if (a == "--width") width = std::stoi(next());
        else if (a == "--cover") cover_out = next();
        else if (a == "--no-reduce") hp.use_reduce = false;
        else if (a == "--no-split") hp.use_split = false;
        else if (a == "--no-dp") hp.use_dp = false;
        else if (a == "--no-mc") hp.use_mc = false;
        else if (a == "--no-strong") hp.use_strong = false;
    }
    Graph g = read_graph(file);
    Timer t;
    std::vector<int> C;
    std::string status = "ok", extra;
    try {
        if (alg == "greedy") { C = greedy_wh(g); status = "heur"; }
        else if (alg == "rmis") { C = random_mis(g, runs, seed); status = "heur"; extra = "runs=" + std::to_string(runs); }
        else if (alg == "dp") {
            Deadline dl(T);
            if (!dp_vertex_cover(g, width, C, &dl)) status = "toowide";
            TreeDec td;
            min_degree_td(g, 1 << 30, td);
            extra = "width=" + std::to_string(td.width);
        } else if (alg == "dc") {
            Deadline dl(T);
            C = divide_and_conquer(g, dl);
        } else if (alg == "hybrid") {
            hp.time_limit = T;
            hp.seed = seed;
            HybridStats st;
            C = hybrid_vc(g, hp, st);
            if (!st.exact) status = "anytime";
            extra = "kernel_n=" + std::to_string(st.kernel_n) + ";comps=" + std::to_string(st.components) +
                    ";dp=" + std::to_string(st.dp_calls) + ";nodes=" + std::to_string(st.nodes) +
                    ";lb=" + std::to_string(st.lower_bound);
        } else {
            std::fprintf(stderr, "algoritmo desconocido: %s\n", alg.c_str());
            return 2;
        }
    } catch (Timeout&) {
        status = "timeout";
        C.clear();
    }
    double dt = t.seconds();
    bool has = (status != "timeout" && status != "toowide");
    if (has && !is_vertex_cover(g, C)) {
        std::fprintf(stderr, "ERROR: %s devolvió algo que no es cubierta en %s\n", alg.c_str(), file.c_str());
        return 1;
    }
    std::printf("%s,%s,%d,%lld,%s,%.4f,%s,%s\n", alg.c_str(), file.c_str(), g.n, g.m(),
                has ? std::to_string(C.size()).c_str() : "", dt, status.c_str(), extra.c_str());
    if (has && !cover_out.empty()) {
        FILE* f = std::fopen(cover_out.c_str(), "w");
        std::fprintf(f, "s vc %d %zu\n", g.n, C.size());
        for (int v : C) std::fprintf(f, "%d\n", v + 1);
        std::fclose(f);
    }
    return 0;
}
