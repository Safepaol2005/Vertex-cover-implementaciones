"""Genera las tablas (.tex) y figuras (.pdf) del paper a partir de results/*.csv."""
from __future__ import annotations

import csv
import os
from collections import defaultdict
from statistics import mean

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.join(os.path.dirname(__file__), "..")
RES = os.path.join(ROOT, "results")
TAB = os.path.join(ROOT, "paper", "tables")
FIG = os.path.join(ROOT, "paper", "figs")
os.makedirs(TAB, exist_ok=True)
os.makedirs(FIG, exist_ok=True)

# Paleta categórica validada (orden fijo; el color sigue al algoritmo, no a su rango)
COL = {"G-MaxDeg": "#2a78d6", "G-FastVC": "#eb6834", "R-MC(g=1,R=16)": "#1baf7a",
       "R-MIS": "#eda100", "2-Approx": "#e87ba4", "DP-TD": "#008300", "DC-WH": "#6250d6",
       "HyVC": "#e34948"}
MARK = {"G-MaxDeg": "o", "G-FastVC": "s", "R-MC(g=1,R=16)": "^", "R-MIS": "v", "2-Approx": "D",
        "DP-TD": "P", "DC-WH": "X", "HyVC": "*"}
LABEL = {"R-MC(g=1,R=16)": "R-MC"}
FAMS = ["ER-n200", "ER-n500", "PartialKTree", "Banded", "BA-n1000", "GreedyTrap", "ER-dense"]
FAMNAME = {"ER-n200": "Erdős–Rényi $n{=}200$", "ER-n500": "Erdős–Rényi $n{=}500$",
           "PartialKTree": "$k$-árbol parcial $n{=}300$", "Banded": "Banda $n{=}400$", "BA-n1000": "Barabási–Albert $n{=}1000$",
           "GreedyTrap": "Trampa de Johnson", "ER-dense": "Erdős–Rényi denso $n{=}60$"}
HEUR = ["G-MaxDeg", "G-FastVC", "2-Approx", "R-MC(g=1,R=16)", "R-MIS"]
EXACT = ["DP-TD", "DC-WH", "HyVC"]

plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.color": "#e6e6e3", "grid.linewidth": 0.6,
                     "axes.edgecolor": "#8a8a85", "lines.linewidth": 2, "legend.frameon": False})


def load(name):
    with open(os.path.join(RES, name)) as f:
        return list(csv.DictReader(f))


def lab(a):
    return LABEL.get(a, a)


def write(name, text):
    with open(os.path.join(TAB, name), "w") as f:
        f.write(text)


def main():
    rows = load("main.csv")
    by_inst = defaultdict(dict)
    meta = {}
    for r in rows:
        by_inst[r["instance"]][r["alg"]] = r
        meta[r["instance"]] = r
    insts = list(meta)
    write("n_instances.tex", f"{len(insts)}~")

    # --------------------------------------------------- familias
    fam_rows = []
    for fam in FAMS:
        I = [i for i in insts if meta[i]["family"] == fam]
        if not I:
            continue
        ns = sorted({int(meta[i]["n"]) for i in I})
        ms = [int(meta[i]["m"]) for i in I]
        tws = [int(meta[i]["tw_ub"]) for i in I]
        fam_rows.append(f"{FAMNAME[fam]} & {len(I)} & {ns[0]}" + (f"--{ns[-1]}" if len(ns) > 1 else "") +
                        f" & {min(ms)}--{max(ms)} & {min(tws)}--{max(tws)} \\\\")
    write("families.tex", "\\begin{table}[!htbp]\\centering\\small\n\\begin{tabular}{@{}lrrrr@{}}\\toprule\n"
          "Familia & Instancias & $n$ & $m$ & ancho (grado mín.) \\\\\\midrule\n" + "\n".join(fam_rows) +
          "\n\\bottomrule\\end{tabular}\n\\caption{Familias de instancias. La última columna es el ancho de la "
          "descomposición obtenida con la heurística de grado mínimo, una cota superior del treewidth.}"
          "\\label{tab:familias}\\end{table}\n")

    # --------------------------------------------------- calidad
    def ratio(i, a):
        r = by_inst[i][a]
        if r["size"] == "" or r["opt"] == "" or r["opt_proven"] != "True":
            return None
        return int(r["size"]) / int(r["opt"])

    algs_q = HEUR + ["HyVC"]
    qrows = []
    overall = defaultdict(list)
    for fam in FAMS:
        I = [i for i in insts if meta[i]["family"] == fam]
        cells = []
        for a in algs_q:
            vals = [ratio(i, a) for i in I]
            ok = [v for v in vals if v is not None]
            overall[a].extend(ok)
            if not ok:
                cells.append("---")
            else:
                s = f"{mean(ok):.3f}"
                if len(ok) < len(I):
                    s += f"$^{{{len(ok)}/{len(I)}}}$"
                cells.append(s)
        qrows.append(f"{FAMNAME[fam]} & " + " & ".join(cells) + " \\\\")
    cells = [f"{mean(overall[a]):.3f}" if overall[a] else "---" for a in algs_q]
    qrows.append("\\midrule Todas & " + " & ".join(cells) + " \\\\")
    write("quality.tex", "\\begin{table}[!htbp]\\centering\\scriptsize\\setlength{\\tabcolsep}{3.5pt}\n"
          "\\begin{tabular}{@{}l" + "r" * len(algs_q) + "@{}}\\toprule\n"
          "Familia & " + " & ".join(f"\\textsf{{{lab(a)}}}" for a in algs_q) + " \\\\\\midrule\n" +
          "\n".join(qrows) + "\n\\bottomrule\\end{tabular}\n"
          "\\caption{Razón media $|C|/\\OPT$ (1.000 = óptimo). El superíndice $a/b$ indica que el algoritmo solo "
          "tiene óptimo certificado en $a$ de las $b$ instancias; la media se calcula sobre esas $a$. "
          "Los exactos \\textsf{DP-TD} y \\textsf{DC-WH} se omiten: cuando terminan su razón es 1 (ver \\cref{tab:exactos}).}\\label{tab:calidad}\\end{table}\n")

    # --------------------------------------------------- exactos: resueltas y tiempo
    erows = []
    totals = {a: [0, 0] for a in EXACT}
    for fam in FAMS:
        I = [i for i in insts if meta[i]["family"] == fam]
        cells = []
        for a in EXACT:
            solved = [i for i in I if by_inst[i][a]["status"] == "ok"]
            totals[a][0] += len(solved)
            totals[a][1] += len(I)
            if solved:
                t = mean(float(by_inst[i][a]["time_s"]) for i in solved)
                cells.append(f"{len(solved)}/{len(I)} & {t:.2f}")
            else:
                cells.append(f"0/{len(I)} & ---")
        erows.append(f"{FAMNAME[fam]} & " + " & ".join(cells) + " \\\\")
    erows.append("\\midrule Total & " + " & ".join(f"{s}/{n} & " for s, n in totals.values()) + " \\\\")
    write("exact.tex", "\\begin{table}[!htbp]\\centering\\small\n\\begin{tabular}{@{}l" + "rr" * 3 + "@{}}\\toprule\n"
          " & \\multicolumn{2}{c}{\\textsf{DP-TD}} & \\multicolumn{2}{c}{\\textsf{DC-WH}} & "
          "\\multicolumn{2}{c}{\\hyvc} \\\\\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\\cmidrule(l){6-7}\n"
          "Familia & resueltas & $\\bar t$ (s) & resueltas & $\\bar t$ (s) & resueltas & $\\bar t$ (s) \\\\\\midrule\n" +
          "\n".join(erows) + "\n\\bottomrule\\end{tabular}\n"
          "\\caption{Algoritmos exactos: instancias resueltas con optimalidad probada (límite 20\\,s; "
          "\\textsf{DP-TD} se rinde si el ancho supera 16) y tiempo medio sobre las resueltas.}"
          "\\label{tab:exactos}\\end{table}\n")

    # --------------------------------------------------- composición de HyVC
    hrows = []
    for fam in FAMS:
        I = [i for i in insts if meta[i]["family"] == fam]
        h = [by_inst[i]["HyVC"] for i in I]
        red = mean(int(x["hy_reduced_n"]) / int(x["n"]) for x in h)
        dpc = sum(int(x["hy_dp_comp"]) for x in h)
        bbc = sum(int(x["hy_bb_comp"]) for x in h)
        nodes = mean(int(x["hy_nodes"]) for x in h)
        hrows.append(f"{FAMNAME[fam]} & {100 * red:.1f}\\% & {dpc} & {bbc} & {nodes:.0f} \\\\")
    write("hyvc_breakdown.tex", "\\begin{table}[!htbp]\\centering\\small\n\\begin{tabular}{@{}lrrrr@{}}\\toprule\n"
          "Familia & núcleo tras reducir & comp.\\ a PD & comp.\\ a B\\&R & nodos B\\&R (media) \\\\\\midrule\n" +
          "\n".join(hrows) + "\n\\bottomrule\\end{tabular}\n"
          "\\caption{Qué sección de \\hyvc{} hace el trabajo: fracción de vértices que sobrevive a las reducciones, "
          "componentes despachadas a programación dinámica o a branch \\& reduce, y nodos de búsqueda.}"
          "\\label{tab:hyvc}\\end{table}\n")

    # --------------------------------------------------- figura 1: razón vs c (ER n=200)
    I = [i for i in insts if meta[i]["family"] == "ER-n200"]
    cs = sorted({float(meta[i]["param"]) for i in I})
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.8))
    ax = axes[0]
    for a in HEUR + ["HyVC"]:
        ys = []
        for c in cs:
            v = [ratio(i, a) for i in I if float(meta[i]["param"]) == c]
            v = [x for x in v if x is not None]
            ys.append(mean(v) if v else float("nan"))
        ax.plot(cs, ys, color=COL[a], marker=MARK[a], markersize=5, label=lab(a))
    ax.set_xlabel("grado medio $c$")
    ax.set_ylabel("$|C|/\\mathrm{OPT}$")
    ax.set_title("Calidad (Erdős–Rényi, $n=200$)", fontsize=9)
    ax.axvline(2.718, color="#8a8a85", lw=0.8, ls="--")
    ax.text(2.78, ax.get_ylim()[1] * 0.98, "$c=e$", fontsize=7, color="#5f5f5a", va="top")
    ax = axes[1]
    for a in EXACT:
        ys = []
        for c in cs:
            v = [float(by_inst[i][a]["time_s"]) if by_inst[i][a]["status"] == "ok" else 20.0
                 for i in I if float(meta[i]["param"]) == c]
            ys.append(mean(v))
        ax.plot(cs, ys, color=COL[a], marker=MARK[a], markersize=5, label=lab(a))
    ax.set_yscale("log")
    ax.axhline(20, color="#8a8a85", lw=0.8, ls=":")
    ax.text(cs[0], 21.5, "límite 20 s", fontsize=7, color="#5f5f5a")
    ax.set_xlabel("grado medio $c$")
    ax.set_ylabel("tiempo (s, escala log)")
    ax.set_title("Tiempo de los exactos", fontsize=9)
    h1, l1 = axes[0].get_legend_handles_labels()
    h2, l2 = axes[1].get_legend_handles_labels()
    seen = {}
    for h, l in list(zip(h1, l1)) + list(zip(h2, l2)):
        seen.setdefault(l, h)
    fig.legend(seen.values(), seen.keys(), loc="lower center", ncol=8, fontsize=7.5,
               bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    fig.savefig(os.path.join(FIG, "er_quality_time.pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(FIG, "er_quality_time.png"), dpi=160, bbox_inches="tight")
    plt.close(fig)

    summary = {"overall_ratio": {a: (mean(v) if v else None) for a, v in overall.items()},
               "exact_totals": totals}

    # --------------------------------------------------- ablación
    if os.path.exists(os.path.join(RES, "ablation.csv")):
        ab = load("ablation.csv")
        variants = list(dict.fromkeys(r["variant"] for r in ab))
        fams = ["ER-n200", "Banded", "BA-n1000"]
        fig, ax = plt.subplots(figsize=(6.6, 2.6))
        width = 0.8 / len(variants)
        vcol = ["#e34948", "#2a78d6", "#008300", "#6250d6", "#1baf7a"]
        abl_rows = []
        for j, v in enumerate(variants):
            ys = []
            for fam in fams:
                R = [r for r in ab if r["family"] == fam and r["variant"] == v]
                ys.append(mean(float(r["time_s"]) for r in R))
                n_ex = sum(r["exact"] == "True" for r in R)
                abl_rows.append((fam, v, mean(float(r["time_s"]) for r in R), n_ex, len(R),
                                 mean(int(r["nodes"]) for r in R)))
            xs = [k + (j - (len(variants) - 1) / 2) * width for k in range(len(fams))]
            ax.bar(xs, ys, width=width * 0.92, color=vcol[j % len(vcol)], label=v)
        ax.set_xticks(range(len(fams)))
        ax.set_xticklabels([FAMNAME[f] for f in fams])
        ax.set_yscale("log")
        ax.set_ylabel("tiempo medio (s, log)")
        ax.legend(ncol=5, fontsize=7, loc="upper center", bbox_to_anchor=(0.5, 1.18))
        fig.tight_layout()
        fig.savefig(os.path.join(FIG, "ablation.pdf"), bbox_inches="tight")
        fig.savefig(os.path.join(FIG, "ablation.png"), dpi=160, bbox_inches="tight")
        plt.close(fig)
        lines = []
        for fam in fams:
            for (f, v, t, ne, nn, nd) in abl_rows:
                if f == fam:
                    lines.append(f"{FAMNAME[f] if v == variants[0] else ''} & {v} & {t:.2f} & {ne}/{nn} & {nd:.0f} \\\\")
            lines.append("\\addlinespace")
        write("ablation.tex", "\\begin{table}[!htbp]\\centering\\small\n\\begin{tabular}{@{}llrrr@{}}\\toprule\n"
              "Familia & Variante & $\\bar t$ (s) & exactas & nodos \\\\\\midrule\n" + "\n".join(lines[:-1]) +
              "\n\\bottomrule\\end{tabular}\n\\caption{Ablación: se apaga una sección de \\hyvc{} a la vez "
              "(límite 10\\,s). Si una variante no termina, aporta 10\\,s a la media.}"
              "\\label{tab:ablacion}\\end{table}\n")
        summary["ablation"] = abl_rows

    # --------------------------------------------------- barrido de gamma
    if os.path.exists(os.path.join(RES, "gamma.csv")):
        gm = load("gamma.csv")
        names = list(dict.fromkeys(r["instance"] for r in gm))
        fig, axes = plt.subplots(1, len(names), figsize=(7.2, 2.2), sharey=False)
        for ax, nm in zip(axes, names):
            R = [r for r in gm if r["instance"] == nm]
            opt = float(R[0]["opt"])
            for runs, colr, mk in (("1", "#2a78d6", "o"), ("16", "#e34948", "s")):
                S = [r for r in R if r["runs"] == runs]
                ax.plot([float(r["gamma"]) + 0.5 for r in S], [int(r["size"]) / opt for r in S],
                        color=colr, marker=mk, markersize=4, label=f"R={runs}")
            ax.set_xscale("log")
            ax.set_xticks([0.5, 1.5, 2.5, 4.5, 8.5, 32.5])
            ax.set_xticklabels(["0", "1", "2", "4", "8", "32"], fontsize=7)
            ax.set_title(nm.replace("_", " "), fontsize=8)
            ax.set_xlabel("$\\gamma$", fontsize=8)
        axes[0].set_ylabel("$|C|/\\mathrm{OPT}$")
        axes[0].legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(os.path.join(FIG, "gamma.pdf"), bbox_inches="tight")
        fig.savefig(os.path.join(FIG, "gamma.png"), dpi=160, bbox_inches="tight")
        plt.close(fig)
        summary["gamma"] = [(r["instance"], r["gamma"], r["runs"], r["size"], r["opt"]) for r in gm]

    import json
    with open(os.path.join(RES, "summary.json"), "w") as f:
        json.dump(summary, f, indent=1, default=str)
    print(json.dumps(summary, indent=1, default=str)[:4000])


if __name__ == "__main__":
    main()
