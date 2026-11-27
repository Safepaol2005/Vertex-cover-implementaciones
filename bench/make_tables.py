"""Analiza results/*_raw.csv y genera las tablas (.tex), figuras (.pdf/.png) y un resumen JSON.

Uso: python3 bench/make_tables.py
"""
from __future__ import annotations

import csv
import json
import os
import re
from collections import defaultdict
from statistics import mean, median

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.join(os.path.dirname(__file__), "..")
RES = os.path.join(ROOT, "results")
TAB = os.path.join(ROOT, "paper", "tables")
FIG = os.path.join(ROOT, "paper", "figs")
os.makedirs(TAB, exist_ok=True)
os.makedirs(FIG, exist_ok=True)
LIMIT = 10.0

COL = {"Voraz-WH": "#2a78d6", "MIS (R=1)": "#eda100", "MIS (R=64)": "#eb6834", "PD-TD": "#008300",
       "DyV-WH": "#6250d6", "HyVC": "#e34948", "WeGotYouCovered": "#1d1d1b"}
MARK = {"Voraz-WH": "o", "MIS (R=1)": "v", "MIS (R=64)": "^", "PD-TD": "P", "DyV-WH": "X", "HyVC": "*",
        "WeGotYouCovered": "s"}
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.color": "#e6e6e3", "grid.linewidth": 0.6,
                     "axes.edgecolor": "#8a8a85", "lines.linewidth": 2, "legend.frameon": False})


def alg_name(alg, extra):
    if alg == "greedy":
        return "Voraz-WH"
    if alg == "rmis":
        return "MIS (R=1)" if "runs=1" in extra.split(";") or extra == "runs=1" else "MIS (R=64)"
    return {"dp": "PD-TD", "dc": "DyV-WH", "hybrid": "HyVC"}[alg]


def load(ds):
    path = os.path.join(RES, f"{ds}_raw.csv")
    rows = defaultdict(dict)
    meta = {}
    if not os.path.exists(path):
        return rows, meta
    for line in open(path):
        p = line.rstrip("\n").split(",")
        if len(p) < 7:
            continue
        alg, f, n, m, size, t, status = p[:7]
        extra = ",".join(p[7:])
        inst = os.path.basename(f).replace(".gr", "")
        name = alg_name(alg, extra)
        rows[inst][name] = {"size": int(size) if size else None, "time": float(t) if t else None,
                            "status": status, "extra": extra}
        if n:
            meta[inst] = {"n": int(n), "m": int(m), "file": f}
    # reejecuciones de HyVC (versión final) sobrescriben las del archivo crudo
    hy = os.path.join(RES, f"{ds}_hy.csv")
    if os.path.exists(hy):
        for line in open(hy):
            p = line.rstrip("\n").split(",")
            if len(p) < 7:
                continue
            alg, f, n, m, size, t, status = p[:7]
            inst = os.path.basename(f).replace(".gr", "")
            rows[inst]["HyVC"] = {"size": int(size) if size else None, "time": float(t) if t else None,
                                  "status": status, "extra": ",".join(p[7:])}
    return rows, meta


def wgyc():
    ref = {}
    p = os.path.join(ROOT, "data", "wgyc_reference.csv")
    for r in csv.DictReader(open(p)):
        inst = "vc-exact_" + r["instance"]
        ref[inst] = {"vc": int(r["wgyc_vc"]) if r["wgyc_vc"] else None,
                     "time": float(r["wgyc_full_time_s"]),
                     "certified": "checkmark" in r["wgyc_solved_in_table"]}
    return ref


def cpsat(ds):
    out = {}
    p = os.path.join(RES, f"{ds}_cpsat.csv")
    if os.path.exists(p):
        for r in csv.DictReader(open(p)):
            inst = os.path.basename(r["file"]).replace(".gr", "")
            out[inst] = {"size": int(r["cpsat_size"]) if r["cpsat_size"] else None,
                         "opt": r["cpsat_optimal"] == "True",
                         "bound": int(r["cpsat_bound"]) if r["cpsat_bound"] else None}
    return out


EXACT = ["PD-TD", "DyV-WH", "HyVC"]
HEUR = ["Voraz-WH", "MIS (R=1)", "MIS (R=64)"]


def reference(ds, rows, W, CP):
    """opt[inst] = (valor, fuente) si está certificado; best[inst] = mejor conocido."""
    opt, best, issues = {}, {}, []
    for inst, R in rows.items():
        vals = []
        for a in EXACT:
            r = R.get(a)
            if r and r["status"] == "ok" and r["size"] is not None:
                vals.append((r["size"], a))
        if W and inst in W and W[inst]["certified"] and W[inst]["vc"]:
            vals.append((W[inst]["vc"], "WeGotYouCovered"))
        if inst in CP and CP[inst]["opt"]:
            vals.append((CP[inst]["size"], "CP-SAT"))
        if len({v for v, _ in vals}) > 1:
            issues.append((inst, vals))
        if vals:
            opt[inst] = (min(v for v, _ in vals), vals[0][1])
        sizes = [r["size"] for r in R.values() if r["size"] is not None]
        if inst in CP and CP[inst]["size"]:
            sizes.append(CP[inst]["size"])
        if sizes:
            best[inst] = min(sizes)
    return opt, best, issues


def pct(x):
    return f"{100 * x:.1f}\\,\\%"


def write(name, text):
    with open(os.path.join(TAB, name), "w") as f:
        f.write(text)


def main():
    W = wgyc()
    summary = {}
    data = {}
    for ds in ("pace2019", "pace2017", "synthetic"):
        rows, meta = load(ds)
        CP = cpsat(ds)
        opt, best, issues = reference(ds, rows, W if ds == "pace2019" else None, CP)
        data[ds] = (rows, meta, opt, best, CP)
        summary[ds] = {"instances": len(rows), "with_opt": len(opt), "inconsistencies": issues}
        # instancias sin óptimo certificado (para CP-SAT)
        with open(os.path.join(RES, f"{ds}_need_ref.txt"), "w") as f:
            for inst in sorted(rows):
                if inst not in opt and inst in meta:
                    f.write(meta[inst]["file"] + "\n")

    # ---------------------------------------------------------- tabla: exactos resueltos
    lines = []
    for ds, label in (("pace2019", "PACE 2019"), ("pace2017", "PACE 2017"), ("synthetic", "Sintéticas")):
        rows, meta, opt, best, CP = data[ds]
        if not rows:
            continue
        N = len(rows)
        cells = []
        for a in EXACT:
            ok = [i for i in rows if rows[i].get(a, {}).get("status") == "ok"]
            t = median(rows[i][a]["time"] for i in ok) if ok else None
            cells.append(f"{len(ok)} & {t:.2f}" if ok else "0 & ---")
            summary[ds][f"solved_{a}"] = len(ok)
        if ds == "pace2019":
            ok = [i for i in rows if i in W and W[i]["time"] <= LIMIT]
            okall = [i for i in rows if i in W and W[i]["time"] < 1800]
            t = median(W[i]["time"] for i in ok)
            cells.append(f"{len(ok)} & {t:.2f}")
            summary[ds]["solved_WGYC_10s"] = len(ok)
            summary[ds]["solved_WGYC_1800s"] = len(okall)
        else:
            cells.append("\\multicolumn{2}{c}{---}")
        cpok = sum(1 for i in rows if i in CP and CP[i]["opt"])
        lines.append(f"{label} ({N}) & " + " & ".join(cells) + " \\\\")
    write("exact_solved.tex",
          "\\begin{table}[!htbp]\\centering\\small\n\\begin{tabular}{@{}l" + "rr" * 4 + "@{}}\\toprule\n"
          " & \\multicolumn{2}{c}{\\alg{PD-TD}} & \\multicolumn{2}{c}{\\alg{DyV-WH}} & \\multicolumn{2}{c}{\\hyvc} & "
          "\\multicolumn{2}{c}{\\textsc{WeGotYouCovered}$^\\dagger$} \\\\"
          "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}\\cmidrule(l){8-9}\n"
          "Conjunto & res. & $\\tilde t$ (s) & res. & $\\tilde t$ (s) & res. & $\\tilde t$ (s) & res. & $\\tilde t$ (s) \\\\\\midrule\n"
          + "\n".join(lines) +
          "\n\\bottomrule\\end{tabular}\n\\caption{Instancias resueltas con optimalidad probada en 10\\,s y mediana "
          "del tiempo sobre las resueltas. $^\\dagger$Tiempos publicados por sus autores~\\citep{hespe2020} en su "
          "propio hardware; se cuentan las instancias que resolvió en 10\\,s o menos.}\\label{tab:exactos}\\end{table}\n")

    # ---------------------------------------------------------- tabla: calidad heurística
    qlines = []
    for ds, label in (("pace2019", "PACE 2019"), ("pace2017", "PACE 2017"), ("synthetic", "Sintéticas")):
        rows, meta, opt, best, CP = data[ds]
        if not rows:
            continue
        cells = []
        for a in HEUR + ["HyVC"]:
            rs = [rows[i][a]["size"] / opt[i][0] for i in opt if rows[i].get(a, {}).get("size") and opt[i][0] > 0]
            hit = sum(1 for i in opt if rows[i].get(a, {}).get("size") == opt[i][0])
            cells.append(f"{mean(rs):.3f} & {hit}" if rs else "--- & ---")
            summary[ds][f"ratio_{a}"] = mean(rs) if rs else None
            summary[ds][f"hits_{a}"] = hit
        qlines.append(f"{label} ({len(opt)}) & " + " & ".join(cells) + " \\\\")
    write("heur_quality.tex",
          "\\begin{table}[!htbp]\\centering\\small\\setlength{\\tabcolsep}{4pt}\n\\begin{tabular}{@{}l" + "rr" * 4 + "@{}}\\toprule\n"
          " & \\multicolumn{2}{c}{\\alg{Voraz-WH}} & \\multicolumn{2}{c}{MIS ($R{=}1$)} & \\multicolumn{2}{c}{MIS ($R{=}64$)} & "
          "\\multicolumn{2}{c}{\\hyvc} \\\\\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}\\cmidrule(l){8-9}\n"
          "Conjunto (con óptimo) & $\\overline{|C|/\\tau}$ & $=\\tau$ & $\\overline{|C|/\\tau}$ & $=\\tau$ & "
          "$\\overline{|C|/\\tau}$ & $=\\tau$ & $\\overline{|C|/\\tau}$ & $=\\tau$ \\\\\\midrule\n" + "\n".join(qlines) +
          "\n\\bottomrule\\end{tabular}\n\\caption{Calidad de las cubiertas en las instancias con óptimo certificado: "
          "razón media $|C|/\\tau$ y número de instancias donde el algoritmo encontró el óptimo.}"
          "\\label{tab:calidad}\\end{table}\n")

    # ---------------------------------------------------------- sintéticas por familia
    rows, meta, opt, best, CP = data["synthetic"]
    fams = [("er1000_c1", "ER $c{=}1$"), ("er1000_c2", "ER $c{=}2$"), ("er1000_c3", "ER $c{=}3$"),
            ("er1000_c4", "ER $c{=}4$"), ("er1000_c6", "ER $c{=}6$"), ("ba2000_m2", "BA $m{=}2$"),
            ("ba2000_m3", "BA $m{=}3$"), ("ba2000_m5", "BA $m{=}5$"), ("band1000_b8", "Banda $b{=}8$"),
            ("band1000_b12", "Banda $b{=}12$"), ("band1000_b16", "Banda $b{=}16$"),
            ("bip1000_c3", "Bipartito $c{=}3$"), ("bip1000_c6", "Bipartito $c{=}6$"), ("trap", "Trampa")]
    flines = []
    fam_summary = {}
    for key, label in fams:
        I = sorted(i for i in rows if i.startswith(key))
        if not I:
            continue
        cells = []
        for a in ["Voraz-WH", "MIS (R=64)"]:
            rs = [rows[i][a]["size"] / (opt[i][0] if i in opt else best[i]) for i in I if rows[i].get(a, {}).get("size")]
            cells.append(f"{mean(rs):.3f}" if rs else "---")
        for a in EXACT:
            ok = [i for i in I if rows[i].get(a, {}).get("status") == "ok"]
            cells.append(f"{len(ok)}/{len(I)}")
        ex = [i for i in I if rows[i].get("HyVC", {}).get("status") == "ok"]
        cells.append(f"{median(rows[i]['HyVC']['time'] for i in ex):.2f}" if ex else "---")
        ks = [int(re.search(r"kernel_n=(\d+)", rows[i]["HyVC"]["extra"]).group(1)) / meta[i]["n"]
              for i in I if "kernel_n=" in rows[i].get("HyVC", {}).get("extra", "")]
        kern = mean(ks) if ks else float("nan")
        cells.append(f"{100 * kern:.0f}\\,\\%")
        nopt = sum(1 for i in I if i in opt)
        flines.append(f"{label} & " + " & ".join(cells) + (" \\\\" if nopt == len(I) else f"$^{{*}}$ \\\\"))
        fam_summary[key] = {"cells": cells, "with_opt": nopt, "n": len(I)}
    summary["synthetic_families"] = fam_summary
    write("synthetic.tex",
          "\\begin{table}[!htbp]\\centering\\small\n\\begin{tabular}{@{}lrrrrrrr@{}}\\toprule\n"
          " & \\multicolumn{2}{c}{$|C|/\\tau$} & \\multicolumn{3}{c}{exactas en 10\\,s} & & \\\\"
          "\\cmidrule(lr){2-3}\\cmidrule(lr){4-6}\n"
          "Familia & Voraz & MIS ($R{=}64$) & PD & DyV & \\hyvc & $\\tilde t_{\\hyvc}$ (s) & núcleo \\\\\\midrule\n"
          + "\n".join(flines) +
          "\n\\bottomrule\\end{tabular}\n\\caption{Instancias sintéticas por familia. ``Núcleo'' es la fracción de vértices "
          "que sobrevive a las reducciones de \\hyvc. $^{*}$Familia con alguna instancia sin óptimo certificado: la razón "
          "se calcula contra la mejor cubierta conocida.}\\label{tab:sinteticas}\\end{table}\n")

    # ---------------------------------------------------------- textos cortos
    rows, meta, opt, best, CP = data["synthetic"]
    t = "trap_k160"
    if t in rows:
        g = rows[t]["Voraz-WH"]["size"]
        write("trap_greedy.tex", f"{g} vértices frente a $k=160$, una razón de {g / 160:.2f}")
    rows, meta, opt, best, CP = data["pace2019"]
    rs = [rows[i]["Voraz-WH"]["size"] / opt[i][0] - 1 for i in opt if rows[i].get("Voraz-WH")]
    if rs:
        write("greedy_avg_gap.tex", f"un {100 * mean(rs):.1f}\\,\\%")
    # anchos PACE 2017
    rows, meta, opt, best, CP = data["pace2017"]
    tw = {r["instance"]: int(r["treewidth"]) for r in csv.DictReader(open(os.path.join(ROOT, "data", "pace2017_tw.csv")))}
    pairs = []
    for i in rows:
        e = rows[i].get("PD-TD", {}).get("extra", "")
        m_ = re.search(r"width=(\d+)", e)
        if m_ and i in tw:
            pairs.append((tw[i], int(m_.group(1))))
    if pairs:
        ratio = [w / t for t, w in pairs if t > 0]
        exact_w = sum(1 for t, w in pairs if w == t)
        write("pace2017_width_gap.tex",
              f"de {median(w for _, w in pairs):.0f} frente a un ancho exacto mediano de {median(t for t, _ in pairs):.0f}: "
              f"en promedio la heurística sobreestima el ancho un {100 * (mean(ratio) - 1):.0f}\\,\\% y lo acierta en "
              f"{exact_w} de {len(pairs)} instancias")
        summary["pace2017_width"] = {"median_heur": median(w for _, w in pairs), "median_true": median(t for t, _ in pairs),
                                     "mean_overest": mean(ratio) - 1, "exact": exact_w, "n": len(pairs),
                                     "true_le20": sum(1 for t, _ in pairs if t <= 20),
                                     "heur_le22": sum(1 for _, w in pairs if w <= 22)}
        fig, ax = plt.subplots(figsize=(3.6, 3.0))
        ax.scatter([t for t, _ in pairs], [w for _, w in pairs], s=10, color="#2a78d6", alpha=0.7, label="instancia")
        mx = max(max(t for t, _ in pairs), max(w for _, w in pairs))
        ax.plot([1, mx], [1, mx], color="#8a8a85", lw=1, ls="--", label="ancho exacto")
        ax.axhline(22, color="#e34948", lw=1, ls=":", label="límite de la PD (22)")
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel("ancho de árbol exacto"); ax.set_ylabel("ancho por grado mínimo")
        ax.legend(fontsize=7, loc="upper left")
        fig.tight_layout()
        fig.savefig(os.path.join(FIG, "pace2017_width.pdf")); fig.savefig(os.path.join(FIG, "pace2017_width.png"), dpi=160)
        plt.close(fig)

    # ---------------------------------------------------------- figura: instancias resueltas vs tiempo (PACE 2019)
    rows, meta, opt, best, CP = data["pace2019"]
    if rows:
        fig, ax = plt.subplots(figsize=(6.4, 3.2))
        curves = {}
        for a in EXACT:
            ts = sorted(rows[i][a]["time"] for i in rows if rows[i].get(a, {}).get("status") == "ok")
            curves[a] = ts
        curves["WeGotYouCovered"] = sorted(W[i]["time"] for i in rows if i in W and W[i]["time"] < 1800)
        for a, ts in curves.items():
            xs = [max(t, 1e-3) for t in ts]
            ax.step(xs, range(1, len(xs) + 1), where="post", color=COL[a], label=f"{a} ({sum(t <= LIMIT for t in ts)} en 10 s)")
        ax.axvline(LIMIT, color="#8a8a85", lw=0.8, ls="--")
        ax.set_xscale("log")
        ax.set_xlabel("tiempo (s, log)")
        ax.set_ylabel("instancias resueltas")
        ax.set_title("PACE 2019: instancias resueltas en función del tiempo", fontsize=9)
        ax.legend(fontsize=7, loc="upper left")
        fig.tight_layout()
        fig.savefig(os.path.join(FIG, "pace2019_cactus.pdf")); fig.savefig(os.path.join(FIG, "pace2019_cactus.png"), dpi=160)
        plt.close(fig)
        # HyVC anytime gap en las no resueltas
        gaps = []
        for i in rows:
            r = rows[i].get("HyVC")
            if r and r["status"] == "anytime" and r["size"]:
                ref = opt[i][0] if i in opt else best[i]
                gaps.append(r["size"] / ref - 1)
        summary["pace2019"]["hyvc_anytime_n"] = len(gaps)
        summary["pace2019"]["hyvc_anytime_gap_mean"] = mean(gaps) if gaps else None
        summary["pace2019"]["hyvc_anytime_gap_max"] = max(gaps) if gaps else None
        # instancias que HyVC resuelve y WGYC no en 10s, y viceversa
        hy = {i for i in rows if rows[i].get("HyVC", {}).get("status") == "ok"}
        wg = {i for i in rows if i in W and W[i]["time"] <= LIMIT}
        summary["pace2019"]["hy_not_wg"] = sorted(hy - wg)
        summary["pace2019"]["wg_not_hy"] = len(wg - hy)
        # clasificación de las no resueltas por HyVC
        uns = [i for i in rows if "HyVC" in rows[i] and rows[i]["HyVC"]["status"] != "ok"]
        def kn(i):
            m_ = re.search(r"kernel_n=(\d+)", rows[i].get("HyVC", {}).get("extra", ""))
            return int(m_.group(1)) if m_ else None
        summary["pace2019"]["unsolved_kernel"] = sorted((i, meta[i]["n"], kn(i)) for i in uns)

    # ---------------------------------------------------------- ablación
    abp = os.path.join(RES, "ablation_raw.csv")
    if os.path.exists(abp):
        ab = defaultdict(dict)
        for line in open(abp):
            p = line.rstrip("\n").split(",")
            if len(p) < 8:
                continue
            flag, f, status = p[0], p[2], p[7]
            ds = "pace2019" if "pace2019" in f else ("pace2017" if "pace2017" in f else "synthetic")
            ab[(ds, flag)][os.path.basename(f)] = (status, float(p[6]) if p[6] else 10.0)
        names = [("full", "\\hyvc{} completo"), ("--no-reduce", "sin reducciones"), ("--no-strong", "sin R4$'$ y R5"),
                 ("--no-dp", "sin PD"), ("--no-split", "sin división"), ("--no-mc", "sin MIS aleatorio")]
        sets = [("pace2017", "PACE 2017 (200)"), ("synthetic", "Sintéticas (35)"), ("pace2019", "PACE 2019$^\\ddagger$")]
        lines = []
        abl = {}
        for flag, label in names:
            cells = []
            for ds, _ in sets:
                rows_ = data[ds][0]
                if flag == "full":
                    files = set(ab[(ds, "--no-mc")].keys())
                    res = {fn: (rows_[fn[:-3]]["HyVC"]["status"], rows_[fn[:-3]]["HyVC"]["time"]) for fn in files}
                else:
                    res = ab[(ds, flag)]
                ok = sum(1 for st, _ in res.values() if st == "ok")
                tt = sum(t if st == "ok" else LIMIT for st, t in res.values())
                cells.append(f"{ok} & {tt:.0f}")
                abl[(ds, flag)] = (ok, tt, len(res))
            lines.append(f"{label} & " + " & ".join(cells) + " \\\\")
        summary["ablation"] = {f"{k[0]}|{k[1]}": v for k, v in abl.items()}
        write("ablation.tex",
              "\\begin{table}[!htbp]\\centering\\small\n\\begin{tabular}{@{}lrrrrrr@{}}\\toprule\n"
              " & \\multicolumn{2}{c}{PACE 2017 (200)} & \\multicolumn{2}{c}{Sintéticas (35)} & \\multicolumn{2}{c}{PACE 2019$^\\ddagger$} \\\\"
              "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\\cmidrule(l){6-7}\n"
              "Variante & res. & $\\sum t$ (s) & res. & $\\sum t$ (s) & res. & $\\sum t$ (s) \\\\\\midrule\n" + "\n".join(lines) +
              "\n\\bottomrule\\end{tabular}\n\\caption{Ablación: se apaga una fase de \\hyvc{} a la vez. ``res.'' es el número de "
              "instancias resueltas en 10\\,s y $\\sum t$ el tiempo total (una instancia no resuelta cuenta 10\\,s). "
              "$^\\ddagger$Las 30 instancias del PACE 2019 que \\hyvc{} completo resuelve.}\\label{tab:ablacion}\\end{table}\n")

    with open(os.path.join(RES, "summary.json"), "w") as f:
        json.dump(summary, f, indent=1, default=str)
    print(json.dumps({k: {kk: vv for kk, vv in v.items() if kk not in ("unsolved_kernel", "inconsistencies")} if isinstance(v, dict) else v
                      for k, v in summary.items()}, indent=1, default=str)[:5000])
    for ds in ("pace2019", "pace2017", "synthetic"):
        if summary[ds]["inconsistencies"]:
            print("INCONSISTENCIAS", ds, summary[ds]["inconsistencies"][:10])


if __name__ == "__main__":
    main()
