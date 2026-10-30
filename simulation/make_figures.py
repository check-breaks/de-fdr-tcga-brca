#!/usr/bin/env python3
"""Publication-style figures. Run from the repo root:  python simulation/make_figures.py [results_dir=results]
Inputs (any missing input is skipped with a message):
  results/simulation_results_raw.csv            (evaluate.py)
  results/simulation_unbalanced_summary.csv     (evaluate_unbalanced.py)
  results/tcga/perm/BRCA1_mut_iter*.csv         (tcga/04_permutation_null.R)
  results/tcga/n_sig_summary.csv, results/tcga/wilcoxon_check_summary.csv   (real BRCA1 result, for the marker)
Outputs: results/figures/fdr_power_v2.png, fdr_unbalanced_v2.png, permutation_null.png, sva_null.png
"""
import sys, glob, pathlib, textwrap
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

R = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else "results"); FIG = R / "figures"; FIG.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
ORDER = ["deseq2", "edger", "limmavoom", "wilcoxon", "ttest", "nbglm_scratch"]
NICE = {"deseq2": "DESeq2", "edger": "edgeR (QL)", "limmavoom": "limma-voom", "wilcoxon": "Wilcoxon", "ttest": "t-test", "nbglm_scratch": "scratch NB GLM"}
COL = dict(zip(ORDER, ["#1f77b4", "#ff7f0e", "#2ca02c", "#8c564b", "#9467bd", "#d62728"]))
cap = lambda fig, t: fig.text(0.01, 0.005, textwrap.fill(t, 135), fontsize=8, va="bottom")
se_ = lambda x: x.std(ddof=1) / np.sqrt(len(x)) if len(x) > 1 else 0.0


def fig_fdr_power():
    f = R / "simulation_results_raw.csv"
    if not f.exists(): return print("skip fdr_power: missing", f)
    d = pd.read_csv(f)
    s = d.groupby(["method", "n_per_group"]).agg(fdr=("fdr", "mean"), fdr_se=("fdr", se_), power=("power", "mean"),
            power_se=("power", se_), reps=("fdr", "size"), any_zero=("n_called", lambda x: bool((x == 0).any()))).reset_index()
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
    for m in [m for m in ORDER if m in set(s.method)]:
        g = s[s.method == m].sort_values("n_per_group")
        for a, y, e in ((ax[0], "fdr", "fdr_se"), (ax[1], "power", "power_se")):
            a.plot(g.n_per_group, g[y], c=COL[m], lw=1.5, label=NICE[m])
            a.errorbar(g.n_per_group, g[y], yerr=[np.minimum(g[e], g[y]), g[e]], fmt="none", ecolor=COL[m], capsize=3, lw=1)
            for _, r in g.iterrows():                       # hollow marker = some replicate called nothing (FDR counted as 0)
                hollow = r.any_zero and y == "fdr"
                a.plot(r.n_per_group, r[y], "o", c=COL[m], mfc="white" if hollow else COL[m], ms=6, zorder=3)
    ax[0].axhline(0.05, ls="--", c="k", lw=1); ax[0].set_title("Realized FDR at nominal 5% (BH)"); ax[0].set_ylabel("FDR"); ax[0].set_ylim(bottom=-0.005)
    ax[1].set_title("Power"); ax[1].set_ylabel("true positive rate")
    for a in ax: a.set_xlabel("samples per group"); a.set_xticks(sorted(s.n_per_group.unique()))
    ax[1].legend(frameon=False, fontsize=8, loc="lower right")
    n = int(s.reps.max())
    cap(fig, f"{n} simulated datasets per point; error bars = ±1 SE (clipped at 0). Hollow marker: at least one dataset called no genes " "(its FDR counted as 0, so low FDR there is not good calibration).")
    plt.tight_layout(rect=(0, 0.07, 1, 1)); plt.savefig(FIG / "fdr_power_v2.png", dpi=170); plt.close()


def fig_unbalanced():
    f = R / "simulation_unbalanced_summary.csv"
    if not f.exists(): return print("skip unbalanced: missing", f)
    d = pd.read_csv(f); d["se_fdr"] = d.fdr_sd / np.sqrt(d.reps)
    d["base"] = d.method.str.rsplit("_", n=1).str[0]; d["design"] = d.method.str.rsplit("_", n=1).str[1]
    rows = [(b, g) for b in ["deseq2", "edger", "limmavoom", "wilcoxon", "ttest"] for g in ["naive", "adjusted"] if ((d.base == b) & (d.design == g)).any()]
    fig, ax = plt.subplots(2, 2, figsize=(11, 8), sharey=True, sharex="col")
    for i, sc in enumerate(["mixed", "confounded"]):
        for j, (y, lab) in enumerate((("fdr", "Realized FDR (nominal 5%)"), ("power", "Power"))):
            a = ax[i, j]
            for k, (nc, colr) in enumerate(((7, "#9ecae1"), (11, "#08519c"))):
                vals, errs = [], []
                for b, g in rows:
                    r = d[(d.scenario == sc) & (d.n_case == nc) & (d.base == b) & (d.design == g)]
                    vals.append(r[y].iloc[0] if len(r) else np.nan)
                    errs.append(r["se_fdr"].iloc[0] if (len(r) and y == "fdr") else 0)
                pos = np.arange(len(rows)) + (k - 0.5) * 0.38
                a.barh(pos, vals, height=0.36, color=colr, label=f"{nc} cases" if (i, j) == (0, 0) else None,
                       xerr=[np.minimum(errs, vals), errs] if y == "fdr" else None, error_kw=dict(lw=1, capsize=2))
            a.set_yticks(range(len(rows))); a.set_yticklabels([f"{NICE[b]} · {g}" for b, g in rows]); a.invert_yaxis()
            if y == "fdr": a.axvline(0.05, ls="--", c="k", lw=1)
            a.set_title(f"{sc}: {lab}", fontsize=10); a.set_xlabel(lab.split(" (")[0])
    for a in ax.flat: a.tick_params(labelbottom=True)
    ax[0, 0].legend(frameon=False, loc="lower right", title="cases vs ~300 controls")
    cap(fig, f"Mean over {int(d.reps.max())} simulated datasets per bar; FDR error bars = ±1 SE. 'confounded': cases enriched in one latent subtype " "(like 7 of 11 BRCA1-mutant tumors being basal-like).")
    plt.tight_layout(rect=(0, 0.06, 1, 1)); plt.savefig(FIG / "fdr_unbalanced_v2.png", dpi=170); plt.close()


def fig_perm():
    files = sorted(glob.glob(str(R / "tcga/perm/BRCA1_mut_iter*.csv")))
    dfs = [x for x in (pd.read_csv(f) for f in files) if "wilcoxon" in x.columns]
    if not dfs: return print("skip permutation: no iteration files")
    d = pd.concat(dfs); d = d[d.design == "adjusted"]; tools = ["deseq2", "edger", "limmavoom", "wilcoxon"]
    real = {"deseq2": np.nan, "edger": np.nan, "limmavoom": np.nan, "wilcoxon": np.nan}
    try:
        n = pd.read_csv(R / "tcga/n_sig_summary.csv", index_col=0).loc["BRCA1_mut_adjusted"]
        real.update({t: float(n[t]) for t in ["deseq2", "edger", "limmavoom"]})
    except Exception as e: print("real BRCA1 (NB tools) not found:", e)
    try:
        w = pd.read_csv(R / "tcga/wilcoxon_check_summary.csv"); real["wilcoxon"] = float(w[(w.analysis == "BRCA1_mut") & (w.design == "adjusted")].mean_calls.iloc[0])
    except Exception as e: print("real BRCA1 Wilcoxon not found:", e)
    rng = np.random.default_rng(0); fig, ax = plt.subplots(1, 2, figsize=(11, 4.6))
    for i, t in enumerate(tools):
        ax[0].scatter(i + rng.uniform(-0.15, 0.15, len(d)), d[t], c=COL[t], alpha=0.8, s=30, zorder=3)
        ax[0].hlines(d[t].median(), i - 0.28, i + 0.28, colors="k", lw=2, zorder=4)
        if not np.isnan(real[t]): ax[0].scatter(i, real[t], marker="D", c="red", s=70, zorder=5, edgecolor="k", label="real BRCA1 mutants" if i == 0 else None)
        ax[1].scatter(i + rng.uniform(-0.15, 0.15, len(d)), d[f"p05_{t}"] * 100, c=COL[t], alpha=0.8, s=30)
    ax[0].set_yscale("symlog", linthresh=1); ax[0].set_ylim(-0.3, None); ax[0].set_yticks([0, 1, 10, 100, 1000]); ax[0].set_yticklabels(["0", "1", "10", "100", "1000"])
    ax[0].set_ylabel("genes called at q < 0.05 (log-ish scale)"); ax[0].set_title(f"Pseudo-mutant null: false calls per draw (n = {len(d)} draws)")
    ax[1].axhline(5, ls="--", c="k", lw=1); ax[1].set_ylabel("% of raw p-values < 0.05"); ax[1].set_title("Calibration (expected 5% under the null)")
    for a in ax: a.set_xticks(range(len(tools))); a.set_xticklabels([NICE[t] for t in tools])
    ax[0].legend(frameon=False, loc="upper right")
    cap(fig, "11 pseudo-mutants (7 Basal / 2 LumA / 2 LumB) drawn from wild-type tumors vs 150 controls, subtype-adjusted. Bars = medians. " "Real BRCA1 result used 300 controls, so the comparison is approximate.")
    plt.tight_layout(rect=(0, 0.07, 1, 1)); plt.savefig(FIG / "permutation_null.png", dpi=170); plt.close()


def fig_sva():
    bf = [pd.read_csv(f) for f in sorted(glob.glob(str(R / "tcga/perm/BRCA1_mut_iter*.csv")))]
    bf = [b for b in bf if "wilcoxon" in b.columns]
    sf = [pd.read_csv(f) for f in sorted(glob.glob(str(R / "tcga/sva/BRCA1_mut_null_iter*.csv")))]
    if not bf or not sf: return print("skip sva: need results/tcga/perm and results/tcga/sva files")
    base = pd.concat(bf); base = base[base.design == "adjusted"]
    j = base.merge(pd.concat(sf), on="iter", suffixes=("_base", "_sva"))
    rf = [pd.read_csv(f) for f in sorted(glob.glob(str(R / "tcga/sva/BRCA1_mut_real_iter*.csv")))]
    real = pd.concat(rf) if rf else None
    tools = ["deseq2", "edger", "limmavoom", "wilcoxon"]; dx = 0.18
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.6))
    for i, t in enumerate(tools):
        for a, col, scale in ((ax[0], t, 1), (ax[1], f"p05_{t}", 100)):
            for _, r in j.iterrows():
                a.plot([i - dx, i + dx], [scale * r[f"{col}_base"], scale * r[f"{col}_sva"]], c=COL[t], alpha=0.45, lw=1)
            a.scatter(np.full(len(j), i - dx), scale * j[f"{col}_base"], c=COL[t], s=24, zorder=3)
            a.scatter(np.full(len(j), i + dx), scale * j[f"{col}_sva"], c=COL[t], marker="s", s=24, zorder=3)
        if real is not None and real[t].notna().any():
            ax[0].scatter(np.full(real[t].notna().sum(), i + dx + 0.12), real[t].dropna(), marker="D", c="red", edgecolor="k", s=55, zorder=5,
                          label="real BRCA1 mutants, with SVs" if i == 0 else None)
    ax[0].set_yscale("symlog", linthresh=1); ax[0].set_ylim(-0.3, None); ax[0].set_yticks([0, 1, 10, 100, 1000]); ax[0].set_yticklabels(["0", "1", "10", "100", "1000"])
    ax[0].set_ylabel("genes called at q < 0.05 (log-ish scale)"); ax[0].set_title(f"False calls per draw (n = {len(j)} paired draws)")
    ax[1].axhline(5, ls="--", c="k", lw=1); ax[1].set_ylabel("% of raw p-values < 0.05"); ax[1].set_title("Calibration (expected 5%)")
    for a in ax: a.set_xticks(range(len(tools))); a.set_xticklabels([NICE[t] for t in tools])
    if real is not None: ax[0].legend(frameon=False, loc="upper right")
    nsv = int(pd.concat(sf).n_sv.max())
    cap(fig, f"Same pseudo-mutant draws without (circles, left) and with {nsv} surrogate variables (squares, right); each line is one draw. "
             "Subtype-adjusted in both. " + ("The estimated number of surrogate variables hit the cap of 10 in every draw."
             if (pd.concat(sf).n_sv == 10).all() else f"Surrogate variables per draw: {', '.join(map(str, pd.concat(sf).n_sv.astype(int)))}."))
    plt.tight_layout(rect=(0, 0.07, 1, 1)); plt.savefig(FIG / "sva_null.png", dpi=170); plt.close()


for fn in (fig_fdr_power, fig_unbalanced, fig_perm, fig_sva):
    try: fn()
    except Exception as e: print(f"{fn.__name__} failed: {e!r}")
print("figures in", FIG)
