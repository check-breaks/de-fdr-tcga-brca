"""Run Python-side methods on every simulated dataset: t-test, Wilcoxon (on log-CPM), scratch NB GLM Wald."""
import sys, pathlib, time
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import numpy as np, pandas as pd
from scipy import stats
from src.scratch.nbglm import fit_nb_glm, wald_test

root = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else "results/sim")
for d in sorted(root.glob("n*_rep*")):
    if (d / "pvals_python.csv").exists(): continue
    if (d / "pvals_python.csv").exists(): continue
    c = pd.read_csv(d / "counts.csv", index_col=0)
    m = pd.read_csv(d / "meta.csv")
    cond = (m.condition == "B").to_numpy().astype(int)
    lib = c.sum(axis=0).to_numpy()
    sf = lib / np.exp(np.mean(np.log(lib)))          # simple depth size factors
    logcpm = np.log2(c.to_numpy() / lib * 1e6 + 1)
    keep = (c.to_numpy() >= 5).sum(1) >= 2            # basic low-count filter (same idea as filterByExpr)
    X = np.column_stack([np.ones(len(cond)), cond]); off = np.log(sf)
    t0 = time.time()
    rows = []
    for i, g in enumerate(c.index):
        if not keep[i]:
            rows.append((g, np.nan, np.nan, np.nan)); continue
        a, b = logcpm[i, cond == 0], logcpm[i, cond == 1]
        p_t = stats.ttest_ind(a, b, equal_var=False).pvalue
        p_w = stats.mannwhitneyu(a, b, method="auto").pvalue if len(a) > 0 else np.nan
        try:
            f = fit_nb_glm(c.iloc[i].to_numpy(), X, offset=off)
            p_nb = wald_test(f, 1)[1]
        except Exception:
            p_nb = np.nan
        rows.append((g, p_t, p_w, p_nb))
    out = pd.DataFrame(rows, columns=["gene", "p_ttest", "p_wilcoxon", "p_nbglm_scratch"])
    out.to_csv(d / "pvals_python.csv", index=False)
    print(d.name, f"{time.time()-t0:.1f}s", "tested", int(keep.sum()))
