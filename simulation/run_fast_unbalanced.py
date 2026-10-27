"""Quick Python-only baselines on the unbalanced sims: Welch t-test (naive) and OLS t-test with subtype
covariate (adjusted), both on log-CPM. Mostly a sanity check of the simulation / evaluation code."""
import sys, pathlib
import numpy as np, pandas as pd
from scipy import stats

root = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else "results/sim_unbal")
for d in sorted(root.glob("*_n*_rep*")):
    c = pd.read_csv(d / "counts.csv", index_col=0); m = pd.read_csv(d / "meta.csv")
    cond = (m.condition == "case").to_numpy().astype(float)
    lib = c.sum(axis=0).to_numpy(); y = np.log2(c.to_numpy() / lib * 1e6 + 1)
    keep = (c.to_numpy() >= 5).sum(1) >= 5
    X = np.column_stack([np.ones(len(cond)), cond, pd.get_dummies(m.subtype, drop_first=True).to_numpy(float)])
    XtXi = np.linalg.inv(X.T @ X); B = XtXi @ X.T @ y.T
    res = y.T - X @ B; df = X.shape[0] - X.shape[1]
    se = np.sqrt((res ** 2).sum(0) / df * XtXi[1, 1])
    p_adj = 2 * stats.t.sf(np.abs(B[1] / se), df)
    p_naive = stats.ttest_ind(y[:, cond == 1], y[:, cond == 0], axis=1, equal_var=False).pvalue
    out = pd.DataFrame({"gene": c.index, "p_ttest_naive": np.where(keep, p_naive, np.nan),
                        "p_ttest_adjusted": np.where(keep, p_adj, np.nan)})
    out.to_csv(d / "pvals_python.csv", index=False); print(d.name)
