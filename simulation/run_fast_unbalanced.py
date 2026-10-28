"""Quick Python-only baselines on the unbalanced sims: Welch t-test and Wilcoxon rank-sum (naive) and OLS t-test /
subtype-residual Wilcoxon (adjusted), all on log-CPM. Mostly a sanity check of the simulation / evaluation code."""
import sys, pathlib
import numpy as np, pandas as pd
from scipy import stats

root = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else "results/sim_unbal")
for d in sorted(root.glob("*_n*_rep*")):
    c = pd.read_csv(d / "counts.csv", index_col=0); m = pd.read_csv(d / "meta.csv")
    cond = (m.condition == "case").to_numpy().astype(float)
    # median-of-ratios size factors (library-size CPM alone is biased when many genes shift between subtypes)
    x = c.to_numpy(float); okg = (x > 0).all(1); lg = np.log(x[okg])
    sf = np.exp(np.median(lg - lg.mean(1, keepdims=True), axis=0))
    depth = sf * c.sum(axis=0).to_numpy(float).mean()
    y = np.log2(x / depth * 1e6 + 1)
    keep = (c.to_numpy() >= 5).sum(1) >= 5
    X = np.column_stack([np.ones(len(cond)), cond, pd.get_dummies(m.subtype, drop_first=True).to_numpy(float)])
    XtXi = np.linalg.inv(X.T @ X); B = XtXi @ X.T @ y.T
    res = y.T - X @ B; df = X.shape[0] - X.shape[1]
    se = np.sqrt((res ** 2).sum(0) / df * XtXi[1, 1])
    p_adj = 2 * stats.t.sf(np.abs(B[1] / se), df)
    p_naive = stats.ttest_ind(y[:, cond == 1], y[:, cond == 0], axis=1, equal_var=False).pvalue
    # Wilcoxon rank-sum (asymptotic): naive on log-CPM; adjusted on residuals after regressing out subtype
    ca, ct = cond == 1, cond == 0
    p_wnaive = stats.mannwhitneyu(y[:, ca], y[:, ct], axis=1, method="asymptotic").pvalue
    X0 = np.column_stack([np.ones(len(cond)), pd.get_dummies(m.subtype, drop_first=True).to_numpy(float)])
    R = (y.T - X0 @ (np.linalg.inv(X0.T @ X0) @ X0.T @ y.T)).T
    p_wadj = stats.mannwhitneyu(R[:, ca], R[:, ct], axis=1, method="asymptotic").pvalue
    out = pd.DataFrame({"gene": c.index, "p_ttest_naive": np.where(keep, p_naive, np.nan),
                        "p_ttest_adjusted": np.where(keep, p_adj, np.nan),
                        "p_wilcoxon_naive": np.where(keep, p_wnaive, np.nan),
                        "p_wilcoxon_adjusted": np.where(keep, p_wadj, np.nan)})
    out.to_csv(d / "pvals_python.csv", index=False); print(d.name)
