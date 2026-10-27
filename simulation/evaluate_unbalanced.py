"""FDR / power at nominal BH 5% for every pvals_*.csv column, by scenario, n_case and method."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from src.scratch.multiple_testing import benjamini_hochberg

root = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else "results/sim_unbal"); alpha = 0.05
rows = []
for d in sorted(root.glob("*_n*_rep*")):
    sc, nc = d.name.split("_")[0], int(d.name.split("_")[1][1:])
    truth = pd.read_csv(d / "truth.csv").set_index("gene"); pv = truth[[]].copy()
    for f in d.glob("pvals_*.csv"):
        pv = pv.join(pd.read_csv(f).set_index("gene"))
    is_de = truth["is_de"].to_numpy()
    for col in pv.columns:
        called = np.nan_to_num(benjamini_hochberg(pv[col].to_numpy()), nan=1.0) <= alpha
        fp = (called & ~is_de).sum()
        rows.append(dict(scenario=sc, n_case=nc, method=col[2:], fdr=fp / max(called.sum(), 1),
                         power=(called & is_de).sum() / is_de.sum(), n_called=int(called.sum()), n_false=int(fp)))
res = pd.DataFrame(rows)
s = res.groupby(["scenario", "n_case", "method"]).agg(reps=("fdr", "size"), fdr=("fdr", "mean"), fdr_sd=("fdr", "std"),
        power=("power", "mean"), n_called=("n_called", "mean"), n_false=("n_false", "mean")).reset_index()
s.to_csv(root.parent / "simulation_unbalanced_summary.csv", index=False)
print(s.round(3).to_string(index=False))
fig, ax = plt.subplots(1, 2, figsize=(11, 4), sharey=True)
for a, sc in zip(ax, ["mixed", "confounded"]):
    g = s[s.scenario == sc].groupby("method").fdr.mean().sort_values()
    a.barh(g.index, g.values); a.axvline(alpha, ls="--", c="k"); a.set_title(f"{sc}: realized FDR (nominal 5%)")
plt.tight_layout(); (root.parent / "figures").mkdir(exist_ok=True)
plt.savefig(root.parent / "figures" / "fdr_unbalanced.png", dpi=150)
