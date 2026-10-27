"""Realized FDR and power at nominal BH-FDR 5%, per method and sample size; writes figures + table."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from src.scratch.multiple_testing import benjamini_hochberg

root = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else "results/sim")
alpha = 0.05
rows = []
for d in sorted(root.glob("n*_rep*")):
    n = int(d.name.split("_")[0][1:])
    truth = pd.read_csv(d / "truth.csv").set_index("gene")
    pv = truth[[]].copy()
    for f in ("pvals_python.csv", "pvals_r.csv"):
        if (d / f).exists():
            pv = pv.join(pd.read_csv(d / f).set_index("gene"))
    for meth in pv.columns:
        q = benjamini_hochberg(pv[meth].to_numpy())
        called = np.nan_to_num(q, nan=1.0) <= alpha
        is_de = truth["is_de"].to_numpy()
        tp, fp = (called & is_de).sum(), (called & ~is_de).sum()
        rows.append(dict(method=meth.replace("p_", ""), n_per_group=n, rep=d.name,
                         fdr=fp / max(called.sum(), 1), power=tp / is_de.sum(), n_called=int(called.sum())))
res = pd.DataFrame(rows)
res.to_csv(root.parent / "simulation_results_raw.csv", index=False)
summ = res.groupby(["method", "n_per_group"]).agg(fdr=("fdr", "mean"), fdr_sd=("fdr", "std"),
                                                  power=("power", "mean"), power_sd=("power", "std")).reset_index()
summ.to_csv(root.parent / "simulation_summary.csv", index=False)
print(summ.round(3).to_string(index=False))

fig, ax = plt.subplots(1, 2, figsize=(10, 4))
for meth, g in summ.groupby("method"):
    ax[0].errorbar(g.n_per_group, g.fdr, g.fdr_sd, marker="o", capsize=3, label=meth)
    ax[1].errorbar(g.n_per_group, g.power, g.power_sd, marker="o", capsize=3, label=meth)
ax[0].axhline(alpha, ls="--", c="k", lw=1); ax[0].set_title("Realized FDR (nominal 5%, BH)")
ax[1].set_title("Power (true positive rate)")
for a in ax: a.set_xlabel("samples per group"); a.set_xticks(sorted(summ.n_per_group.unique()))
ax[0].set_ylabel("FDR"); ax[1].set_ylabel("power"); ax[1].legend(fontsize=8)
plt.tight_layout(); plt.savefig(root.parent / "figures" / "fdr_power.png", dpi=150)
