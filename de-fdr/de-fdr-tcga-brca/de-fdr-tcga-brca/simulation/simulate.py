"""Simulate NB count data with known DE genes.

Mean: lognormal base means. Dispersion: DESeq2-style trend phi(mu)=a0+a1/mu with
lognormal gene-wise scatter. Library sizes vary by sample. 10% of genes are DE
with log2FC drawn from sign * (min_lfc + Exp). Writes counts/meta/truth per scenario.
"""
import argparse, pathlib
import numpy as np, pandas as pd


def simulate(n_per_group, n_genes=3000, frac_de=0.10, seed=0, a0=0.05, a1=1.0,
             disp_sd=0.5, min_lfc=0.6, lfc_scale=0.7):
    rng = np.random.default_rng(seed)
    base = np.clip(np.exp(rng.normal(4.0, 1.8, n_genes)), 0.5, 1e5)
    is_de = np.zeros(n_genes, bool)
    is_de[rng.choice(n_genes, int(frac_de * n_genes), replace=False)] = True
    lfc = np.zeros(n_genes)
    k = is_de.sum()
    lfc[is_de] = rng.choice([-1, 1], k) * (min_lfc + rng.exponential(lfc_scale, k))
    n = 2 * n_per_group
    cond = np.repeat([0, 1], n_per_group)
    sf = np.exp(rng.normal(0, 0.2, n))
    mu = base[:, None] * 2.0 ** (lfc[:, None] * cond[None, :]) * sf[None, :]
    phi = (a0 + a1 / base) * np.exp(rng.normal(0, disp_sd, n_genes))
    r = 1.0 / phi[:, None]
    counts = rng.negative_binomial(r, r / (r + mu))
    genes = [f"g{i:05d}" for i in range(n_genes)]
    samples = [f"s{j}" for j in range(n)]
    return (pd.DataFrame(counts, genes, samples),
            pd.DataFrame({"sample": samples, "condition": np.where(cond == 1, "B", "A")}),
            pd.DataFrame({"gene": genes, "is_de": is_de, "true_lfc": lfc, "base_mean": base, "phi": phi}))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results/sim")
    ap.add_argument("--n_per_group", nargs="+", type=int, default=[3, 5, 10])
    ap.add_argument("--reps", type=int, default=5)
    ap.add_argument("--n_genes", type=int, default=3000)
    a = ap.parse_args()
    for n in a.n_per_group:
        for rep in range(a.reps):
            d = pathlib.Path(a.out) / f"n{n}_rep{rep}"
            d.mkdir(parents=True, exist_ok=True)
            c, m, t = simulate(n, a.n_genes, seed=1000 * n + rep)
            c.to_csv(d / "counts.csv"); m.to_csv(d / "meta.csv", index=False); t.to_csv(d / "truth.csv", index=False)
    print("simulated", len(a.n_per_group) * a.reps, "datasets ->", a.out)
