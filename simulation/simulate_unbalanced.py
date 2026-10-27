"""Unbalanced case/control simulation with a latent 'subtype' structure, mimicking the TCGA-BRCA design.

Controls are a mixture of 3 latent groups (like PAM50 subtypes). Groups 1 and 2 differ from group 0 in many
genes (subtype effect, NOT related to case status). Cases carry a true effect in 10% of genes.
Scenarios:
  mixed       cases drawn from the same group mix as controls            -> subtype is NOT a confounder
  confounded  cases enriched in group 1 (like 7/11 BRCA1 mutants basal)  -> subtype IS a confounder
Truth = genes with a case effect; group-effect genes are NOT DE for case status, so a naive design that
picks them up produces false positives.
"""
import argparse, pathlib
import numpy as np, pandas as pd

CTRL_MIX = np.array([0.55, 0.20, 0.25])
CASE_MIX = {"mixed": CTRL_MIX, "confounded": np.array([0.20, 0.64, 0.16])}


def simulate(scenario, n_case, n_ctrl, n_genes=3000, frac_de=0.10, seed=0, a0=0.05, a1=1.0,
             disp_sd=0.5, min_lfc=0.6, lfc_scale=0.7, frac_group_genes=0.30, group_sd=1.0):
    rng = np.random.default_rng(seed)
    base = np.clip(np.exp(rng.normal(4.0, 1.8, n_genes)), 0.5, 1e5)
    is_de = np.zeros(n_genes, bool)
    is_de[rng.choice(n_genes, int(frac_de * n_genes), replace=False)] = True
    lfc = np.zeros(n_genes); k = is_de.sum()
    lfc[is_de] = rng.choice([-1, 1], k) * (min_lfc + rng.exponential(lfc_scale, k))
    # subtype effects (groups 1, 2 vs group 0), on genes chosen independently of is_de
    gl = np.zeros((n_genes, 3))
    for g in (1, 2):
        aff = rng.random(n_genes) < frac_group_genes
        gl[aff, g] = rng.normal(0, group_sd, aff.sum())
    n = n_case + n_ctrl
    cond = np.r_[np.ones(n_case, int), np.zeros(n_ctrl, int)]
    grp = np.r_[rng.choice(3, n_case, p=CASE_MIX[scenario]), rng.choice(3, n_ctrl, p=CTRL_MIX)]
    sf = np.exp(rng.normal(0, 0.2, n))
    log2mu = np.log2(base)[:, None] + gl[:, grp] + lfc[:, None] * cond[None, :]
    mu = 2.0 ** log2mu * sf[None, :]
    phi = (a0 + a1 / base) * np.exp(rng.normal(0, disp_sd, n_genes))
    r = 1.0 / phi[:, None]
    counts = rng.negative_binomial(r, r / (r + mu))
    genes = [f"g{i:05d}" for i in range(n_genes)]; samples = [f"s{j}" for j in range(n)]
    return (pd.DataFrame(counts, genes, samples),
            pd.DataFrame({"sample": samples, "condition": np.where(cond == 1, "case", "ctrl"),
                          "subtype": [f"G{x}" for x in grp]}),
            pd.DataFrame({"gene": genes, "is_de": is_de, "true_lfc": lfc, "base_mean": base}))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results/sim_unbal")
    ap.add_argument("--n_case", nargs="+", type=int, default=[7, 11])
    ap.add_argument("--n_ctrl", type=int, default=300)
    ap.add_argument("--reps", type=int, default=10)
    ap.add_argument("--n_genes", type=int, default=3000)
    a = ap.parse_args()
    for sc in CASE_MIX:
        for nc in a.n_case:
            for rep in range(a.reps):
                d = pathlib.Path(a.out) / f"{sc}_n{nc}_rep{rep}"; d.mkdir(parents=True, exist_ok=True)
                c, m, t = simulate(sc, nc, a.n_ctrl, a.n_genes, seed=7919 * nc + 31 * rep + (sc == "confounded"))
                c.to_csv(d / "counts.csv"); m.to_csv(d / "meta.csv", index=False); t.to_csv(d / "truth.csv", index=False)
    print("simulated ->", a.out)
