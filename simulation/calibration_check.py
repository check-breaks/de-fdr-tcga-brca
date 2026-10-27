"""Null calibration of the scratch NB GLM: true alpha vs ML-estimated alpha.
Shows false-positive rate is ~5% with the true dispersion, and inflated with gene-wise ML dispersion
at small n (the motivation for dispersion shrinkage in DESeq2/edgeR)."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import numpy as np
from src.scratch.nbglm import fit_nb_glm, wald_test

rng = np.random.default_rng(0); mu, alpha, G = 100, 0.1, 1500; r = 1 / alpha
for n in (3, 5, 10, 20):
    X = np.column_stack([np.ones(2 * n), np.repeat([0, 1], n)]); pt, pe = [], []
    for _ in range(G):
        y = rng.negative_binomial(r, r / (r + mu), size=2 * n)
        pt.append(wald_test(fit_nb_glm(y, X, alpha=alpha), 1)[1]); pe.append(wald_test(fit_nb_glm(y, X), 1)[1])
    print(f"n/group={n:2d}  FPR@0.05  true alpha: {np.mean(np.array(pt)<.05):.3f}   estimated alpha: {np.mean(np.array(pe)<.05):.3f}")
