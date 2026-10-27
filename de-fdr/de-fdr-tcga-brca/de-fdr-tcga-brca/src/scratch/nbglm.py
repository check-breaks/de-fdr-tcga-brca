"""Negative binomial GLM (NB2, log link, offset) fit by IRLS, from scratch.

Var(Y) = mu + alpha * mu^2.  Dispersion alpha is estimated by maximum
likelihood given the mean model (alternating with IRLS). This is a simplified
analogue of what DESeq2/edgeR do: it has NO dispersion shrinkage, which is
exactly why it should be anti-conservative at small n (a nice thing to show).
"""
import numpy as np
from scipy import stats, optimize
from scipy.special import gammaln


def nb_loglik(y, mu, alpha):
    r = 1.0 / alpha
    return np.sum(
        gammaln(y + r) - gammaln(r) - gammaln(y + 1)
        + r * np.log(r / (r + mu)) + y * np.log(mu / (r + mu))
    )


def _irls(y, X, offset, alpha, beta=None, max_iter=50, tol=1e-8, ridge=1e-6):
    n, p = X.shape
    if beta is None:
        beta = np.zeros(p)
        beta[0] = np.log(max(y.mean(), 0.1)) - offset.mean()
    for _ in range(max_iter):
        eta = np.clip(X @ beta + offset, -30, 30)
        mu = np.exp(eta)
        w = mu / (1.0 + alpha * mu)
        z = (eta - offset) + (y - mu) / mu
        XtW = X.T * w
        A = XtW @ X + ridge * np.eye(p)
        new = np.linalg.solve(A, XtW @ z)
        if np.max(np.abs(new - beta)) < tol:
            beta = new
            break
        beta = new
    eta = np.clip(X @ beta + offset, -30, 30)
    mu = np.exp(eta)
    w = mu / (1.0 + alpha * mu)
    cov = np.linalg.inv((X.T * w) @ X + ridge * np.eye(p))
    return beta, mu, cov


def fit_nb_glm(y, X, offset=None, alpha=None, n_outer=4):
    """Fit one gene. If alpha is None, estimate it by ML. Returns dict."""
    y = np.asarray(y, float)
    X = np.asarray(X, float)
    offset = np.zeros_like(y) if offset is None else np.asarray(offset, float)
    est_alpha = alpha is None
    a = 0.1 if est_alpha else alpha
    beta = None
    for _ in range(n_outer if est_alpha else 1):
        beta, mu, cov = _irls(y, X, offset, a, beta)
        if not est_alpha:
            break
        res = optimize.minimize_scalar(
            lambda la: -nb_loglik(y, mu, np.exp(la)),
            bounds=(np.log(1e-6), np.log(50.0)), method="bounded",
        )
        a = float(np.exp(res.x))
    beta, mu, cov = _irls(y, X, offset, a, beta)
    return {"beta": beta, "se": np.sqrt(np.diag(cov)), "alpha": a,
            "mu": mu, "cov": cov, "loglik": nb_loglik(y, mu, a)}


def wald_test(fit, coef):
    z = fit["beta"][coef] / fit["se"][coef]
    return z, 2 * stats.norm.sf(abs(z))


def lrt(y, X_full, X_red, offset=None):
    """Likelihood-ratio test of the columns dropped in X_red.

    Dispersion is estimated once from the FULL model and held fixed for both fits
    (as DESeq2/edgeR do). Re-estimating it under the reduced model lets alpha
    absorb the effect being tested and makes the LRT conservative.
    """
    f = fit_nb_glm(y, X_full, offset)
    r = fit_nb_glm(y, X_red, offset, alpha=f["alpha"])
    stat = max(2 * (f["loglik"] - r["loglik"]), 0.0)
    df = X_full.shape[1] - X_red.shape[1]
    return stat, stats.chi2.sf(stat, df)
