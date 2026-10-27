"""Validation of from-scratch implementations. Run: pytest -q  (or python tests/test_scratch.py)"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import numpy as np
from scipy import stats, optimize
from src.scratch.multiple_testing import benjamini_hochberg
from src.scratch.permutation import permutation_test
from src.scratch.nbglm import fit_nb_glm, wald_test, lrt, nb_loglik


def test_bh_matches_scipy():
    rng = np.random.default_rng(1)
    p = np.concatenate([rng.uniform(size=900), rng.beta(0.1, 5, size=100)])
    np.testing.assert_allclose(benjamini_hochberg(p), stats.false_discovery_control(p, method="bh"))


def test_bh_known_example():
    p = np.array([0.01, 0.04, 0.03, 0.005])
    np.testing.assert_allclose(benjamini_hochberg(p), [0.02, 0.04, 0.04, 0.02])


def test_bh_handles_nan():
    out = benjamini_hochberg([0.01, np.nan, 0.5])
    assert np.isnan(out[1]) and not np.isnan(out[0])


def test_permutation_close_to_welch_on_normal_data():
    rng = np.random.default_rng(2)
    x, y = rng.normal(0, 1, 30), rng.normal(0.8, 1, 30)
    _, p_perm = permutation_test(x, y, n_perm=20000, seed=3)
    p_t = stats.ttest_ind(x, y, equal_var=False).pvalue
    assert abs(np.log10(p_perm) - np.log10(p_t)) < 0.35


def test_permutation_null_is_uniform():
    rng = np.random.default_rng(4)
    ps = [permutation_test(rng.normal(size=8), rng.normal(size=8), n_perm=499, seed=i)[1] for i in range(300)]
    assert stats.kstest(ps, "uniform").pvalue > 0.01


def _design(n):
    g = np.repeat([0, 1], n // 2)
    return np.column_stack([np.ones(n), g]), g


def test_nbglm_group_means_closed_form():
    # For a 2-group factor design with log link and no offset, the MLE of beta
    # does not depend on alpha: exp(b0)=mean(group0), exp(b0+b1)=mean(group1).
    rng = np.random.default_rng(5)
    X, g = _design(12)
    y = rng.negative_binomial(5, 5 / (5 + np.where(g == 1, 80, 40)))
    fit = fit_nb_glm(y, X)
    np.testing.assert_allclose(np.exp(fit["beta"][0]), y[g == 0].mean(), rtol=1e-4)
    np.testing.assert_allclose(np.exp(fit["beta"].sum()), y[g == 1].mean(), rtol=1e-4)


def test_nbglm_matches_bruteforce_mle_with_offset_and_covariate():
    rng = np.random.default_rng(6)
    n = 40
    X = np.column_stack([np.ones(n), rng.integers(0, 2, n), rng.normal(size=n)])
    off = rng.normal(0, 0.3, n)
    mu = np.exp(X @ [3.0, 0.5, -0.3] + off)
    alpha = 0.2
    y = rng.negative_binomial(1 / alpha, (1 / alpha) / (1 / alpha + mu))
    fit = fit_nb_glm(y, X, offset=off, alpha=alpha)
    nll = lambda b: -nb_loglik(y, np.exp(X @ b + off), alpha)
    ref = optimize.minimize(nll, np.zeros(3), method="BFGS", options={"gtol": 1e-8}).x
    np.testing.assert_allclose(fit["beta"], ref, atol=1e-3)


def test_nbglm_alpha_recovery():
    rng = np.random.default_rng(7)
    n, alpha = 400, 0.3
    X = np.ones((n, 1))
    mu = 50.0
    y = rng.negative_binomial(1 / alpha, (1 / alpha) / (1 / alpha + mu), size=n)
    assert abs(fit_nb_glm(y, X)["alpha"] - alpha) < 0.06


def test_wald_and_lrt_agree_at_moderate_n():
    rng = np.random.default_rng(8)
    X, g = _design(30)
    y = rng.negative_binomial(8, 8 / (8 + np.where(g == 1, 120, 60)))
    fit = fit_nb_glm(y, X)
    _, p_w = wald_test(fit, 1)
    _, p_l = lrt(y, X, X[:, :1])
    assert abs(np.log10(p_w) - np.log10(p_l)) < 0.5


def test_against_statsmodels_if_available():
    try:
        import statsmodels.api as sm
    except ImportError:
        return  # skipped (statsmodels not installed here)
    rng = np.random.default_rng(9)
    X, g = _design(20)
    y = rng.negative_binomial(6, 6 / (6 + np.where(g == 1, 90, 50)))
    mine = fit_nb_glm(y, X, alpha=0.25)
    ref = sm.GLM(y, X, family=sm.families.NegativeBinomial(alpha=0.25)).fit()
    np.testing.assert_allclose(mine["beta"], ref.params, atol=1e-4)
    np.testing.assert_allclose(mine["se"], ref.bse, rtol=1e-3)


if __name__ == "__main__":
    fails = 0
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            try:
                fn(); print("PASS", name)
            except Exception as e:
                fails += 1; print("FAIL", name, "->", repr(e)[:200])
    sys.exit(fails)
