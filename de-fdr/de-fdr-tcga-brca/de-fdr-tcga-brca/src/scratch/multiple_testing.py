"""Multiple-testing corrections implemented from scratch."""
import numpy as np


def benjamini_hochberg(pvals):
    """BH-adjusted p-values (q-values). NaNs are ignored and returned as NaN.

    Sort p ascending, compute p_(i) * m / i, then enforce monotonicity by a
    reverse cumulative minimum, and cap at 1.
    """
    p = np.asarray(pvals, dtype=float)
    out = np.full(p.shape, np.nan)
    ok = ~np.isnan(p)
    pv = p[ok]
    m = pv.size
    if m == 0:
        return out
    order = np.argsort(pv)
    ranked = pv[order] * m / np.arange(1, m + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adj = np.empty(m)
    adj[order] = np.minimum(ranked, 1.0)
    out[ok] = adj
    return out


def bonferroni(pvals):
    p = np.asarray(pvals, dtype=float)
    return np.minimum(p * np.sum(~np.isnan(p)), 1.0)
