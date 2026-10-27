"""Two-sample permutation test, implemented from scratch."""
import numpy as np


def permutation_test(x, y, n_perm=10000, stat="mean_diff", seed=0):
    """Two-sided permutation test for a difference between two groups.

    Uses the +1 correction, p = (1 + #{|T*| >= |T_obs|}) / (1 + n_perm), so a
    p-value can never be exactly zero. Returns (observed_stat, p_value).
    """
    rng = np.random.default_rng(seed)
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    pooled = np.concatenate([x, y])
    nx = x.size

    def _stat(a, b):
        if stat == "mean_diff":
            return a.mean() - b.mean()
        if stat == "median_diff":
            return np.median(a) - np.median(b)
        raise ValueError(stat)

    obs = _stat(x, y)
    count = 0
    for _ in range(n_perm):
        perm = rng.permutation(pooled)
        if abs(_stat(perm[:nx], perm[nx:])) >= abs(obs) - 1e-12:
            count += 1
    return obs, (1 + count) / (1 + n_perm)
