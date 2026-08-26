import numpy as np

from p9.geometry import dm_matrix, geometric_grid, kappa_difference_bounds
from test_kappa import random_class_member, kappa_exact


def test_kappa_difference_bounds_are_rigorous_and_tight():
    x = geometric_grid()
    rng = np.random.default_rng(3)
    for L in [0.5, 1.5, 5.0]:
        lo, hi = kappa_difference_bounds(x, L)
        assert np.all(hi > lo)
        worst = -np.inf
        for _ in range(200):
            u = random_class_member(x, L, rng)
            kap = np.concatenate([[np.log10(u[0])], kappa_exact(x, u, x[1:])])
            d = np.diff(kap)
            assert np.all(d >= lo - 1e-12) and np.all(d <= hi + 1e-12)
            worst = max(worst, np.max((d - lo) / (hi - lo)))
        assert worst > 0.5  # adversarial members approach the bounds
        # at low z the width (in mag) must be small
        w = 5 * (hi - lo)
        assert np.max(w[x[:-1] < 0.1]) < 0.03 * L, w[x[:-1] < 0.1].max()
