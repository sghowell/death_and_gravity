"""Checks of the kappa-interpolation machinery (FORMULATION v2 §3.2)."""
import numpy as np
import pytest

from p9.geometry import (dm_matrix, geometric_grid, kappa_interp_slack, segment_index,
                         theta_bounds, kappa_second_derivative_bound)


def random_class_member(x, L, rng, u0=1.0):
    """Piecewise-linear u on grid x with |u_{k+1}-u_k| <= L h_k min(u_k,u_{k+1}), adversarial slopes."""
    u = np.empty(len(x)); u[0] = u0
    for k in range(len(x) - 1):
        h = x[k + 1] - x[k]
        s = rng.choice([-1.0, 1.0]) * rng.uniform(0.7, 1.0)
        # choose u_{k+1} = u_k (1 + s L h) if s>0 (then min = u_k), or u_k/(1 + |s| L h) if s<0
        u[k + 1] = u[k] * (1 + s * L * h) if s > 0 else u[k] / (1 + abs(s) * L * h)
    return u


def kappa_exact(x_nodes, u, xs):
    A = dm_matrix(x_nodes, np.expm1(xs))
    D = A @ u
    return np.log10(D / np.expm1(xs))


def test_geometric_grid_shape():
    x = geometric_grid()
    assert x[0] == 0.0 and abs(x[-1] - np.log1p(2.5)) < 1e-12
    assert np.all(np.diff(x) > 0)
    assert 60 < len(x) < 140


def test_constant_u_on_geometric_grid():
    x = geometric_grid()
    z = np.array([0.011, 0.05, 0.3, 1.0, 2.33])
    A = dm_matrix(x, z)
    assert np.allclose(A @ np.ones(len(x)), z, rtol=1e-12)


def test_theta_bounds_contain_truth():
    x = geometric_grid()
    rng = np.random.default_rng(1)
    L = 1.5
    for _ in range(50):
        u = random_class_member(x, L, rng)
        xs = rng.uniform(0.01, x[-1], 20)
        A = dm_matrix(x, np.expm1(xs))
        D = A @ u
        # u(xs) by interpolation
        k = segment_index(x, xs); t = (xs - x[k]) / (x[k + 1] - x[k])
        uxs = (1 - t) * u[k] + t * u[k + 1]
        theta = uxs * np.expm1(xs) / D
        for xx, th in zip(xs, theta):
            lo, hi = theta_bounds(xx, L)
            assert lo - 1e-12 <= th <= hi + 1e-12


@pytest.mark.parametrize("L", [0.5, 1.5, 5.0])
def test_interpolation_slack_is_rigorous(L):
    """Monte Carlo over adversarial class members: |kappa(x) - interp| <= e_k everywhere."""
    x = geometric_grid()
    e = kappa_interp_slack(x, L)
    assert np.all(np.isfinite(e)) and np.all(e > 0)
    rng = np.random.default_rng(2)
    worst = 0.0
    for _ in range(200):
        u = random_class_member(x, L, rng)
        kap_nodes = np.concatenate([[np.log10(u[0])], kappa_exact(x, u, x[1:])])
        xs = np.sort(rng.uniform(0.0, x[-1], 400))
        xs = xs[xs > 1e-6]
        kap = kappa_exact(x, u, xs)
        k = segment_index(x, xs); t = (xs - x[k]) / (x[k + 1] - x[k])
        interp = (1 - t) * kap_nodes[k] + t * kap_nodes[k + 1]
        ratio = np.abs(kap - interp) / e[k]
        worst = max(worst, ratio.max())
    assert worst <= 1.0, f"slack bound violated: worst ratio {worst}"
    # and the bound should not be absurdly loose everywhere
    assert worst > 0.05


def test_slack_small_where_it_matters():
    x = geometric_grid()
    e = kappa_interp_slack(x, 1.5)
    # 5*e in magnitudes should be well below SN errors (0.1 mag) everywhere
    assert np.max(5 * e) < 0.02
