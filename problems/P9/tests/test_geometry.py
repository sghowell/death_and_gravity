import numpy as np
from scipy.integrate import quad

from p9 import C_KM_S
from p9.geometry import dh_matrix, dm_matrix, grid, lcdm_u_nodes


def test_constant_u_gives_dm_equals_z():
    x = grid(50)
    z = np.array([0.013, 0.3, 0.51, 1.0, 2.33, 2.5])
    A = dm_matrix(x, z)
    assert np.allclose(A @ np.ones(51), z, rtol=1e-12, atol=1e-12)
    assert np.all(A >= -1e-15)


def test_piecewise_linear_integral_matches_quadrature():
    x = grid(20)
    rng = np.random.default_rng(0)
    u = 100 + 10 * rng.standard_normal(21)
    z = np.array([0.02, 0.4, 1.7])
    A = dm_matrix(x, z)

    def u_of_x(xx):
        return np.interp(xx, x, u)

    for zi, row in zip(z, A):
        exact = quad(lambda xx: u_of_x(xx) * np.exp(xx), 0, np.log1p(zi), points=list(x), limit=200)[0]
        assert abs(row @ u - exact) < 1e-9 * exact


def test_dh_matrix_interpolates_nodes():
    x = grid(10)
    A = dh_matrix(x, np.expm1(x))
    assert np.allclose(A, np.eye(11))


def test_lcdm_nodes_reproduce_distances_to_discretization_accuracy():
    x = grid(50)
    om, h_rd = 0.3, 100.0
    u = lcdm_u_nodes(x, om, h_rd)
    z = np.array([0.5, 1.0, 2.33])
    A = dm_matrix(x, z)
    exact = np.array([quad(lambda zz: 1 / np.sqrt(om * (1 + zz) ** 3 + 1 - om), 0, zi)[0] for zi in z])
    exact *= C_KM_S / (100 * h_rd)
    assert np.allclose(A @ u, exact, rtol=2e-4)
