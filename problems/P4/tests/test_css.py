"""Stage A tests: sonic-point data, analytic branches, constraint identity."""
import numpy as np
import pytest
from scipy.integrate import solve_ivp

from p4 import css, taylor
from p4.tps import TPS

V0_LIST = [0.05, 0.112439401388092, 0.3, 0.5, -0.7, -0.9]


@pytest.mark.parametrize("V0", V0_LIST)
def test_sonic_data_vanishing_determinant(V0):
    y0 = css.sonic_data(V0)
    q = css.coeffs(*y0)
    Delta, P, Q = css.det_and_numerators(q)
    scale = abs(q["a"] * q["d"]) + abs(q["b"] * q["c"])
    assert abs(Delta) < 1e-12 * scale                       # KHA99 (214)
    assert abs(q["a"] * q["f"] - q["c"] * q["e"]) < 1e-12 * scale  # KHA99 (215)
    assert abs(P) < 1e-12 * scale and abs(Q) < 1e-12 * scale
    assert abs(css.constraint(y0)) < 1e-13                  # KHA99 (211)
    # relative velocity to the x = const line is the sound speed
    assert abs(css.v_rel(y0[1], V0) - 1 / np.sqrt(3)) < 1e-14


@pytest.mark.parametrize("V0", V0_LIST)
def test_first_order_branches_distinct_and_consistent(V0):
    y0, branches = css.sonic_branches(V0)
    assert len(branches) == 2
    assert abs(branches[0][3] - branches[1][3]) > 1e-3       # distinct V1
    q0 = css.coeffs(*y0)
    ell = np.array([q0["c"], -q0["a"]])
    for y1 in branches:
        ser = [TPS(np.array([y0[i], y1[i], 0.0])) for i in range(4)]
        r3, r4 = css.fluid_residuals(*ser, ser[2].deriv(), ser[3].deriv())
        assert abs(r3.c[0]) < 1e-12 and abs(r4.c[0]) < 1e-12          # O(x^0)
        assert abs(ell[0] * r3.c[1] + ell[1] * r4.c[1]) < 1e-10        # ell . O(x^1)
        q = css.coeffs(*ser)
        assert abs((ser[0] * q["FA"]).c[0] - y1[0]) < 1e-12
        assert abs((ser[1] * q["FN"]).c[0] - y1[1]) < 1e-12


def test_discriminant_sign_and_friedmann_point():
    assert css.first_order_discriminant(0.3) > 0
    assert css.first_order_discriminant(-0.3) < 0
    assert css.first_order_discriminant(-0.8) > 0
    assert abs(css.first_order_discriminant(-1 / np.sqrt(3))) < 1e-12
    # saddle on (0, 1/sqrt3): opposite-sign eigenvalues of the desingularised flow
    mu = css.sonic_eigenvalues(0.3)
    assert np.real(mu[0]) * np.real(mu[1]) < 0


def test_constraint_identity_along_solutions():
    """The momentum constraint (KHA99 211) is preserved by the 4D CSS flow."""
    for V0, branch in ((0.3, 0), (0.3, 1), (0.112439401388092, 1)):
        coef = taylor.background_series(V0, branch, K=30)
        for sign in (-1, 1):
            y1 = taylor.eval_series(coef, sign * 0.02)
            sol = solve_ivp(css.rhs_plain, (sign * 0.02, sign * 1.2), y1, method="DOP853",
                            rtol=1e-12, atol=1e-14, dense_output=True)
            assert sol.status == 0
            for x in np.linspace(sign * 0.02, sign * 1.2, 8):
                y = sol.sol(x)
                q = css.coeffs(*y)
                # invariant up to integration error (the constraint surface is repelling
                # for the backward-x flow, so errors are amplified inward)
                assert abs(q["FA"] - q["G"]) < 1e-8 * (1 + abs(q["FA"]) + abs(q["G"]))


def test_series_matches_ode():
    coef = taylor.background_series(0.112439401388092, 1, K=40)
    y1 = taylor.eval_series(coef, -0.005)
    sol = solve_ivp(css.rhs_plain, (-0.005, -0.08), y1, method="DOP853", rtol=1e-13, atol=1e-15)
    assert np.max(np.abs(sol.y[:, -1] - taylor.eval_series(coef, -0.08))) < 1e-12


def test_scaled_and_reduced_rhs_agree():
    x = -0.7
    y = np.array([1.9, 3.1, 0.25, -0.1])
    yh = css.to_scaled(x, y)
    d_plain = css.rhs_plain(x, y)
    d_scaled = css.rhs_scaled(x, yh)
    # d/dx of scaled variables from the plain derivatives
    expect = np.array([d_plain[0] * np.exp(-2 * x) - 2 * yh[0], d_plain[1] * np.exp(x) + yh[1],
                       d_plain[2] * np.exp(-2 * x) - 2 * yh[2], d_plain[3] * np.exp(-x) - yh[3]])
    assert np.allclose(d_scaled, expect, rtol=1e-12, atol=1e-12)
    # reduced system on the constraint surface
    N, W, V = 3.1, 0.25, -0.1
    A = 1 + css.Am1_constraint(N, W, V)
    assert abs(css.constraint([A, N, W, V])) < 1e-14
    assert np.allclose(css.rhs_plain3(x, [N, W, V]), css.rhs_plain(x, [A, N, W, V])[1:], rtol=1e-12)
