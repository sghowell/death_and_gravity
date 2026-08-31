"""S4-1 (Theorem C prerequisites): exact 2nd-order reduction of the linearised system.

Exact identity tests for p4.validated.tc_reduce (fmpq_mpoly arithmetic; no floats), a
float consistency check of the reduced 2x2 quotient system against the full linearised
system, and the known-answer reproduction of kappa_1 and kbar by the reduced problem.
"""
import numpy as np
import pytest
from p4.validated.tc_reduce import Reduction


@pytest.fixture(scope="module")
def red():
    return Reduction()


@pytest.fixture(scope="module")
def rp(red):
    return red.reduce_pair()


def test_gauge_identity_exact(red):
    """(S Delta~)^3 [P g' - (DQ - Psi - kappa P_s) g] == 0 as polynomials, for
    g = (A', N' + kappa N, W', V'): the gauge vector solves the 4D linearised system
    identically -- no background-constraint elimination is even needed (closes the
    'ball identity to order 38' caveat of s2-theorem-b.md section 3.7/3.8)."""
    assert all(r.is_zero() for r in red.gauge_residual())


def test_gauge_on_constraint_surface(red):
    """c(g) == 0 modulo the background constraint A = Anum/S."""
    assert red.gauge_constraint_residual().is_zero()


def test_pair_annihilates_gauge(red):
    """lchi . ghat == 0 and leta . ghat == 0 identically (the quotient functionals
    chi, eta vanish on the gauge solution)."""
    g = red.gauge_vector()
    z = red.ctx.constant(0)
    assert sum((red.lchi[i] * g[i] for i in range(4)), z).is_zero()
    assert sum((red.leta[i] * g[i] for i in range(4)), z).is_zero()


def test_closure_identities_exact(red, rp):
    """dn e L_chi = Achi lchi + Bchi leta + Cchi c  (and the eta analogue) modulo the
    background constraint: the (chi, eta) pair closes on the constraint surface --
    THE exact 2nd-order reduction (components Ap, Wp, Vp match by construction; the
    Np component is the nontrivial residual, eliminated exactly)."""
    assert red.elim(rp["chi"][3]).is_zero()
    assert red.elim(rp["eta"][3]).is_zero()


def test_weight_identities(red):
    """Delta~ = 4 S W D and D = 3(1+NV)^2 - (N+V)^2, i.e. the sonic weight
    w = -D = 3 (N+V)^2 (1/3 - v_rel^2) is GHJS's w times the positive 3(N+V)^2."""
    N, V, W, S = red.N, red.V, red.W, red.S
    assert (red.Delta - 4 * S * W * red.Dpoly).is_zero()
    assert (red.Dpoly - (3 * (1 + N * V) ** 2 - (N + V) ** 2)).is_zero()


def test_D2_factorisation(red, rp):
    """D2 = S^2 Delta~^2 dn e with dn = (k - A) S, e = Delta~ (Q_1 + k N): the reduced
    system's only kappa-dependent denominators are (A - kappa) and (F_N + kappa)."""
    assert (rp["D2"] - red.SD * red.SD * red.dn * red.e).is_zero()


def test_reduced_matches_full_system():
    """Float: the reduced 2x2 integration reproduces (chi, eta) of the full linearised
    system at x = -8 (kappa = 2.5, real, no detour)."""
    from p4.validated.tc_reduced_eigs import direct_final, reduced_final
    cr, er = reduced_final(2.5)
    cd, ed = direct_final(2.5)
    assert abs(cr - cd) / abs(cd) < 1e-8
    assert abs(er - ed) / abs(ed) < 1e-4


def test_reduced_eigenvalues_known_answer():
    """Float known-answer test (loose tolerance): the reduced scalar problem's zeros
    reproduce kappa_1 and the gauge zero kbar to <= 1e-9 (certified values from
    Theorem B, s2-theorem-b.md section 3.6/3.7)."""
    from p4 import taylor
    from p4.validated.tc_reduced_eigs import V0, find_zero
    bg = taylor.background_series(V0, 1, K=36)
    z1 = find_zero(2.81055, 2.81056, x_end=-9.0, bg=bg)
    assert abs(z1 - 2.8105525488271472) < 1e-9
    assert abs(z1.imag) < 1e-9
    zg = find_zero(0.3556, 0.3558, x_end=-9.0, bg=bg)
    assert abs(zg - 0.355699203710964) < 1e-9


def test_survey_sign_quantities():
    """Float: the route-(b) sign facts at three interior points -- weight w = -D > 0,
    w'/w < 0, det C < 0 (so a0's kappa^2-density has one sign), and the exact
    identities -D = 3(N+V)^2 w_GHJS and D detC = 3 - V^2 (detC = (3-V^2)/D < 0)."""
    from p4 import css, taylor
    from p4.validated.tc_reduced_eigs import V0
    bg = taylor.background_series(V0, 1, K=36)
    for x in (-0.05, -0.5, -2.0):
        A, N, W, V = taylor.eval_series(bg, x) if x > -0.06 else _u(x)
        D = 3 * N**2 * V**2 - N**2 + 4 * N * V - V**2 + 3
        du = css.rhs_plain(x, [A, N, W, V])
        Dp = (6 * N * V**2 - 2 * N + 4 * V) * du[1] + (6 * N**2 * V + 4 * N - 2 * V) * du[3]
        q = css.coeffs(A, N, W, V)
        detC = float(np.linalg.det(np.linalg.solve(
            np.array([[q["a"], q["b"]], [q["c"], q["d"]]]),
            np.array([[q["sa"], q["sb"]], [q["sc"], q["sd"]]]))))
        vrel = (1 + N * V) / (N + V)
        assert -D > 0 and Dp / D < 0 and detC < 0
        assert abs(-D - 3 * (N + V) ** 2 * (1 / 3 - vrel**2)) < 1e-8 * abs(D)
        assert abs(D * detC - (3 - V * V)) < 1e-8


def _u(x):
    from p4 import css, shoot
    from p4.validated.tc_reduced_eigs import V0
    sh = shoot.shoot(V0, 1, x_end=-9.0, keep_sol=True)
    y = sh.sol.sol(x)
    N, W, V = y[0] * np.exp(-x), y[1] * np.exp(2 * x), y[2] * np.exp(x)
    return np.array([1 + css.Am1_constraint(N, W, V), N, W, V])
