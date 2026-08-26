"""S2 tests for the validated shooting integrator (Theorem A, item A3): building blocks.

Fast tests (< 1 min total): the polynomial shooting system, the Taylor/variational
coefficients, the y-tail certificate, the augmented level equations, and a short
validated integration checked against S1's DOP853 trajectory.  The full Krawczyk
run lives in test_validated_shoot_matching.py.
"""
import numpy as np
import pytest
from flint import arb
from p4 import css, shoot
from p4.validated import lintail, matching, recursion, shootsys as ss, sonic, tmint, variational as va
from p4.validated.arbseries import precision
from p4.validated.tailbound import certify_tail

V0_EC = 0.112439401388092


@pytest.fixture(scope="module")
def s1():
    """S1 profile (DOP853, scaled reduced system) and its V0-derivative by central differences."""
    h = 1e-6
    s0 = shoot.shoot(V0_EC, 1, x_end=-5.0, keep_sol=True, delta=0.05)
    sp = shoot.shoot(V0_EC + h, 1, x_end=-5.0, keep_sol=True, delta=0.05)
    sm = shoot.shoot(V0_EC - h, 1, x_end=-5.0, keep_sol=True, delta=0.05)
    return lambda x: (s0.sol.sol(x), (sp.sol.sol(x) - sm.sol.sol(x)) / (2 * h))


def test_shoot_system_matches_S1_rhs_and_jacobian(s1):
    sys4 = ss.shoot_system()
    x = -0.5
    u, _ = s1(x)
    with precision(256):
        u0 = [arb(float(v)) for v in u] + [arb(x).exp()]
        co = ss.taylor_coefficients(sys4, u0, 12)
        Y = ss.variational_coefficients(sys4, co, 4)
    assert np.allclose([float(c) for c in co[1][:3]], css.rhs_scaled3(x, u), rtol=1e-12)
    assert abs(float(co[1][3]) - np.exp(x)) < 1e-15
    eps = 1e-6
    for j in range(3):
        e = np.zeros(3)
        e[j] = eps
        fd = (css.rhs_scaled3(x, u + e) - css.rhs_scaled3(x, u - e)) / (2 * eps)
        assert np.allclose([float(Y[1][i, j]) for i in range(3)], fd, rtol=1e-6, atol=1e-8)


def test_y_tail_certificate_bounds_explicit_orders(s1):
    """The affine-contraction bound for the tail of y = du/dV0 checked against explicitly
    computed orders K+1..44 (regular point x = -2)."""
    sys4 = ss.shoot_system()
    sys7 = va.augment(sys4, skip=(3,))
    eqs4, eqs7 = ss.regular_level_equations(4), ss.regular_level_equations(7)
    x = -2.0
    u, y = s1(x)
    with precision(384):
        z0 = [arb(float(v)) for v in u] + [arb(x).exp()] + [arb(float(v)) for v in y]
        co = ss.taylor_coefficients(sys7, z0, 44, blocks=tmint.BLOCKS)
        K = 28
        co4 = [c[:4] for c in co[:K + 1]]
        D, E = recursion.structure_matrices(sys4, eqs4, co4, dm=1)
        cu = certify_tail(sys4, eqs4, co4, D, E)
        cy = lintail.linear_tail_certificate(sys7, eqs7, co[:K + 1], 4, cu)
        assert cu.ok and float(cu.nu) >= 0.08
        tail = arb(0)
        for n in range(K + 1, 45):
            yn = max(abs(c) for c in co[n][4:])
            assert yn <= cy.eps * cy.nu ** (-n)
            tail += yn * cy.nu ** n
        assert tail <= cy.eps


def test_augmented_sonic_level_equations_need_ell_prime():
    """The y-rows of the sonic level equations need the d(ell)/dV0 . F_n terms (they cancel the
    u_{n+1} contribution): without them the order-K residual is nonzero, with them all vanish."""
    w, m, K = 1e-9, 3, 12
    ex = sonic.sonic_expansion("0.112439401388092", K=K, width=w, m=m)
    sys8 = va.augment(ex.sys)
    by = va.derivative_balls(ex.point["coefs"], ex.interval["coefs"], m, w)
    coefs8 = va.augmented_coefs(ex.balls, by)
    ell = ex.interval["ell"]
    bad0 = va.check_level_residuals(sys8, va.augmented_level_equations(ex.eqs, 4), coefs8, 2)
    eqs8 = va.augmented_level_equations(ex.eqs, 4, extra={3: [(ell[0][1], 2, 0), (ell[1][1], 3, 0)]})
    assert bad0 == [K] and va.check_level_residuals(sys8, eqs8, coefs8, 2) == []


def test_sonic_initial_state_and_short_integration(s1):
    """Sonic-side state at x0 = -0.05 (Taylor model in V0, width 2e-9) integrated to x = -0.3:
    the point set contains S1's trajectory and derivative, the interval set the point set."""
    with precision(384):
        st = matching.sonic_initial_state("0.112439401388092", 1e-9, x0=-0.05)
        assert max(st.pt.widths()) < 1e-30 and max(st.iv.widths()) < 1e-6
        it = tmint.Integrator(K=28, hmax=0.02)
        it.integrate(st, -0.3)
        up, yp, yi, ui = st.u_point(), st.y_point(), st.y_interval(), st.u_interval()
    u, y = s1(-0.3)
    for i in range(3):
        assert (up[i] + arb(0, 1e-9)).contains(arb(float(u[i]))), i     # S1 accuracy ~1e-10
        assert (yp[i] + arb(0, 1e-5)).contains(arb(float(y[i]))), i     # central differences ~1e-6
        assert yi[i].contains(yp[i]) and ui[i].contains(up[i])
    assert max(st.pt.widths()) < 1e-20
    assert max(float(b.rad()) for b in yi) < 1e-4
    assert all(l["nu"] >= 0.005 for l in st.log) and len(st.log) < 100
