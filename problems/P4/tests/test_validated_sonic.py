"""S2 tests for the validated sonic-point expansion (Theorem A, item A1)."""
import numpy as np
import pytest
from flint import arb
from p4 import css, shoot, taylor
from p4.validated import sonic
from p4.validated.arbseries import precision

V0_EC = "0.112439401388092"
K = 40


@pytest.fixture(scope="module")
def ex_point():
    return sonic.sonic_expansion(V0_EC, K=60)


@pytest.fixture(scope="module")
def ex_interval():
    return sonic.sonic_expansion(V0_EC, K=K, width=5e-11, m=3)


def test_zeroth_and_first_order_data(ex_point):
    A0, N0, W0, V0 = ex_point.balls[0]
    assert A0.contains(arb("1.8614267226")) or abs(float(A0) - 1.8614267226) < 1e-9
    assert abs(float(N0) - 2.0113168913) < 1e-9 and abs(float(W0) - 0.3256888894) < 1e-9
    A1, N1, W1, V1 = [float(c) for c in ex_point.balls[1]]
    assert np.allclose([A1, N1, W1, V1], [-0.37029024, -0.71542382, -0.48249731, 0.48472202], atol=2e-8)
    roots = ex_point.info["roots"]
    assert float(roots[0]) < 0 < float(roots[1])                   # two real branches, EC has V1 > 0
    # S1 closed-form quadratic is proportional to the extracted one
    c2, c1, c0 = css.first_order_quadratic(float(V0_EC))
    e2, e1, e0 = [float(ex_point.info[k]) for k in ("c2", "c1", "c0")]
    assert abs(e2 * c1 - e1 * c2) < 1e-12 * abs(e2 * c1) and abs(e2 * c0 - e0 * c2) < 1e-12 * abs(e2 * c0)
    assert float(ex_point.info["Delta1"]) > 0                        # simple zero of Delta~ on the EC branch
    assert ex_point.info["constraint_exponent"].contains(arb(0))     # gamma = 0: constraint propagates


def test_balls_contain_S1_float_coefficients(ex_interval):
    """(i) the ball coefficients (V0 interval of width 1e-10) contain S1's taylor.py floats."""
    ref = taylor.background_series(float(V0_EC), 1, K=K + 1)
    for n in range(K + 1):
        for i in range(4):
            assert ex_interval.balls[n][i].contains(arb(ref[i, n])), (n, i)
    rad = ex_interval.radii()
    assert rad[:, 0].max() < 3e-10 and rad[:, K].max() < 1e-7        # width-driven, not blown up
    assert max(max(r) for r in ex_interval.rem) < 1e-9                # Lagrange remainders (m = 3)


def test_point_coefficients_match_S1_to_float_accuracy(ex_point):
    ref = taylor.background_series(float(V0_EC), 1, K=61)
    fl = ex_point.floats()
    assert np.abs(fl[:, :41] - ref[:, :41]).max() < 1e-12
    assert ex_point.radii()[:, 60].max() < 1e-30


def test_series_encloses_S1_ode_solution_at_minus_0p05(ex_interval):
    """(ii) truncated ball series at x = -0.05 vs S1's DOP853 integration from x = -0.1.

    The reference is the integrator's *endpoint* value (rtol 1e-13); DOP853's dense-output
    interpolant is only accurate to ~1e-10 over the large steps a smooth solution allows."""
    x = -0.05
    sh = shoot.shoot(float(V0_EC), 1, x_end=x)
    assert sh.reason == "end" and sh.x_stop == x
    ref = css.from_scaled(x, sh.y_end)                                # (A, N, W, V)
    with precision(256):
        vals = ex_interval.eval(x, with_tail=False)
    for v, r in zip(vals, ref):
        assert (v + arb(0, 1e-11)).contains(arb(float(r)))            # 1e-11: ODE tolerance of the reference


def test_tail_certificate_and_majorant(ex_point):
    """(iii) certificate at K = 40 checked against the exactly computed orders 41..60."""
    ex40 = sonic.sonic_expansion(V0_EC, K=40)
    cert = ex40.certify()
    assert cert.ok and float(cert.nu) >= 0.08 and float(cert.Z1) < 1
    nu, eps = cert.nu, cert.eps
    with precision(256):
        tail = arb(0)
        for n in range(41, 61):
            un = max(abs(c) for c in ex_point.balls[n])
            assert (un <= eps * nu ** (-n))                            # |u_n| <= eps nu^-n
            tail += un * nu**n
        assert tail <= eps                                             # sum_{n>40} |u_n| nu^n <= eps
        assert float(cert.tail_bound(arb("0.05"))) < 1e-40
        # enclosure with tail at x = -0.05 contains the K = 60 evaluation
        v40 = ex40.eval(-0.05)
        v60 = ex_point.eval(-0.05, with_tail=False)
        assert all(a.overlaps(b) for a, b in zip(v40, v60))


def test_interval_certificate_and_enclosure_at_minus_0p05(ex_interval):
    """The tail certificate holds uniformly over the V0 interval (same nu as at the point)."""
    cert = ex_interval.certify()
    assert cert.ok and float(cert.nu) >= 0.1 and float(cert.tail_bound(arb("0.05"))) < 1e-50
    with precision(256):
        A, N, W, V = ex_interval.eval(-0.05)
    assert A.contains(arb("1.8787434870")) and V.contains(arb("0.0884703313"))
    assert max(float(b.rad()) for b in (A, N, W, V)) < 3e-10


def test_other_branch_and_general_width():
    other = sonic.sonic_expansion(V0_EC, K=6, branch="other")
    assert float(other.balls[1][3]) < 0
    wide = sonic.sonic_expansion(V0_EC, K=12, width=1e-6, m=4)
    assert wide.radii()[:, 12].max() < 1e-3 and max(max(r) for r in wide.rem) < 1e-6
