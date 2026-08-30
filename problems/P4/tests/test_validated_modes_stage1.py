"""Theorem B, Stage 1: linearised system (linsys) and its sonic-point series for complex kappa (linsonic)."""
import numpy as np
import pytest
from flint import acb, arb
from p4 import perturb, taylor
from p4.validated import centre, linsonic, sonic
from p4.validated.arbseries import precision
from p4.validated.linsys import LinSystem, abs_up, plain_from_scaled

V0_EC = "0.112439401388092"
KAPPA1, KGAUGE = "2.8105525488", "0.3556992037"
K = 40


@pytest.fixture(autouse=True)
def _prec():
    """Ball evaluations outside the constructors run at the ambient precision."""
    with precision(256):
        yield


@pytest.fixture(scope="module")
def bg():
    ex = sonic.sonic_expansion(V0_EC, K=K + 1)
    ex.certify()
    return ex


@pytest.fixture(scope="module")
def L():
    return LinSystem()


def _contains(ball, z, tol):
    return (ball + acb(arb(0, tol), arb(0, tol))).contains(acb(z.real, z.imag))


def _S1_matrix(x, y_plain, kap):
    """S1's dual-number linearisation (perturb.linearized_rhs) as the 3x3 matrix q' = A q, plus A_p."""
    A, N, W, V = y_plain
    e = np.exp(x)
    bgs = [N * e, W / e**2, V / e]
    cols = []
    for j in range(3):
        p = np.zeros(3, complex)
        p[j] = 1
        u = np.concatenate([bgs, [p[0].real, p[0].imag, p[1].real, p[1].imag, p[2].real, p[2].imag]])
        d = perturb.linearized_rhs(x, u, np.array([kap]))
        cols.append([d[3] + 1j * d[4], d[5] + 1j * d[6], d[7] + 1j * d[8]])
    Ap = perturb.A_p_from_constraint(x, np.array(bgs), np.array([1 + 0j]), np.array([2j]), np.array([0.5 + 0j]),
                                     np.array([kap]))[0]
    return np.array(cols).T, Ap


def test_linear_system_matches_S1_on_float_background(L):
    """(c) the acb coefficient matrix A(x; kappa) and A_p contain S1's float values (S1 sonic series bg)."""
    bgf = taylor.background_series(float(V0_EC), 1, K=40)
    kap = 2.81 + 0.7j
    for x in (-0.05, -0.02, -0.3):
        y = taylor.eval_series(bgf, x)
        Af, Apf = _S1_matrix(x, y, kap)
        Ab = L.coefficient_matrix([arb(float(v)) for v in y], acb(kap.real, kap.imag))
        for r in range(3):
            for c in range(3):
                assert _contains(Ab[r, c], Af[r, c], 1e-10 * (1 + abs(Af[r, c])))
                assert float(abs_up(Ab[r, c] - Ab[r, c].mid())) < 1e-60
        Apb = L.A_p([arb(float(v)) for v in y], acb(kap.real, kap.imag), [acb(1), acb(0, 2), acb(0.5)])
        assert _contains(Apb, Apf, 1e-12)


def test_linear_system_on_certified_backgrounds(L, bg):
    """A(x; kappa) as balls on the certified A1 Taylor model (x = -0.05, -0.03) and on the certified
    A2 centre family (x = -3, -4) contains S1's linearisation evaluated on the midpoints."""
    ce = centre.centre_expansion("5.82098013", nhat="1.2365999612", K=30)
    ce.certify()
    kap = 2.81 - 0.3j
    for x in (-0.05, -0.03, -3.0, -4.0):
        u = bg.eval(x) if x > -1 else plain_from_scaled(x, *ce.eval(x))
        Af, _ = _S1_matrix(x, [float(v) for v in u], kap)
        Ab = L.coefficient_matrix(u, acb(kap.real, kap.imag))
        for r in range(3):
            for c in range(3):
                assert _contains(Ab[r, c], Af[r, c], 1e-10 * (1 + abs(Af[r, c])))
                assert float(abs_up(Ab[r, c] - Ab[r, c].mid())) < 1e-25


def test_sonic_series_at_kappa1_contains_S1(bg):
    """(c) the acb series at kappa = kappa1 (point) contains S1's taylor.perturbation_series floats
    (orders <= 38: S1's top orders suffer its TPS truncation), and the constraint-eliminated series
    satisfies row 1 of the 4D linearisation order by order (ball identity)."""
    ex = linsonic.linear_sonic_expansion(bg, KAPPA1, width=0.0, K=K)
    bgf = taylor.background_series(float(V0_EC), 1, K=42)
    pc = taylor.perturbation_series(bgf, float(KAPPA1), K=K)
    F = ex.floats()
    balls = [ex.Ap] + [[b[i] for b in ex.balls] for i in range(3)]
    for n in range(39):
        for i in range(4):
            assert _contains(balls[i][n], pc[i, n], 1e-10 * (1 + abs(pc[i, n]))), (i, n, F[i, n], pc[i, n])
    assert ex.radii().max() < 1e-38
    assert ex.info["row1_contains_zero"] and ex.info["row1_max_abs"] < 1e-35
    assert abs(F[:, 0] - np.array([1, 0, 0.09347531, -0.37577332])).max() < 1e-8


@pytest.mark.parametrize("kc", [KAPPA1, KGAUGE])
def test_box_expansion_certified_and_contains_interior_point(bg, kc):
    """Box kappa_c +- 1e-6 (both parts): tail certificate on |x| <= 0.05; the Taylor-model balls
    contain the point series at an interior complex kappa (coefficients and values at x = -0.05)."""
    ex = linsonic.linear_sonic_expansion(bg, kc, width=1e-6, m=5, K=K)
    cert = ex.certify()
    assert cert.ok and float(cert.nu) >= 0.05 and float(cert.tail_bound(0.05)) < 1e-40
    kin = acb(arb(kc) + arb("5e-7"), arb("3e-7"))
    pt = linsonic.linear_sonic_expansion(bg, kin, width=0.0, K=K)
    pt.certify()
    for n in range(K + 1):
        assert ex.Ap[n].contains(pt.Ap[n])
        for i in range(3):
            assert ex.balls[n][i].contains(pt.balls[n][i])
    vb, vp = ex.eval(-0.05), pt.eval(-0.05)
    assert all(vb[i].contains(vp[i]) for i in range(4))
    assert all(float(abs_up(vb[i] - vb[i].mid())) < 1e-4 for i in range(4))
    assert ex.info["row1_contains_zero"]
    if kc == KGAUGE:
        assert ex.kappa.contains(acb(linsonic.gauge_eigenvalue(bg)))


def test_tail_bound_consistent_with_higher_orders(bg):
    """|q_n| <= eps nu^{-n} for the exactly computed orders 41..50 (K = 50 run vs K = 40 certificate)."""
    big = sonic.sonic_expansion(V0_EC, K=51)
    big.certify()
    ex40 = linsonic.linear_sonic_expansion(bg, KAPPA1, width=0.0, K=K)
    cert = ex40.certify()
    ex50 = linsonic.linear_sonic_expansion(big, KAPPA1, width=0.0, K=50)
    for n in range(K + 1, 51):
        for i in range(3):
            assert abs_up(ex50.balls[n][i]) <= cert.coef_bound(n)


def test_resonance_is_detected():
    """A box containing the first resonance kappa = -0.0990 - 1.0990 (M_1 singular) is rejected."""
    bg3 = sonic.sonic_expansion(V0_EC, K=6)
    with pytest.raises(ZeroDivisionError):
        linsonic.linear_sonic_expansion(bg3, "-1.198", width=5e-3, m=1, K=4, check=False)
