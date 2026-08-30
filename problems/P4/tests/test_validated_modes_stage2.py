"""Theorem B, Stage 2: validated propagation of the linearised system along the certified
background (linscaled, linstep, lintube, linprop), the centre condition (lincentre) and the
matching function E(kappa) with its kappa-derivative (linmatch).

Fast tests use short tubes; the full pipeline (tube -0.05 -> -3 -> -8, ~10 min) is built once by
the module fixture ``pipeline`` unless P4_TUBE_CACHE points to a tube saved by ``Tube.save``."""
import os
import time

import numpy as np
import pytest
from flint import acb, acb_mat, arb
from p4 import perturb
from p4.validated import centre, lincentre, linmatch, linprop, linscaled, linsonic, linstep, lintube, sonic
from p4.validated import shootsys as ss
from p4.validated.arbseries import precision
from p4.validated.linsys import LinSystem, abs_up, plain_from_scaled

V0_EC = "0.1124394013880983"          # c + 6.3e-15: centre of the certified V0* enclosure (A3, results/a3_midpoint.json)
W_V0 = 1e-16                          # V0* in [c* +/- 1e-16]  (K(X) - m = 6.3e-15 +/- 7.5e-17)
A_STAR = "-0.2123656467659762832750918714540807905889"
MU_STAR = "8.901323275379966931515526907200000000000"
KAPPA1, KGAUGE = "2.8105525488", "0.3556992037"


@pytest.fixture(autouse=True)
def _prec():
    with precision(256):
        yield


def certified_centre():
    a = arb(A_STAR) + arb("7.031e-12") + arb(0, 8.87e-16)
    mu = arb(MU_STAR) + arb("2.3757e-9") + arb(0, 4.61e-14)
    ce = centre.centre_expansion(mu * (2 * a).exp(), nhat=(-a).exp(), K=30)
    assert ce.certify().ok
    return ce


def rad(z):
    return float(abs_up(z - z.mid()))


# ---------------------------------------------------------------------------------------------
# fast: the scaled systems and the centre family
# ---------------------------------------------------------------------------------------------
def test_scaled_systems_match_stage1_and_each_other():
    """3D scaled system = Lambda^{-1} A Lambda - Gamma of Stage 1; the 4D system restricted to the
    linearised constraint surface reproduces the 3D one (exact polynomial forms, ball check)."""
    L = LinSystem()
    S3, S4 = linscaled.reduced_system(L), linscaled.full_system(L)
    sys4 = ss.shoot_system()
    x = arb("-1")
    z = [arb("1.1"), arb("2.3"), arb("-0.35"), x.exp()]
    f, _ = ss.rhs_enclosure(sys4, z, lintube.BLOCKS4)
    kap = acb("2.3", "0.7")
    u = plain_from_scaled(x, *z[:3])
    A = L.coefficient_matrix(u, kap)
    T, Lam, Gam = x.exp(), [1 / x.exp(), x.exp() ** 2, x.exp()], [-1, 2, 1]
    At = [[A[r, c] * Lam[c] / Lam[r] - (Gam[r] if r == c else 0) for c in range(3)] for r in range(3)]
    Pj, Gj = S3.box_matrices(z, f[:3])
    P = acb_mat([[acb(Pj[0][r][c]) + kap * Pj[1][r][c] for c in range(3)] for r in range(3)])
    G = acb_mat([[acb(Gj[0][r][c]) + kap * Gj[1][r][c] + kap**2 * Gj[2][r][c] for c in range(3)] for r in range(3)])
    As = P.solve(G)
    assert all(As[r, c].overlaps(At[r][c]) and rad(As[r, c]) < 1e-60 for r in range(3) for c in range(3))
    Pj4, Gj4 = S4.box_matrices(z, f[:3])
    A4 = acb_mat([[acb(Pj4[0][r][c]) for c in range(4)] for r in range(4)]).solve(
        acb_mat([[acb(Gj4[0][r][c]) + kap * Gj4[1][r][c] for c in range(4)] for r in range(4)]))
    args = [arb(0)] + u
    Cq = [ss.eval_box(L.dC[l + 1], args) for l in range(3)]
    Sv = ss.eval_box(L.S, args)
    for c in range(3):
        qt = [acb(int(i == c)) for i in range(3)]
        Ap = u[0] * sum((Cq[l] * qt[l] * Lam[l] for l in range(3)), acb(0)) / ((kap - u[0]) * Sv)
        p = [Ap / (T * T)] + qt
        for r in range(3):
            d4 = sum((A4[r + 1, l] * p[l] for l in range(4)), acb(0))
            d3 = sum((As[r, l] * qt[l] for l in range(3)), acb(0))
            assert (d4 - d3).contains(acb(0)) and rad(d4 - d3) < 1e-60


def test_centre_regular_family_certified_and_on_constraint():
    """4D leading matrix at the centre: exponents {-3,-3,0,0} for every kappa (incl. kappa = 1 where
    the 3D reduced form degenerates); the regular family certifies with nu = 0.06 and its members
    satisfy the linearised momentum constraint at t = e^{-3} (ball identity)."""
    ce = certified_centre()
    S4, S3 = linscaled.full_system(), linscaled.reduced_system()
    assert lincentre._rows_vanish_at_centre(S4) == [1, 3]
    z0 = list(ce.balls[0]) + [arb(0)]
    L = LinSystem()
    for kap in (acb(0), acb(1), acb("2.8105525488"), acb(1, 0.5), acb(5, 7)):
        Pj, Gj = S4.box_matrices(z0, [arb(0)] * 3)
        P = np.array([[complex(Pj[0][r][c]) for c in range(4)] for r in range(4)])
        G = np.array([[complex(Gj[0][r][c] + kap * Gj[1][r][c]) for c in range(4)] for r in range(4)])
        ev = np.sort_complex(np.linalg.eigvals(np.linalg.solve(P, G)))
        assert np.allclose(ev, [-3, -3, 0, 0], atol=1e-8)
        rf = lincentre.RegularFamily(S4, ce, kap, K=50)
        ok, eps, det = rf.certify(0.06)
        assert ok and det["Z1"] + det["Z2"] < 0.6 and det["g"] < 6 and max(det["eps"]) < 1e-14
        x = arb(-3)
        T = x.exp()
        u = plain_from_scaled(x, *ce.eval(x))
        args = [arb(0)] + u
        for r in rf.eval(x):
            q = [r[1] / T, r[2] * T * T, r[3] * T]
            Ap = r[0] * T * T
            C = (kap - u[0]) * ss.eval_box(L.S, args) * Ap - u[0] * sum((ss.eval_box(L.dC[l + 1], args) * q[l] for l in range(3)), acb(0))
            assert C.contains(acb(0)) and rad(C) < 1e-6
    with pytest.raises((ZeroDivisionError, ValueError)):
        lincentre.RegularFamily(S3, ce, acb(1), K=10)     # 3D reduced form: P_0 = 0 at kappa = 1


def test_short_tube_propagation_encloses_stage1_series():
    """Validated propagation of the linearised sonic solution from -0.05 to -0.08 along a short
    certified tube (point V0 background) lands inside the Stage-1 series' own enclosure at -0.08."""
    tube = lintube.Tube.build(V0_EC, W_V0, x_c=-0.08)
    S3 = linscaled.reduced_system()
    bg = sonic.sonic_expansion(V0_EC, K=41)
    bg.certify()
    ex = linsonic.linear_sonic_expansion(bg, KAPPA1, width=0.0, K=40)
    ex.certify()

    def scaled(x):
        q, T = ex.eval(x), arb(x).exp()
        return [q[1] * T, q[2] / (T * T), q[3] / T]
    st = linstep.LohnerSet(scaled(-0.05))
    sets, log = linprop.propagate(tube, S3, KAPPA1, [st])
    ref = scaled(log[-2]["x"])
    for a, b in zip(sets[0].hull(), ref):
        assert a.overlaps(b) and rad(a) < 1e-9
    assert all(l["bound"] < 1e-11 and l["Dsup"] < 1e-30 for l in log[:-1])
    Phi, _ = linprop.fundamental_matrix(tube, S3, KAPPA1)
    q0, qe = scaled(-0.05), scaled(log[-2]["x"])
    for r in range(3):
        assert sum((Phi[r, c] * q0[c] for c in range(3)), acb(0)).overlaps(qe[r])


# ---------------------------------------------------------------------------------------------
# the full pipeline (slow: ~10 min without a cached tube)
# ---------------------------------------------------------------------------------------------
@pytest.fixture(scope="module")
def pipeline():
    with precision(256):
        cache = os.environ.get("P4_TUBE_CACHE")
        ce = certified_centre()
        if cache and os.path.exists(cache):
            tube = lintube.Tube.load(cache)
        else:
            tube = lintube.Tube.build(V0_EC, W_V0, x_c=-3.0)
        if tube.steps[-1].region != "centre":
            tube.extend_centre(ce, -8.0, hmax=0.05)
        bg = linmatch.box_background(V0_EC, W_V0)
        M = linmatch.Matcher(tube, ce, bg)
        prob = perturb.Problem(0.112439401388092, x_end=float(M.x_d))
        return M, prob


def _close(ball, z, tol):
    return (ball + acb(arb(0, tol), arb(0, tol))).contains(acb(z.real, z.imag))


def test_E_fin_matches_S1_and_E_has_S1_signs(pipeline):
    """(d) E_fin = e^{x_d} A_p(x_d) encloses S1's float E at x_end = x_d (S1 accuracy ~1e-8) at
    kappa = 2.5, 3.0, 1 +/- 0.5i; the analytic centre determinant E has the same sign pattern
    (real for real kappa, conjugate-symmetric); ball radii <= 1e-6 and ~2 s per kappa once the
    tube data are derived (the first call derives them: ~1 min)."""
    M, prob = pipeline
    res = [M.E(k) for k in ("2.5", "3.0", acb(1, "0.5"), acb(1, "-0.5"))]
    S1 = prob.E(np.array([2.5, 3.0, 1 + 0.5j, 1 - 0.5j]))
    for r, e in zip(res, S1):
        assert _close(r["E_fin"], e, 1e-7), (r["E_fin"], e)
        assert rad(r["E_fin"]) < 1e-6 and rad(r["E"]) < 1e-6 and r["t_prop"] < 180
    assert max(r["t_prop"] for r in res[1:]) < 20
    e25, e30 = res[0]["E"].real, res[1]["E"].real
    assert (e25 < 0 < e30) or (e25 > 0 > e30)
    assert (res[0]["E_fin"].real < 0 < res[1]["E_fin"].real)
    c1, c2 = res[2]["E"], res[3]["E"]
    assert c1.real.overlaps(c2.real) and c1.imag.overlaps(-c2.imag) and not c1.imag.contains(arb(0))


def test_E_zero_at_kappa1_and_gauge_and_box_contains_point(pipeline):
    """(d) E(kappa1) and E(kappa_gauge) contain 0 (S1's zeros, accurate to 1e-9), with the enclosed
    dE/dkappa the zero of E is located within |E|/|dE| of them; the kappa-box enclosure (Lipschitz
    perturbation bound, w = 1e-8) contains the point value."""
    M, prob = pipeline
    pt = M.E(KAPPA1, deriv=True)
    assert pt["E"].contains(acb(0)) and rad(pt["E"]) < 2e-9 and rad(pt["E_fin"]) < 5e-9
    assert _close(pt["E_fin"], prob.E(np.array([float(KAPPA1)]))[0], 1e-8)
    assert abs(pt["dE"]) > 0.02 and rad(pt["dE"]) < 1e-5 and abs(pt["E"] / pt["dE"]) < 1e-7
    g = M.E(KGAUGE, deriv=True)
    assert g["E"].contains(acb(0)) and abs(g["E"] / g["dE"]) < 1e-5
    bx = M.E(KAPPA1, 1e-8)
    assert bx["E"].contains(pt["E"]) and bx["E_fin"].contains(pt["E_fin"]) and rad(bx["E"]) < 1e-2


def test_apparent_singularity_needs_4D(pipeline):
    """(c) real kappa = 1.2 in (A(-3), A_0): the reduced 3D system has (kappa - A(x)) = 0 near x = -1.47
    and its validated propagation breaks down there (singular box matrix or exploding bound), while
    the 4D system (A_p kept) propagates through and E_fin encloses S1's E; at kappa = 2.5 the 4D and
    3D results agree.  (kappa = 1.5 is not even certified at the sonic point by the 3D form.)"""
    M, prob = pipeline
    try:
        r3 = M.E("1.2")
        broken = r3["width_end"] > 1e3 or not np.isfinite(r3["width_end"])
    except (ZeroDivisionError, ValueError, RuntimeError):
        broken = True
    assert broken
    r4 = M.E("1.2", system="4D")
    assert _close(r4["E_fin"], prob.E(np.array([1.2]))[0], 2e-6) and rad(r4["E"]) < 1e-6 * abs(r4["E"])
    r3, r4b = M.E("2.5"), M.E("2.5", system="4D")
    assert r3["E"].overlaps(r4b["E"]) and r3["E_fin"].overlaps(r4b["E_fin"])
