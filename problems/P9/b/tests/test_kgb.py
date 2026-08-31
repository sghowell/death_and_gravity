"""P9(b) S1: KGB independent route (symbolic) + DPSV numeric known-answer."""
import pytest

sp = pytest.importorskip("sympy")

from derivation import actions as A
from derivation import ibp, kgb
from derivation import tools as T


@pytest.fixture(scope="module")
def bg():
    return kgb.background()


def test_box_background():
    pd, pdd = T.phb.diff(T.t), T.phb.diff(T.t, 2)
    box0 = sp.expand(T.bgsubs(kgb.box_phi())).coeff(T.eps, 0)
    assert sp.simplify(box0 + pdd + 3 * T.H * pd) == 0


def test_kgb_background_rho_p(bg):
    _, _, rho, p, _ = bg
    pd, pdd = T.phb.diff(T.t), T.phb.diff(T.t, 2)
    Xb = pd**2 / 2
    assert sp.cancel(rho - (2 * Xb * T.K1 - T.K0 - 3 * T.G1 * T.H * pd**3)) == 0
    assert sp.cancel(p - (T.K0 + T.G1 * pd**2 * pdd)) == 0


def test_dpsv_40_41_shift_charge(bg):
    """Bianchi: rhodot + 3H(rho+p) = +phidot (Jdot + 3HJ) exactly, with
    J = phidot (K_X - 3 H phidot G_X)  [DPSV (41) under G_DPSV = -G_ours]."""
    _, _, rho, p, J = bg
    pd = T.phb.diff(T.t)
    lhs = T.bgsubs(sp.diff(rho, T.t)) + 3 * T.H * (rho + p)
    rhs = pd * (T.bgsubs(sp.diff(J, T.t)) + 3 * T.H * J)
    assert sp.cancel(sp.together(sp.expand(lhs - rhs))) == 0


def test_kgb_quadratic_matches_frozen_action(bg):
    """The covariant KGB eps^2 action equals frozen L86 with
    M^2 = MPl^2, alpha_T = alpha_H = 0 and the DERIVED
      alpha_K = [2X(K_X + 2X K_XX) - 12 phid X H (G_X + X G_XX)]/(M^2 H^2),
      alpha_B^GLV = phid X G_X/(H M^2)  (i.e. alpha_B^BS = -2 phid X G_X/(H M^2)),
    matching BS14 (A.7) restricted to G_3 = -G (known-answer)."""
    E_dn, E_z, rho, p, _ = bg
    pd = T.phb.diff(T.t)
    Hd = sp.Derivative(T.H, T.t)
    K0_sol = sp.solve(sp.Eq(3 * T.MPl2 * T.H**2, rho), T.K0)[0]
    Hd_sol = sp.solve(sp.Eq(T.MPl2 * (2 * Hd + 3 * T.H**2), -p), Hd)[0]
    def bsub(e):
        return sp.expand(T.bgsubs(e.subs(Hd, Hd_sol).subs(T.K0, K0_sol)))
    M2e, aK_e, aB_BS_e = kgb.expected_alphas()
    L86e = T.subs_fun(sp.expand(T.bgsubs(A.L86())),
                      {T.M2: M2e, T.aK: aK_e, T.aBg: -aB_BS_e / 2,
                       T.aT: sp.Integer(0), T.aHh: sp.Integer(0)})
    diff = bsub(kgb.quad() - L86e)
    Lc, F = ibp.canon(diff, [T.zf, T.dnf, T.psf])
    assert ibp.check_canon(diff, Lc, F)
    assert sp.cancel(sp.together(bsub(Lc))) == 0


# --- numeric DPSV reproduction (digest D2; known-answer, loose float check) --------
def test_dpsv_stable_phantom_crossing():
    np = pytest.importorskip("numpy")
    from p9b import kgb_check as KC
    res = KC.run(pd0=0.5, om0=0.9, n_efolds=8.0)
    rep = KC.crossing_report(res)
    assert rep is not None, "no w = -1 crossing found"
    assert rep["w_before"] > -1 > rep["w_after"]
    assert rep["rho_de_c"] > 0
    assert rep["D_min"] > 0.1
    assert rep["cs2_min"] > 0.1
    # Friedmann constraint monitored: rho_phi + rho_m - 3 H^2 == 0 along the run
    import numpy as _np
    err = _np.max(_np.abs(KC.rho_phi(res["pd"], res["h"]) + res["rhom"]
                          - 3 * res["h"]**2))
    assert err < 1e-6
