"""P9(b) S1: exact identity tests for the derivation chain (deliverables i-vii).

Run:  uv run pytest problems/P9/b/tests -q     (from the repo root)
Every assertion is an exact sympy zero-remainder unless marked numeric.
"""
import pytest

sp = pytest.importorskip("sympy")

from derivation import actions as A
from derivation import background as B
from derivation import dictionary as D
from derivation import geometry as G
from derivation import ibp
from derivation import tools as T


# --- geometry: GLV14 (73)-(75) as known answers ------------------------------------
@pytest.fixture(scope="module")
def geo():
    return G.build()


def test_glv74_delta_K(geo):
    """delta K^i_j = (zetadot - H dN) delta^i_j - d^i d_j psi_l  (psi_l = a^2 psi;
    GLV14 (74) uses the lower-index shift potential N_i = d_i psi_l while their
    (73) writes N^i = delta^{ij} d_j psi — an internal a^2 mismatch in the
    paper's notation, physically irrelevant since psi is auxiliary)."""
    dK = [[sp.expand(T.bgsubs(geo["Kud"][i][j]).coeff(T.eps, 1)) for j in range(3)]
          for i in range(3)]
    zd = sp.Derivative(T.zf, T.t)
    diag = sp.expand(zd * T.cx - T.H * T.dnf * T.cx)
    assert sp.simplify(dK[0][0] - diag - T.k**2 * T.psf * T.cx) == 0
    assert sp.simplify(dK[1][1] - diag) == 0 and sp.simplify(dK[2][2] - diag) == 0
    assert all(dK[i][j] == 0 for i in range(3) for j in range(3) if i != j)


def test_glv74_delta_sqrth(geo):
    assert sp.simplify(sp.expand(T.bgsubs(geo["sqh"])).coeff(T.eps, 1)
                       - 3 * T.a**3 * T.zf * T.cx) == 0


def test_glv75_ricci3(geo):
    z = T.eps * T.zf * T.cx
    R3 = sp.expand(T.bgsubs(geo["R3"]))
    e1 = sp.expand(-(4 / T.a**2) * sp.diff(z, T.x, 2) / T.eps)
    e2 = sp.expand(-(2 / T.a**2) * (sp.diff(z, T.x)**2
                                    - 4 * z * sp.diff(z, T.x, 2)) / T.eps**2)
    assert sp.simplify(R3.coeff(T.eps, 1) - e1) == 0
    assert sp.simplify(R3.coeff(T.eps, 2) - e2) == 0


def test_adm_inverse_and_r4_background(geo):
    for m in range(4):
        for n in range(4):
            s = T.cut(sum(T.mul(geo["g4"][m][l], geo["g4i"][l][n]) for l in range(4)))
            assert sp.simplify(s - (1 if m == n else 0)) == 0
    R40 = T.bgsubs(G.ricci4()).coeff(T.eps, 0)
    assert sp.simplify(R40 - 6 * (sp.Derivative(T.H, T.t) + 2 * T.H**2)) == 0


# --- background: GLV14 (88)-(89), matter eom, family B -----------------------------
@pytest.fixture(scope="module")
def tads():
    return B.tadpole_equations()


def test_glv88_glv89(tads):
    E_dn, E_z, _ = tads
    fd, fdd = sp.Derivative(T.f_, T.t), sp.Derivative(T.f_, (T.t, 2))
    Hd = sp.Derivative(T.H, T.t)
    eq88 = T.cc + T.Lam - (3 * T.Ms2 * (T.f_ * T.H**2 + fd * T.H) - B.rhom_of_P)
    eq89 = T.Lam - T.cc - (T.Ms2 * (2 * T.f_ * Hd + 3 * T.f_ * T.H**2
                                    + 2 * fd * T.H + fdd) + B.pm_of_P)
    assert sp.cancel(sp.expand(E_dn + T.a**3 * eq88)) == 0
    assert sp.cancel(sp.expand(E_z + 3 * T.a**3 * eq89)) == 0


def test_matter_eom_and_conservation(tads):
    _, _, E_ds = tads
    sgd = T.sgb.diff(T.t)
    # E_ds = -d/dt(a^3 sigma' P1)  => sigma'' (P1 + 2 Yb P2) = -3 H sigma' P1
    eom = sp.solve(E_ds, sp.Derivative(T.sgb, (T.t, 2)))[0]
    assert sp.cancel(sp.together(eom + 3 * T.H * sgd * B.cm2_of_P)) == 0
    # energy conservation: rhodot + 3H(rho+p) == 0 on the eom
    cons = T.bgsubs(sp.diff(B.rhom_of_P, T.t)) + 3 * T.H * (B.rhom_of_P + B.pm_of_P)
    cons = cons.subs(sp.Derivative(T.sgb, (T.t, 2)), eom)
    assert sp.cancel(sp.together(sp.expand(cons))) == 0


def test_family_B_de_conservation():
    # d/dt rho_DE + 3H(rho_DE + p_DE) == 0 given dust conservation (definitions)
    rd = sp.Function("rho_m")(T.t)
    rho = B.rho_DE.subs(T.rhom, rd)
    f = B.f_DE.subs({T.rhom: rd, T.pm: 0})
    cons = sp.diff(rho, T.t) + 3 * T.H * f
    cons = cons.subs(sp.Derivative(rd, T.t), -3 * T.H * rd)
    assert sp.simplify(cons) == 0


def test_kessence_unitary_operators():
    c_val, lam_val, M24_val = B.kessence_unitary()
    Xb = T.phb.diff(T.t)**2 / 2
    assert sp.cancel(c_val - Xb * T.K1) == 0
    assert sp.cancel(lam_val - (Xb * T.K1 - T.K0)) == 0
    assert sp.cancel(M24_val - Xb**2 * T.K2) == 0


# --- GPV dictionary (deliverable v) and f(R) limit (vi) ----------------------------
def test_gpv_dictionary_exact():
    """[L87 + Lmat]_eps2 == [L86|dict + Lmat]_eps2 + remainder + d/dt F, with the
    remainder exactly a^3 (3/2 rho_m zeta dN - 9/4 p_m zeta^2) (x-averaged)."""
    L87q = sp.expand(T.bgsubs(A.quad(A.L87_jet())))
    L86d = T.subs_fun(sp.expand(T.bgsubs(A.L86())), D.DICT)
    def csub(e):
        for fn, val in ((T.cc, D.C_SOL), (T.Lam, D.LAM_SOL)):
            for o in (2, 1):
                e = e.subs(sp.Derivative(fn, (T.t, o)), T.bgsubs(sp.diff(val, (T.t, o))))
            e = e.subs(fn, val)
        return sp.expand(T.bgsubs(e))
    diff = csub(L87q - L86d)
    Lc, F = ibp.canon(diff, [T.zf, T.dnf, T.psf])
    assert ibp.check_canon(diff, Lc, F)
    resid = B.fluid(sp.cancel(sp.together(csub(Lc))))
    expected = (sp.Rational(3, 2) * T.a**3 * T.rhom * T.zf * T.dnf
                - sp.Rational(9, 4) * T.a**3 * T.pm * T.zf**2)
    assert sp.cancel(sp.together(sp.expand(resid - expected))) == 0


def test_fR_limit():
    aK_l, aBS_plus_aM, aT_l, aH_l = D.frow_limit()
    assert aK_l == 0 and aT_l == 0 and aH_l == 0
    assert aBS_plus_aM == 0          # alpha_B^BS = -alpha_M  (digest D10)
