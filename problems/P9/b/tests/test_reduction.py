"""P9(b) S1: reduction known-answer tests — GLV14 (77), (79)-(80), (83)-(85),
BS14 (3.12)-(3.13), the convention map alpha_B^BS = -2 alpha_B^GLV, I1, and
the k-essence dictionary."""
import pytest

sp = pytest.importorskip("sympy")

from derivation import actions as A
from derivation import ibp, rows
from derivation import reduction as R
from derivation import tools as T


# --- vacuum single-field reduction (GLV14 77, 79-80, 83; BS14 3.12) ----------------
@pytest.fixture(scope="module")
def vacuum():
    L2 = sp.expand(T.bgsubs(A.L86()))
    E_dn = sp.expand(sp.diff(L2, T.dnf))
    E_ps = sp.expand(sp.diff(L2, T.psf))
    sol = sp.solve([E_dn, E_ps], [T.dnf, T.psf], dict=True)
    assert len(sol) == 1
    Lred = sp.expand(sp.cancel(sp.together(
        L2.subs([(T.dnf, sol[0][T.dnf]), (T.psf, sol[0][T.psf])]))))
    Lc, F = ibp.canon(Lred, [T.zf])
    assert ibp.check_canon(Lred, Lc, F)
    return sol[0], Lc


def test_glv77_lapse_solution(vacuum):
    sol, _ = vacuum
    zd = sp.Derivative(T.zf, T.t)
    assert sp.cancel(sp.together(sol[T.dnf]) - zd / (T.H * (1 + T.aBg))) == 0


def test_glv79_80_and_bs312_normalisation(vacuum):
    _, Lc = vacuum
    zd = sp.Derivative(T.zf, T.t)
    Akin = sp.cancel(sp.expand(Lc).coeff(zd, 2))
    Cpot = sp.cancel(sp.expand(Lc).coeff(T.zf, 2).subs(zd, 0))
    # (80) kinetic: a^3/4 * M2 alpha/(1+aB)^2; equals a^3 Q_S/2 under the BS map
    assert not Akin.has(T.k)
    assert sp.cancel(Akin - T.a**3 * T.M2 * R.alpha / (4 * (1 + T.aBg)**2)) == 0
    assert sp.cancel(Akin - T.a**3 * R.QS_BS / 2) == 0     # (3.12)-(3.13) Q_S
    # (80) gradient: a k^2/4 * Lgrad; no k^0 mass term (adiabatic mode exact)
    Ck = sp.cancel(Cpot * 4 * T.a**2 / (T.a**3 * T.k**2))
    assert not Ck.has(T.k)
    assert sp.cancel(sp.expand(Ck - R.Lgrad_glv80)) == 0
    assert sp.cancel(Cpot - sp.expand(Cpot).coeff(T.k, 2) * T.k**2) == 0


def test_glv83_no_ghost_form():
    # (83): kinetic positivity <=> alpha = alpha_K + 6 alpha_B^2 > 0; under the
    # map this is D = alpha_K^BS-basis + (3/2)(alpha_B^BS)^2 (FORMULATION D2)
    aBS = -2 * T.aBg
    assert sp.expand(R.alpha - (T.aK + sp.Rational(3, 2) * aBS**2)) == 0


# --- two-field system with matter --------------------------------------------------
@pytest.fixture(scope="module")
def sum_prod():
    return {m: R.sum_product(m) for m in ("ah0", "gen")}


def test_bs313_exact_for_alphaH_zero(sum_prod):
    """UV speeds factorise exactly as {c_m^2, c_s^2(3.13)} when alpha_H = 0.
    This simultaneously verifies BS14 (3.13), GLV14 (85)|aH=0, and the
    convention map alpha_B^BS = -2 alpha_B^GLV (deliverable iv)."""
    S0, P0 = sum_prod["ah0"]
    tgt = R.cs2_bs313.subs(sp.Derivative(T.aHh, T.t), 0).subs(T.aHh, 0)
    assert sp.cancel(sp.together(sp.expand(S0 - T.cm2 - tgt))) == 0
    assert sp.cancel(sp.together(sp.expand(P0 - T.cm2 * tgt))) == 0


def test_glv85_product_rule_general(sum_prod):
    """c_1^2 c_2^2 = c_m^2 * c_s^2(85) exactly, for all five alphas free."""
    _, P = sum_prod["gen"]
    assert sp.cancel(sp.together(sp.expand(P - T.cm2 * R.cs2_glv85))) == 0


def test_dust_exact_speed_and_mixing_term(sum_prod):
    """Dust (c_m^2 -> 0, p_m = 0): exact DE eigen-speed differs from (85) by
    + alpha_H^2 rho_m/(M^2 H^2 alpha)  — kinetic matter mixing; equality holds
    iff alpha_H = 0."""
    cs2d = R.cs2_dust_exact()
    delta = sp.cancel(sp.together(sp.expand(
        cs2d - R.cs2_glv85.subs({T.cm2: 0, T.pm: 0}))))
    expected = T.rhom * T.aHh**2 / (R.alpha * T.H**2 * T.M2)
    assert sp.cancel(sp.together(sp.expand(delta - expected))) == 0


# --- I1 (deliverable iii) and k-essence (vii) --------------------------------------
def test_I1_exact():
    lhs, rhs = rows.I1_lhs(), rows.I1_rhs()
    assert sp.cancel(sp.together(sp.expand(lhs - rhs))) == 0


def test_kessence_dictionary():
    cs2, cs2_expected, rho_plus_p, aK_kess = rows.kessence_check()
    assert sp.cancel(sp.together(sp.expand(cs2 - cs2_expected))) == 0
    Xb = T.phb.diff(T.t)**2 / 2
    # alpha_K = (2 X P_X + 4 X^2 P_XX)/(M^2 H^2)  [GLV14 text after (85)]
    assert sp.cancel(aK_kess - (2 * Xb * T.K1 + 4 * Xb**2 * T.K2)
                     / (T.MPl2 * T.H**2)) == 0
    assert sp.cancel(rho_plus_p - 2 * Xb * T.K1) == 0


def test_row_targets_no_alphaB():
    """The 8 no-alpha_B rows' identity targets (S2 gate): exact expressions."""
    Hd = sp.Derivative(T.H, T.t)
    K = rows.cs2_times_aK_row()
    assert sp.cancel(K - (-2 * Hd / T.H**2 - T.rhom / (T.MPl2 * T.H**2))) == 0
    KT = rows.cs2_times_aK_row(aT_on=True)
    assert sp.cancel(sp.expand(KT - K + 2 * T.aT)) == 0        # -2 alpha_T offset
    KM = rows.cs2_times_aK_row(aM_on=True)
    aM = sp.Derivative(T.M2, T.t) / (T.H * T.M2)
    assert sp.cancel(sp.together(sp.expand(
        KM - (2 * aM - 2 * Hd / T.H**2 - T.rhom / (T.M2 * T.H**2))))) == 0
    KH = rows.cs2_times_aK_row(aH_on=True)
    aHd = sp.Derivative(T.aHh, T.t)
    exp_KH = (2 * T.aHh + 2 * aHd / T.H - 2 * (1 + T.aHh) * Hd / T.H**2
              - (1 + 2 * T.aHh) * T.rhom / (T.MPl2 * T.H**2))
    assert sp.cancel(sp.together(sp.expand(KH - exp_KH))) == 0
