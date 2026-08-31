"""P9(b) S2: P9b-0 sign-lemma machine checks, exact identity targets for all
8 no-alpha_B ladder rows, row-monotonicity substitutions, and the float
near-witnesses for the 7 (E)-candidate rows on the fiducial F*.

Conventions: FORMULATION.md v1.1 (exact dust DE eigen-speed pin for alpha_H
rows).  All "target" expressions are alpha_K * c_s^2 with c_s^2 the derived
exact dust eigen-speed and D = alpha_K (alpha_B = 0)."""
import pytest

sp = pytest.importorskip("sympy")

from derivation import rows
from derivation import tools as T
from p9b import fiducial

t = T.t
Hd = sp.Derivative(T.H, t)
aM_expr = sp.Derivative(T.M2, t) / (T.H * T.M2)
aHd = sp.Derivative(T.aHh, t)

FLAGS = {"K": (0, 0, 0), "KT": (1, 0, 0), "KM": (0, 1, 0), "KH": (0, 0, 1),
         "KTM": (1, 1, 0), "KTH": (1, 0, 1), "KMH": (0, 1, 1),
         "KTMH": (1, 1, 1)}


def derived_target(row):
    aT_on, aM_on, aH_on = (bool(v) for v in FLAGS[row])
    return rows.cs2_times_aK_row(aT_on=aT_on, aM_on=aM_on, aH_on=aH_on)


def master_instance(row):
    """fiducial.MASTER_TARGET written back in derivation variables."""
    aT_on, aM_on, aH_on = (bool(v) for v in FLAGS[row])
    M2row = T.M2 if aM_on else T.MPl2
    sub = {fiducial.s_aT: T.aT if aT_on else 0,
           fiducial.s_aM: aM_expr if aM_on else 0,
           fiducial.s_aH: T.aHh if aH_on else 0,
           fiducial.s_dH: aHd / T.H if aH_on else 0,
           fiducial.s_hd: Hd / T.H**2,
           fiducial.s_rm: T.rhom / (M2row * T.H**2)}
    return fiducial.MASTER_TARGET.subs(sub)


@pytest.mark.parametrize("row", list(FLAGS))
def test_master_target_matches_derived_rows(row):
    """The single master expression used by ALL numerics equals the derived
    alpha_K c_s^2 (exact dust eigen-speed) for every no-alpha_B row."""
    d = sp.cancel(sp.together(sp.expand(derived_target(row)
                                        - master_instance(row))))
    assert d == 0


def test_combined_row_offsets_exact():
    """Combined rows = pure-row offsets, exactly (S1 note §6 'sums')."""
    K, KT = derived_target("K"), derived_target("KT")
    KM, KH = derived_target("KM"), derived_target("KH")
    # {K,T,M} = {K,M} - 2 alpha_T ;  {K,T,H} = {K,H} - 2 alpha_T
    assert sp.cancel(sp.expand(derived_target("KTM") - KM + 2 * T.aT)) == 0
    assert sp.cancel(sp.expand(derived_target("KTH") - KH + 2 * T.aT)) == 0
    assert sp.cancel(sp.expand(derived_target("KTMH")
                               - derived_target("KMH") + 2 * T.aT)) == 0
    # {K,T} = {K} - 2 alpha_T (re-pin) and {K,M,H} vs its closed form
    assert sp.cancel(sp.expand(KT - K + 2 * T.aT)) == 0
    exp_KMH = (2 * T.aHh + 2 * aHd / T.H + 2 * (1 + T.aHh) * (aM_expr - Hd / T.H**2)
               - (1 + 2 * T.aHh) * T.rhom / (T.M2 * T.H**2))
    assert sp.cancel(sp.together(sp.expand(derived_target("KMH") - exp_KMH))) == 0


@pytest.mark.parametrize("sup,kill,sub", [
    ("KTH", "H", "KT"), ("KMH", "M", "KH"), ("KTMH", "T", "KMH"),
    ("KTM", "M", "KT"), ("KT", "T", "K")])
def test_row_monotonicity_substitution(sup, kill, sub):
    """Setting an extra alpha == 0 in a superset row's target gives the subset
    row's target exactly (the inheritance lemma's algebraic content; 0 is an
    admissible bounded free function for every alpha_i, i != K)."""
    e = derived_target(sup)
    if kill == "T":
        e = T.subs_fun(e, {T.aT: sp.Integer(0)})
    elif kill == "H":
        e = T.subs_fun(e, {T.aHh: sp.Integer(0)})
    else:                                   # alpha_M off: M^2 frozen at M_Pl^2
        e = T.subs_fun(e, {T.M2: T.MPl2})
    assert sp.cancel(sp.together(sp.expand(e - derived_target(sub)))) == 0


# --- P9b-0: sign-lemma machine checks ----------------------------------------------
def test_p9b0_hypotheses_inconsistent_assumption_engine():
    """I1 + D1 + D2 formally inconsistent: with alpha_K > 0 (D2(i)/D3),
    c_s^2 >= 0 (D2(ii)) and fhat < 0 (D1 phantom side, H^2 M_Pl^2 > 0),
    the identity alpha_K c_s^2 = fhat is False by sign alone."""
    aK = sp.Symbol("alpha_K", positive=True)
    cs2 = sp.Symbol("c_s2", nonnegative=True)
    fh = sp.Symbol("fhat", negative=True)
    assert sp.Eq(aK * cs2, fh) is sp.false          # sympy proves infeasibility
    assert (aK * cs2 - fh).is_positive is True      # obstruction strictly > 0


def test_p9b0_sign_decomposition_certificate():
    """Certificate: alpha_K c_s^2 - fhat = [alpha_K c_s^2] + [-fhat], a sum of
    a nonnegative term (pos * nonneg) and a strictly positive term, hence > 0;
    but I1 says it is 0.  Each piece's sign is machine-verified."""
    aK = sp.Symbol("alpha_K", positive=True)
    cs2 = sp.Symbol("c_s2", nonnegative=True)
    fh = sp.Symbol("fhat", negative=True)
    g = aK * cs2 - fh
    assert sp.expand(g - (aK * cs2) - (-fh)) == 0   # exact decomposition
    assert (aK * cs2).is_nonnegative is True
    assert (-fh).is_positive is True


def test_p9b0_phantom_side_forces_negative_cs2_fractions():
    """Exact Fraction check: on the phantom side (fhat < 0) the unique c_s^2
    solving I1 is fhat/alpha_K < 0 for every 0 < alpha_K <= 10^3, including
    the boundary alpha_K = 10^3 (the alpha_K -> infinity escape is closed by
    the frozen bound; even unbounded alpha_K only sends c_s^2 -> 0^-)."""
    from fractions import Fraction as F
    for aK in (F(1, 10**6), F(1), F(999), F(1000)):
        for fh in (F(-1, 10**9), F(-3, 7), F(-20)):
            cs2 = fh / aK
            assert cs2 < 0
            assert aK * cs2 == fh                    # I1 holds exactly


# --- near-witnesses on F* (float level; S3 owns Arb certification) -----------------
EXPECTED_MIN_CS2 = {"KT": 0.0761, "KM": 0.3036, "KH": 0.2803, "KH2": 1.8761,
                    "KTM": 0.5036, "KTH": 0.4803, "KMH": 1.3383, "KTMH": 1.5383}


@pytest.mark.parametrize("row", list(fiducial.WITNESSES))
def test_near_witness(row):
    r = fiducial.evaluate(row)
    # D1 on F*: crossing interior to J, f changes sign, transversal, rho_DE > 0
    assert r["crossing_interior"]
    assert r["fhat_left"] < 0 < r["fhat_right"]
    assert r["dfhat_dx_at_crossing"] > 0.1
    assert r["min_Ode"] > 0
    # D2 with the exact eigen-speed pin: c_s^2 = target/alpha_K > 0; D > 0
    assert r["min_cs2"] > 0.02
    assert r["D"] > 0
    assert r["min_cT2"] >= 0                        # tensor sector
    assert r["min_m2"] > 0                          # M*^2 > 0
    # frozen D3 bounds
    assert 0 < r["alpha_K"] <= 1e3
    assert r["max_abs_aT"] <= 10 and r["max_abs_aM"] <= 10
    assert r["max_abs_aH"] <= 10
    # loose pin of the recorded margin (regression guard)
    assert abs(r["min_cs2"] - EXPECTED_MIN_CS2[row]) < 0.02
