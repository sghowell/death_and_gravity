"""P9(b) S3: Arb-certified existence witnesses for the 15 rows S != {alpha_K}.

Layers: (1) the braiding-extended master target MASTER_TARGET_B (the single
expression driving every S3 certificate) equals the DERIVED exact dust DE
eigen-speed times D for all alpha_B rows, exactly (sympy); (2) the
inheritance-lemma algebra extended to the braiding rows; (3) coarse Arb
re-certification of every row (the real certification path at lower
precision); (4) validation of the stored certificate JSONs; (5) float
regression of the frozen witnesses against the certified bounds."""
import json
from fractions import Fraction as F
from pathlib import Path

import pytest

sp = pytest.importorskip("sympy")
flint = pytest.importorskip("flint")

from derivation import rows
from derivation import tools as T
from p9b import certify, fiducial
from p9b import witnesses as W

CERT_DIR = Path(__file__).resolve().parents[1] / "certificates"
ROWS = list(W.REG)
B_ROWS = [r for r in ROWS if "B" in r]

t = T.t
H = T.H


def derived_target(row):
    B_on, T_on, M_on, H_on = W.flags(row)
    return rows.cs2_times_aK_row(aB_on=B_on, aT_on=T_on, aM_on=M_on, aH_on=H_on)


def scalarize(e, aM_on):
    """Derivation variables -> the scalar variables of MASTER_TARGET_B."""
    M2row = T.M2 if aM_on else T.MPl2
    e = e.subs(sp.Derivative(T.aBg, t), -H * W.s_dB / 2)
    e = e.subs(sp.Derivative(T.aHh, t), H * W.s_dH)
    e = e.subs(sp.Derivative(T.M2, t), W.s_aM * H * T.M2)
    e = e.subs(sp.Derivative(H, t), W.s_hd * H**2)
    e = e.subs({T.aBg: -W.s_aB / 2, T.aHh: W.s_aH, T.aT: W.s_aT,
                T.rhom: W.s_rm * M2row * H**2})
    return sp.cancel(sp.together(sp.expand(e)))


def master_instance(row):
    B_on, T_on, M_on, H_on = W.flags(row)
    m = W.MASTER_TARGET_B
    if not B_on:
        m = m.subs({W.s_aB: 0, W.s_dB: 0})
    if not T_on:
        m = m.subs(W.s_aT, 0)
    if not M_on:
        m = m.subs(W.s_aM, 0)
    if not H_on:
        m = m.subs({W.s_aH: 0, W.s_dH: 0})
    return m


@pytest.mark.parametrize("row", B_ROWS)
def test_master_b_matches_derived_rows(row):
    """MASTER_TARGET_B == derived D c_s^2 (exact dust eigen-speed, rev-v1.1
    pin) for every braiding row, exactly.  (The 8 alpha_B = 0 instances are
    covered by S2's test_master_target_matches_derived_rows.)"""
    _, _, M_on, _ = W.flags(row)
    d = sp.cancel(sp.together(sp.expand(
        scalarize(derived_target(row), M_on) - master_instance(row))))
    assert d == 0


def test_master_b_reduces_to_s2_master():
    """At alpha_B = 0 the S3 master equals the S2 master exactly."""
    sub = {W.s_aB: 0, W.s_dB: 0, W.s_aT: fiducial.s_aT, W.s_aM: fiducial.s_aM,
           W.s_aH: fiducial.s_aH, W.s_dH: fiducial.s_dH,
           W.s_hd: fiducial.s_hd, W.s_rm: fiducial.s_rm}
    d = sp.expand(W.MASTER_TARGET_B.subs(sub) - fiducial.MASTER_TARGET)
    assert d == 0


@pytest.mark.parametrize("sup,kill,sub", [
    ("KB", "B", "K"), ("KBT", "T", "KB"), ("KBM", "M", "KB"),
    ("KBH", "H", "KB"), ("KBTM", "M", "KBT"), ("KBTH", "H", "KBT"),
    ("KBMH", "H", "KBM"), ("KBTMH", "T", "KBMH")])
def test_inheritance_lemma_b_rows(sup, kill, sub):
    """Inheritance-lemma algebra extended to braiding rows: setting an extra
    alpha == 0 (an admissible D3 function) in the superset row's derived
    target gives the subset row's target exactly.  With S2's five no-B cases
    this chains every row down to a two-operator row."""
    e = derived_target(sup)
    if kill == "B":
        e = T.subs_fun(e, {T.aBg: sp.Integer(0)})
    elif kill == "T":
        e = T.subs_fun(e, {T.aT: sp.Integer(0)})
    elif kill == "H":
        e = T.subs_fun(e, {T.aHh: sp.Integer(0)})
    else:
        e = T.subs_fun(e, {T.M2: T.MPl2})
    key = "K" if sub == "K" else sub
    tgt = rows.cs2_times_aK_row() if key == "K" else derived_target(sub)
    assert sp.cancel(sp.together(sp.expand(e - tgt))) == 0


# --- Arb re-certification (coarse) and stored-certificate validation ---------------
@pytest.mark.parametrize("row", ROWS)
def test_recertify_coarse(row):
    """Re-run the full Arb certification path per row at coarse precision
    (64 bits): every check must still certify strictly."""
    prec0 = flint.ctx.prec
    flint.ctx.prec = 64
    try:
        cert = certify.certify_row(row)
    finally:
        flint.ctx.prec = prec0
    assert cert["verdict"] == "CERTIFIED (E)"
    assert all(c["pass"] for c in cert["checks"].values())
    c = cert["checks"]["cs2_ge_cmin"]
    assert c["min_lower_endpoint"] > float(W.REG[row]["c_min"]) - 1e-15
    assert cert["checks"]["fhat_neg_left"]["max_upper_endpoint"] < 0
    assert cert["checks"]["fhat_pos_right"]["min_lower_endpoint"] > 0


@pytest.mark.parametrize("row", ROWS)
def test_certificate_json_valid(row):
    """The committed certificate matches the frozen registry and passes."""
    w = W.REG[row]
    cert = json.loads((CERT_DIR / f"P9b-1.{row}.json").read_text())
    assert cert["claim"] == f"P9b-1.{row}"
    assert cert["verdict"] == "CERTIFIED (E)"
    assert cert["prec_bits"] == certify.PREC_BITS
    assert cert["witness"] == W.witness_desc(row)
    assert all(c["pass"] for c in cert["checks"].values())
    c = cert["checks"]["cs2_ge_cmin"]
    assert c["bound"] == str(w["c_min"]) and c["relation"] == ">"
    assert c["min_lower_endpoint"] > float(w["c_min"]) - 1e-15
    assert c["n_subintervals"] >= 1 and c["max_depth"] <= certify.MAX_DEPTH
    assert c["interval"] == [str(w["J"][0]), str(w["J"][1])]
    if w["b1"] is not None:
        assert cert["checks"]["D_ge_Dmin"]["bound"] == str(w["D_min"])
        assert cert["checks"]["aB_ne_2"]["max_upper_endpoint"] < 2
    else:
        assert cert["checks"]["D_ge_Dmin"]["exact"].startswith("D = alpha_K")
    if w["aM"] is not None:
        assert cert["checks"]["M2_pos"]["pass"]
    assert cert["checks"]["fhat_neg_left"]["max_upper_endpoint"] < 0
    assert cert["checks"]["fhat_pos_right"]["min_lower_endpoint"] > 0
    assert cert["crossing"]["a_c"] == "3/4" and cert["crossing"]["ivt_gap"] == "1/50"


# --- float regression of the frozen witnesses --------------------------------------
@pytest.mark.parametrize("row", ROWS)
def test_witness_float_margins(row):
    """Float grid scan agrees with (and cannot beat) the rigorous bounds, and
    the witness sits inside every frozen D1/D2/D3 requirement."""
    w = W.REG[row]
    s = W.scan(row, n=1001)
    assert s["crossing_interior"]
    assert s["fhat_left"] < 0 < s["fhat_right"]
    assert s["min_Ode"] > 0
    assert s["min_cs2"] > float(w["c_min"])
    assert s["min_D"] >= float(w["D_min"]) - 1e-12
    assert s["min_cT2"] >= 0
    assert s["min_m2"] > 0
    assert s["max_abs_aB"] <= 0.5 < 2            # braiding bound, aB != 2
    # frozen D3 bounds (exact rationals in the registry)
    assert 0 < W.ALPHA_K <= 1000
    for key in ("aT", "aM", "aH", "b1"):
        if w[key] is not None:
            assert abs(w[key]) <= 10
    # the stored rigorous lower bound cannot exceed the float grid minimum
    cert = json.loads((CERT_DIR / f"P9b-1.{row}.json").read_text())
    assert cert["checks"]["cs2_ge_cmin"]["min_lower_endpoint"] <= s["min_cs2"] + 1e-12
