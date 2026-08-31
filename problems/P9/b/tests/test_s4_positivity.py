"""P9(b) S4: P9b-2 POS (Melville-Noller) certificates.

Layers: (1) exact sympy backing for the KT row-level impossibility (the row
identity used by KT_IMPOSSIBILITY); (2) the frozen applicability map;
(3) exact rational POS arithmetic (RHS of MN (10), side condition);
(4) stored-certificate JSON validation for the 7 applicable rows;
(5) coarse (64-bit) Arb re-run of the full S4 certification path;
(6) independent float grid re-checks of frozen verdicts and alternative
witnesses (numpy arithmetic, no Arb on this path)."""
import json
from fractions import Fraction as F
from pathlib import Path

import numpy as np
import pytest

sp = pytest.importorskip("sympy")
flint = pytest.importorskip("flint")

from p9b import certify, positivity as P
from p9b import witnesses as W

CERT_DIR = Path(__file__).resolve().parents[1] / "certificates"
APPLICABLE = list(P.POS_APPLICABLE)
ALT_ROWS = sorted(P.ALT)
NA_ROWS = sorted(P.POS_NOT_APPLICABLE)


def cert(row):
    return json.loads((CERT_DIR / f"P9b-2.{row}.json").read_text())


# --- (1) exact algebra behind the KT impossibility --------------------------------
def test_kt_row_identity_exact():
    """Row {K,T}: MASTER_TARGET_B reduces to fhat - 2 alpha_T, with the
    family-B identity fhat = -2 hd - rm (dust, M^2 = MPl^2)."""
    m = W.MASTER_TARGET_B.subs({W.s_aB: 0, W.s_dB: 0, W.s_aM: 0,
                                W.s_aH: 0, W.s_dH: 0})
    assert sp.expand(m - (-2 * W.s_hd - W.s_rm - 2 * W.s_aT)) == 0
    Om, Ode, w1 = sp.symbols("Om Ode w1")
    hd = -sp.Rational(3, 2) * (Om + w1 * Ode)   # family B, dust
    rm, fhat = 3 * Om, 3 * Ode * w1             # m2 = 1 on this row
    assert sp.expand((-2 * hd - rm) - fhat) == 0


def test_kt_certificate_records_impossibility():
    c = cert("KT")
    assert c["alternative"]["found"] is False
    imp = c["alternative"]["impossible"]
    assert "EMPTY" in imp["statement"] and len(imp["argument"]) == 4
    assert c["status"].startswith("(ii+)")


# --- (2) applicability map ---------------------------------------------------------
def test_applicability_partition():
    assert set(APPLICABLE) == {"KB", "KT", "KM", "KBT", "KBM", "KTM", "KBTM"}
    assert set(NA_ROWS) == {r for r in W.REG if "H" in r} and len(NA_ROWS) == 8
    assert set(APPLICABLE) | set(NA_ROWS) == set(W.REG)


@pytest.mark.parametrize("row", NA_ROWS)
def test_na_rows_have_no_certificate(row):
    assert not (CERT_DIR / f"P9b-2.{row}.json").exists()
    assert "beyond-Horndeski" in P.POS_NOT_APPLICABLE[row]
    with pytest.raises(AssertionError):
        P.certify_pos_row(row)


# --- (3) exact POS arithmetic ------------------------------------------------------
def test_pos_rhs_exact():
    assert P.pos_rhs(F(0)) == 0
    assert P.pos_rhs(F(-1, 10)) == F(-2, 9)
    assert P.pos_rhs(F(1, 10)) == F(2, 11)
    with pytest.raises(AssertionError):
        P.pos_rhs(F(-1))            # side condition 1 + alpha_T > 0


# --- (4) stored certificates -------------------------------------------------------
EXPECT_VERDICT = {"KM": "satisfies_exact"}


@pytest.mark.parametrize("row", APPLICABLE)
def test_certificate_json_valid(row):
    c = cert(row)
    assert c["claim"] == f"P9b-2.{row}" and c["conditional_on"] == "POS"
    assert c["prec_bits"] == certify.PREC_BITS
    assert "alpha_B^BS <= 2 alpha_T / (1 + alpha_T)" in c["pos_pin"]["inequality"]
    assert "1904.05874" in c["pos_pin"]["source"]
    assert c["frozen_witness"] == W.witness_desc(row)
    assert c["verdict"] == EXPECT_VERDICT.get(row, "violates_everywhere")
    assert c["pos_check"]["pass"] is True
    w = W.REG[row]
    if c["verdict"] == "violates_everywhere":
        assert c["violating_subinterval"] == [str(w["J"][0]), str(w["J"][1])]
        assert c["pos_check"]["max_upper_endpoint"] < 0
        assert c["pos_check"]["interval"] == c["violating_subinterval"]
    else:
        assert "0 <= 0" in c["pos_check"]["exact"]
    aT = w["aT"] if w["aT"] is not None else F(0)
    assert c["pos_rhs"].endswith(str(P.pos_rhs(aT)))


@pytest.mark.parametrize("row", ALT_ROWS)
def test_certificate_alternative_valid(row):
    a = cert(row)["alternative"]
    assert a["found"] is True
    s = P.ALT[row]
    assert a["witness"] == P.alt_desc(s)
    ck = a["checks"]
    assert all(v["pass"] for v in ck.values())
    assert ck["cs2_ge_cmin"]["bound"] == str(s["c_min"])
    assert ck["cs2_ge_cmin"]["min_lower_endpoint"] > float(s["c_min"]) - 1e-15
    assert ck["pos_strict"]["min_lower_endpoint"] > 0
    assert ck["fhat_neg_left"]["max_upper_endpoint"] < 0
    assert ck["fhat_pos_right"]["min_lower_endpoint"] > 0
    if s["b1"] is not None:
        assert ck["aB_ne_2"]["max_upper_endpoint"] < 2
        assert ck["D_ge_Dmin"]["min_lower_endpoint"] > float(s["D_min"]) - 1e-15
    if s["aM"] is not None:
        assert ck["M2_pos"]["min_lower_endpoint"] > 0


def test_km_certificate_saturation():
    c = cert("KM")
    assert "alternative" not in c
    assert c["status"].startswith("(i)") and "saturation" in c["status"]


# --- (5) coarse Arb re-run ---------------------------------------------------------
@pytest.mark.parametrize("row", APPLICABLE)
def test_recertify_coarse(row):
    prec0 = flint.ctx.prec
    flint.ctx.prec = 64
    try:
        c = P.certify_pos_row(row)
    finally:
        flint.ctx.prec = prec0
    stored = cert(row)
    assert c["verdict"] == stored["verdict"]
    assert c["status"] == stored["status"]
    if row in P.ALT:
        assert all(v["pass"] for v in c["alternative"]["checks"].values())


# --- (6) independent float re-checks -----------------------------------------------
def float_eval(spec, n=2001):
    """numpy re-implementation (independent of Arb) of ext_eval on a grid."""
    from p9b import fiducial as fid
    a1, a2 = float(spec["J"][0]), float(spec["J"][1])
    a = np.linspace(a1, a2, n)
    aT = float(spec["aT"]) if spec["aT"] is not None else 0.0
    aM = float(spec["aM"]) if spec["aM"] is not None else 0.0
    b = (float(spec["b1"]) * a + float(spec["b0"])
         if spec["b1"] is not None else np.zeros_like(a))
    db = float(spec["b1"]) * a if spec["b1"] is not None else np.zeros_like(a)
    g = 1.0 - b / 2.0
    rm = 3.0 * fid.Om(a) / a**aM
    hd = fid.hdot_over_H2(a)
    tgt = (-2 * g**2 * (1 + aT) + 2 * g * (1 + aM - hd) + db - rm)
    D = 1.0 + 1.5 * b**2
    q = float(P.pos_rhs(spec["aT"] if spec["aT"] is not None else F(0))) - b
    return dict(a=a, cs2=tgt / D, D=D, q=q, b=b)


@pytest.mark.parametrize("row", APPLICABLE)
def test_frozen_verdict_float(row):
    e = float_eval(P.spec_frozen(row))
    c = cert(row)
    if c["verdict"] == "violates_everywhere":
        assert e["q"].max() < 0
        assert e["q"].max() <= c["pos_check"]["max_upper_endpoint"] + 1e-12
    else:
        assert np.all(e["q"] == 0)


@pytest.mark.parametrize("row", ALT_ROWS)
def test_alternative_float_margins(row):
    s = P.ALT[row]
    e = float_eval(s)
    assert e["cs2"].min() > float(s["c_min"])
    assert e["q"].min() > 0                       # strict POS on the grid
    assert e["D"].min() >= float(s["D_min"]) - 1e-12
    assert np.abs(e["b"]).max() <= 10 and np.all(e["b"] < 2)
    assert np.any(e["b"] != 0) if s["b1"] is not None else np.all(e["b"] == 0)
    aT = s["aT"] if s["aT"] is not None else F(0)
    assert 1 + aT > 0 and abs(aT) <= 10
    ck = cert(row)["alternative"]["checks"]
    assert ck["cs2_ge_cmin"]["min_lower_endpoint"] <= e["cs2"].min() + 1e-12
    assert ck["pos_strict"]["min_lower_endpoint"] <= e["q"].min() + 1e-12


@pytest.mark.parametrize("row", APPLICABLE)
def test_ext_eval_matches_certify_row_eval(row):
    """On frozen witnesses the S4 evaluator agrees (ball overlap) with the
    S3 certify.row_eval on a sample subinterval."""
    lo, hi = F(7, 10), F(71, 100)
    a = P.ext_eval(P.spec_frozen(row), lo, hi)
    b = certify.row_eval(row, lo, hi)
    for k in ("cs2", "D", "E2", "Ode", "m2", "fhat", "aB"):
        assert a[k].overlaps(b[k])
