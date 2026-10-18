"""S4 part (i) at P_B uniformly in l >= 20: regimes T, A (part1_T), S (throat), L (lowfreq), the bulk facts,
the cover and the controls (certificates/part1/; criteria in notes/s4-preregistration.md).

Fast (~30 s): the sympy identities and bulk facts; T's exact splitting against S2's E and its derivative
formulas against finite differences; one T zeta-box, the A box, one L omega-box, both S closing boxes and
one S zeta-box per side, recomputed; affine-form enclosures against sampled points; K2 at one point; the
stored certificates pass and the cover file is consistent.
P2_SLOW=1: one stored certificate box per regime replayed exactly, N1 at l = 20, N3, N5.
"""
import json
import os
import random
from fractions import Fraction

import pytest
from flint import acb, arb

from p2 import part1, part1_T, throat
from p2 import lowfreq as LF
from p2.affine import AF
from p2.arbutil import precision
from p2.m3 import Checker as S2Checker
from p2.m3 import fa
from p2.params import Point

CERT = part1.CERT_DIR
SLOW = pytest.mark.skipif(os.environ.get("P2_SLOW") is None, reason="P2_SLOW=1")
W_BOUND = Fraction(99, 12500)


@pytest.fixture(scope="module")
def pt():
    return Point("PB", 256)


@pytest.fixture(scope="module")
def geo(pt):
    with precision(128):
        return part1_T.setup(pt)


def test_bulk(pt):
    with precision(128):
        b = part1.bulk(pt)
    assert b["pass"] and b["identity_sympy_ok"] and b["splitting_V_eq_U_over_h2_plus_Wt_sympy_ok"]
    assert Fraction(b["W_bound"]) == W_BOUND


def test_T_splitting_and_derivatives(pt, geo):
    _, rays, tails = geo
    with precision(128):
        ch = part1_T.TChecker(pt, rays, tails, part1_T.DESIGN_T)
        l = 37
        h = arb(1) / (l + arb(1) / 2)
        w = acb(arb("0.93") * pt.OmegaP / h, arb("-0.017"))
        zeta, nu = w.real * h, w.imag
        mc = S2Checker(pt, l, rays, {"s1": 7.0, "p": 4, "phi_inf": -1.4375, "eta": 0.5})
        for sv in ("-3.3", "-1.1", "0.2", "2.2", "3.9", "25.3"):
            s = arb(sv)
            side = 1 if s > 0 else -1
            a = abs(float(s.mid()))
            zone = 0 if a <= 2 else (1 if a <= 4 else 2)
            r = rays[side].at(s)
            p = ch.parts(s, r, zeta, zeta, nu, -side, zone, side)
            G, Gp, _ = ch.D.gamma(s, zone, side)
            Es = (p["E0"] + h * p["E1"] + h * h * p["E2"]) / (h * h)
            E2, _, _ = mc.E_parts(s, r, w * w, G / h, Gp / h, -side)
            assert abs(Es - E2) < 1e-20 * abs(E2)
            e = arb("1e-7")
            pp = ch.parts(s + e, rays[side].at(s + e), zeta, zeta, nu, -side, zone, side)
            pm = ch.parts(s - e, rays[side].at(s - e), zeta, zeta, nu, -side, zone, side)
            for k, kd in (("E0", "E0p"), ("E1", "E1p"), ("E2", "E2p"), ("E0p", "E0pp")):
                assert abs((pp[k] - pm[k]) / (2 * e) - p[kd]) < 1e-9


def test_T_box_and_A_box(pt, geo):
    _, rays, tails = geo
    with precision(128):
        for design, zb in ((part1_T.DESIGN_T, (Fraction(4, 5), Fraction(17, 20))),
                           (part1_T.DESIGN_A, (Fraction(3, 2), Fraction(3, 2)))):
            ch = part1_T.TChecker(pt, rays, tails, design)
            r = part1_T.certify_box(ch, pt, zb, design)
            assert r["pass"] and r["n_fails"] == 0
        assert part1_T.monotonicity_A(pt, part1_T.DESIGN_A)["pass"]


def test_L_box(pt):
    with precision(128):
        ch = LF.LChecker(pt, W_BOUND)
        r = ch.box((Fraction(159, 20), Fraction(8), Fraction(-11, 20), Fraction(-1, 2)))
    assert r["pass"]


def test_S_closing_and_one_box_each_side(pt):
    g = throat.FGeo(pt)
    with precision(128):
        for side in (1, -1):
            assert throat.certify_closing(pt, g, side, W_BOUND)["pass"]
        assert throat.certify_zbox(pt, g, -1, Fraction(3, 10), Fraction(8, 25), W_BOUND)["pass"]
        assert throat.certify_zbox(pt, g, 1, Fraction(22, 25), Fraction(9, 10), W_BOUND)["pass"]


def test_affine_forms_enclose():
    rng = random.Random(1)
    with precision(128):
        Rw, Rh = arb("0.01"), arb("0.002")
        x = AF(acb("0.7", "-0.2"), acb("1.3", "0.4"), acb("-2", "1"), "1e-6", Rw, Rh)
        y = AF(acb("-0.3", "0.9"), acb("0.2", "-0.5"), acb("3", "0"), "2e-6", Rw, Rh)
        for f, F in ((lambda a, b: a * b, lambda a, b: a * b), (lambda a, b: a / b, lambda a, b: a / b),
                     (lambda a, b: a + b * b, lambda a, b: a + b * b)):
            z = f(x, y)
            for _ in range(200):
                dw = complex(rng.uniform(-1, 1), rng.uniform(-1, 1)) * 0.0099 / 1.42
                dh = rng.uniform(-1, 1) * 0.00199
                xv = complex(0.7, -0.2) + complex(1.3, 0.4) * dw + complex(-2, 1) * dh
                yv = complex(-0.3, 0.9) + complex(0.2, -0.5) * dw + complex(3, 0) * dh
                zv = F(xv, yv)
                assert abs(zv - complex(float(z.c0.real.mid()), float(z.c0.imag.mid()))) <= float(z.dev().upper()) * 1.0001


def test_K2_one_point():
    r = part1.check_K2(ls=(20,), points=[(5, -0.55)])
    assert r["pass"], r


@pytest.mark.skipif(not (CERT / "cover_PB.json").exists(), reason="no assembled certificates yet")
def test_certificates_pass_and_cover():
    for n in ("bulk_PB", "T_PB", "A_PB", "S_PB_plus", "S_PB_c", "L_PB_plus", "L_PB_c", "cover_PB"):
        st = json.loads((CERT / f"{n}.json").read_text())
        assert st["pass"], n
    cv = part1.cover()
    st = json.loads((CERT / "cover_PB.json").read_text())
    assert cv["pass"] and cv["tiling"] == st["tiling"] and cv["inequalities"] == st["inequalities"]


@SLOW
@pytest.mark.skipif(not (CERT / "cover_PB.json").exists(), reason="no assembled certificates yet")
def test_replay_one_box_per_regime(pt, geo):
    _, rays, tails = geo
    g = throat.FGeo(pt)
    with precision(128):
        st = json.loads((CERT / "T_PB.json").read_text())
        bx = st["boxes"][0]
        ch = part1_T.TChecker(pt, rays, tails, part1_T.DESIGN_T)
        r = part1_T.certify_box(ch, pt, tuple(Fraction(x) for x in bx["zeta_over_OmegaP"]), part1_T.DESIGN_T,
                                leaves=bx["leaves"])
        assert r["mu_cert"] == bx["mu_cert"] and r["pass"]
        for side, tag in ((1, "plus"), (-1, "c")):
            st = json.loads((CERT / f"S_PB_{tag}.json").read_text())
            x = st["leaves"][len(st["leaves"]) // 2]
            r = throat.certify_zbox(pt, g, side, *(Fraction(v) for v in x["zeta_over_OmegaP"]), W_BOUND,
                                    design={"design": x["design"], "X_a": x["X_a"], "k": x["k"], "target": x["target"]})
            assert r["pass"] and [e["mu_cells"] for e in r["eps_boxes"]] == [e["mu_cells"] for e in x["eps_boxes"]]
        st = json.loads((CERT / "L_PB_plus.json").read_text())
        x = st["leaves_detail"][len(st["leaves_detail"]) // 3]
        r = LF.LChecker(pt, W_BOUND).box(tuple(Fraction(v) for v in x["box"]))
        assert r["pass"] and r["plus"]["Re_y_bound"] == x["Re_y_bound"]


@SLOW
def test_controls_N1_l20_N3_N5(pt):
    with precision(128):
        n12 = part1.control_N1_N2(pt)
        assert n12["N1"]["20"]["pass"] and n12["N2"]["pass"]
        n3 = part1.control_N3(pt, W_BOUND)
        assert all(r["pass"] for r in n3.values())
        assert all(r["OUT_half_depth"]["certified"] for r in n3.values())
        n5 = part1.control_N5(pt, W_BOUND)
        assert all(r["pass"] for r in n5.values())
