"""S4-4 tests: the near-sonic layer engine (tc_layer).

Fast by default (~20 s, dominated by the Atoms fixture): pivot-form identities
on synthetic Laurent polynomials, small/large-form agreement, box-geometry
containment, the float check of the September F values, the layer algebra of
Sigma against the independent direct chain at points, and one end-to-end
adaptive cover at a single xi.  The full-piece certifications (one Laurent
build is ~55 s) are gated behind P4_SLOW=1; the gap-step test also needs the
gap tube (P4_GAP_TUBE=<file>, or problems/P4/data/cache/tube_gap.json).
"""
import json
import glob
import os
import random

import pytest
from flint import acb, arb

import p4.validated.tc_energy as te
from p4.validated import tc_layer as tl
from p4.validated.arbseries import precision
from p4.validated.tc_tm import LC, TM

SLOW = pytest.mark.skipif(os.environ.get("P4_SLOW") is None,
                          reason="P4_SLOW=1 for the full-piece runs")
P4 = os.path.join(os.path.dirname(__file__), "..")
LAYER_PARTS = os.path.join(P4, "results", "theorem_c", "parts", "layer")
INF = float("inf")


@pytest.fixture(scope="module")
def atoms():
    return te.Atoms(full=True)


@pytest.fixture(scope="module")
def sonic_ex():
    return te.sonic_expansion_cached(320)


# ------------------------------------------------------------ pivot forms
def _synthetic(offs, seed=0, R=1e-3):
    rnd = random.Random(seed)
    kp = []
    for o in offs:
        c = [arb(rnd.uniform(0.5, 2.0) * rnd.choice((-1, 1)))] + \
            [arb(rnd.uniform(-3, 3)) for _ in range(3)]
        kp.append(LC(TM(c, R), o))
    return kp


def _direct(kp, l, k):
    acc = acb(0)
    for x in reversed(kp):
        acc = acc * k + x.tm.eval(l, l) * arb(l) ** x.off
    return acc


@pytest.mark.parametrize("offs", [[1] * 9 + [2],          # nB-like (layer)
                                  [0] * 9 + [1],          # nB'-like
                                  [3] * 8,                # dB-like (uniform)
                                  [0, 1, 1, 2, 3]])       # mixed
def test_pform_identity(offs):
    """P = ell^m z^p Q(z, w) = ell^m z^d Q^L(zeta, w) at points."""
    with precision(256):
        kp = _synthetic(offs)
        f = tl.PForm(kp, "syn")
        assert all(v >= 0 for v in f.e.values()) and all(v >= 0 for v in f.g.values())
        for l, kap in ((3e-4, acb(0.5, 400.0)), (1e-4, acb(2.0, 3e4)),
                       (7e-4, acb(0.0, 250.0))):
            P = _direct(kp, l, kap)
            z, w = arb(l) * kap, 1 / kap
            L = (l, l)
            Ps = arb(l) ** f.m * z ** f.p * f.small(L, z, w)
            Pl = arb(l) ** f.m * z ** f.d * f.large(L, 1 / z, w)
            for v in (Ps, Pl):
                assert abs(v - P).upper() <= 1e-30 * abs(P).upper() + 1e-60


def test_pform_pivot_structure():
    f = tl.PForm(_synthetic([1] * 9 + [2]), "nB")
    assert (f.d, f.p, f.m) == (9, 8, -7) and f.pivot_on_layer and f.top_on_layer
    g = tl.PForm(_synthetic([3] * 8), "dB")
    assert (g.d, g.p, g.m) == (7, 7, -4)
    assert all(v == 0 for v in g.g.values())


def test_ratio_small_large_agree():
    with precision(256):
        N = tl.PForm(_synthetic([0] * 9 + [1], seed=3), "N")
        D = tl.PForm(_synthetic([1] * 9 + [2], seed=4), "D")
        L = (2e-4, 2.1e-4)
        for zc in (acb(0.01, 4.0), acb(0.0, 7.0)):
            z = acb(zc.real + arb(0, 1e-6), zc.imag + arb(0, 1e-3))
            w = arb(3e-4) / z
            a = tl.ratio_hull(N, D, 1, L, None, z, None, w, False, True)
            b = tl.ratio_hull(N, D, 1, L, None, None, 1 / z, w, True, True)
            assert a.overlaps(b)


def test_pw_hull_refuses_negative_power_at_zero():
    with pytest.raises(ArithmeticError):
        tl.pw_hull((0.0, 1e-6), -1)
    assert tl.pw_hull((0.0, 1e-6), 2).contains(arb(0))


# ------------------------------------------------------------- box geometry
class _LD:
    laurent = True
    xi_of_ell = staticmethod(lambda l: l)
    ell_of_xi = staticmethod(lambda x: x)

    def xi_hull(self, L):
        return tl.iv(L[0], L[1])


@pytest.mark.parametrize("box", [(0.0, 1e-6, 0.0, 0.03, 0.0, 1e-3, 224.0, INF),
                                 (0.0, 1e-6, 1.0, 1.6, 0.5, 0.6, 224.0, INF),
                                 (0.0, 1e-6, 3.4, 4.6, 1e3, INF, 224.0, INF),
                                 (8e-5, 9.5e-5, 0.0, 6.0, 0.02, 0.05, 224.0, 500.0)])
def test_box_geometry_contains(box):
    rnd = random.Random(7)
    ld = _LD()
    with precision(128):
        g = tl.box_geometry(ld, box)
        assert g is not None
        la, lb, e0, e1, s0, s1, t0, t1 = box
        n = 0
        while n < 200:
            l = rnd.uniform(la, lb)
            tau = t0 * (min(t1, 1e12) / t0) ** rnd.random()
            s = l * tau
            if not (s0 <= s <= s1) or l <= 0:
                continue
            n += 1
            eps = rnd.uniform(e0, e1)
            k = acb(eps, tau)
            assert g["L"][0] <= l <= g["L"][1]
            assert g["w"].contains(1 / k)
            assert g["itau"].contains(arb(1) / arb(tau))
            assert g["lam"].contains(arb(l))
            if g["use_large"]:
                assert g["zeta"].contains(1 / (arb(l) * k))
            else:
                assert g["z"].contains(arb(l) * k)


# ------------------------------------------------- float check (design data)
def test_float_check_september(atoms, sonic_ex):
    """The September corner3 values (notes/s4-commutator.md section 8):
    S xi^2/tau^2 = F(s, eps) at xi = 1e-4."""
    pr = tl.point_pairs(atoms, sonic_ex, 1e-4)
    want = {(0.0, 0.03): -0.904, (0.0, 1.0): -3.21, (0.0, 1e4): -15.46,
            (3.0, 0.3): -40.9, (3.0, 1.0): -43.0, (3.0, 3.0): -49.8}
    for (eps, s), v in want.items():
        F = tl.point_F(atoms, sonic_ex, 1e-4, eps, s, pr=pr)
        assert abs(float(F.mid()) - v) < 0.01 * max(1.0, abs(v)), (eps, s, F)


@pytest.fixture(scope="module")
def point_ld(atoms, sonic_ex):
    """LayerData at the single point xi = 2e-3 (constant TMs, no Laurent
    structure): the layer algebra is then checked against the direct chain."""
    xi = 2e-3
    pr = tl.point_pairs(atoms, sonic_ex, xi)
    with precision(320):
        ld = tl.LayerData(pr, False, xik=arb(xi), label="point")
    return xi, pr, ld


@pytest.mark.parametrize("eps,tau", [(0.0, 224.0), (0.5, 1000.0), (3.0, 3e4),
                                     (6.0, 2e6), (0.01, 1e9)])
def test_sigma_layer_algebra_vs_direct(point_ld, eps, tau):
    xi, pr, ld = point_ld
    with precision(320):
        box = (0.0, 0.0, eps, eps, xi * tau, xi * tau, tau, tau)
        S_layer = tl.sigma_box(ld, box)
        S_dir, _ = tl.point_S_direct(pr, eps, tau)
        S_dir = S_dir * arb(xi) ** 2 / arb(tau) ** 2
        assert S_layer.overlaps(S_dir)
        assert float(S_layer.rad()) < 1e-6 * max(1.0, abs(float(S_layer.mid())))


@pytest.mark.parametrize("tau_min", [224.0, 14.0])
def test_certify_single_xi(point_ld, tau_min):
    """End-to-end adaptive cover at one xi: all eps in [0, 6], tau >= tau_min."""
    xi, pr, ld = point_ld
    res = tl.run_layer_piece(ld, 0.0, 0.0, tau_min=tau_min)
    assert res["ok"] and res["certified_bound"] < 0, res["fails"]
    assert res["params"]["tau_min"] == tau_min


# ----------------------------------------------------------- full pieces
@SLOW
def test_closing_piece_certifies(atoms, sonic_ex):
    ld = tl.sonic_layer_data(atoms, sonic_ex, 1e-6, label="sonic/s00")
    res = tl.run_layer_piece(ld, 0.0, 1e-6)
    assert res["ok"] and res["certified_bound"] <= -0.25, res["fails"]
    assert res["forms"]["nB"]["p"] == 8 and res["forms"]["nB"]["d"] == 9


@SLOW
def test_mid_sonic_piece_certifies(atoms, sonic_ex):
    lo, hi = 7.949684720339076e-05, 9.539621664406892e-05
    ld = tl.sonic_layer_data(atoms, sonic_ex, hi, label="sonic/s25")
    res = tl.run_layer_piece(ld, lo, hi)
    assert res["ok"] and res["certified_bound"] <= -0.5, res["fails"]


@SLOW
def test_gap_step_certifies(atoms):
    path = os.environ.get("P4_GAP_TUBE") or os.path.join(P4, "data", "cache", "tube_gap.json")
    if not os.path.exists(path):
        pytest.skip("gap tube cache not available")
    from p4.validated.lintube import Tube
    tube = Tube.load(path)
    steps = [sd for sd in tube.steps if float(sd.x) - sd.h >= -0.056 - 1e-12
             and float(sd.x) <= -0.008 + 1e-12]
    sd = steps[88]
    ld = tl.tube_layer_data(atoms, sd, label="gap/gap-088")
    res = tl.run_layer_piece(ld, 0.0, float(sd.h))
    assert res["ok"] and res["certified_bound"] <= -0.5, res["fails"]


# ---------------------------------------------------------------- artifacts
def test_layer_parts_consistency():
    files = sorted(glob.glob(os.path.join(LAYER_PARTS + "*", "*", "*.json")))
    if not files:
        pytest.skip("no layer parts yet")
    for f in files:
        d = json.load(open(f))
        assert d["schema"] == "s44-layer/1"
        tmin = 224.0 if os.path.basename(os.path.dirname(os.path.dirname(f))) == "layer" \
            else float(os.path.basename(os.path.dirname(os.path.dirname(f)))[len("layer_t"):])
        assert d["scope"]["tau"][0] == tmin and d["scope"]["eps"] == [0.0, 6.0]
        if d["status"] == "pass":
            assert d["result"]["params"]["tau_min"] == tmin
        if d["status"] == "pass":
            assert d["result"]["ok"] and d["certified_bound"] < 0
