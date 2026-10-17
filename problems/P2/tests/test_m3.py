"""S2 part-(i) probe, route M3 (certificates/m3/*.json; criteria in notes/s2-preregistration.md).

Fast (~5 s): geometry certificate A0 recomputed; l = 50, 100 windows replayed with the stored
ground states (A1); the l = 100 P2-2 certificate replayed, with a control just below a_100;
stored A2/A2'/A3/overlap flags and the pre-registered fit.
P2_SLOW=1: l = 20 replay, the A2 negative control re-run at l = 100, and the A2' near-mode
control (S2 review, section 5) re-run at l = 50.
"""
import json
import os
from fractions import Fraction
from pathlib import Path

import pytest

from flint import arb, fmpq

from p2 import lg_bounds, m3
from p2.arbutil import precision
from p2.params import Point

CERT = Path(m3.CERT_DIR)
SLOW = pytest.mark.skipif(os.environ.get("P2_SLOW") is None, reason="P2_SLOW=1")


@pytest.fixture(scope="module")
def geo():
    pt = Point("PB", 256)
    with precision(m3.PREC):
        g, rays, tails = m3.geometry_certificate(pt)
    return pt, g, rays, tails


def _replay(geo, l):
    pt, _, rays, tails = geo
    st = json.loads((CERT / f"m3_PB_l{l}.json").read_text())
    nodes = [(Fraction(x), Fraction(y)) for x, y in st["ground_state_nodes"]]
    with precision(m3.PREC):
        c = m3.certify_l(pt, l, rays, tails, nodes=nodes, lam_float=st["float_lambda_min"])
    return st, c


def test_geometry_A0(geo):
    _, g, _, _ = geo
    assert g["A0_pass"] and g["identity_Re_rstar_holds"] and g["lower_left_sector_free"]
    assert all(v["ok"] for v in g["upper_right_special_segments"].values())


@pytest.mark.parametrize("l", [50, 100])
def test_window_replay(geo, l):
    st, c = _replay(geo, l)
    assert c["A1_pass"] and st["A1_pass"]
    assert c["mu_cert"] == st["mu_cert"]
    assert c["c_nonneg"] and c["jost_decay_condition"]


@SLOW
def test_window_replay_l20(geo):
    st, c = _replay(geo, 20)
    assert c["A1_pass"] and c["mu_cert"] == st["mu_cert"]


@pytest.mark.parametrize("l", [20, 50, 100])
def test_stored_controls_and_overlap(l):
    st = json.loads((CERT / f"m3_PB_l{l}.json").read_text())
    assert st["A2_control"]["A2_pass"]                 # the PS-mode box is NOT certified
    assert st["A3_mu_cert_le_float"]
    assert st["P2_2_overlap"]["a_l_le_top"] and st["P2_2_overlap"]["lg_criterion_certified"]


def test_lg_l100_replay():
    """The l = 100 P2-2 certificate (S2; S1 code): the criterion replays at the stored exact
    (a_100, r_X); a_100 is the a_l of the l = 100 window and lies below its top. Control: 0.2%
    below a_100 the criterion is certainly above 2, so the stored threshold is not vacuous."""
    st = json.loads((CERT / "lg_PB_l100.json").read_text())
    a0 = fmpq(*[int(x) for x in st["W_l_exact"].split("/")])
    rX = fmpq(*[int(x) for x in st["matching_point_exact"].split("/")])
    pt = Point("PB")
    lg = lg_bounds.LG(pt, 100)
    res = lg.criterion(a0, rX)
    assert res is not None and res["ok"] and res["value"] < 2
    win = json.loads((CERT / "m3_PB_l100.json").read_text())
    assert win["P2_2_overlap"]["a_l_exact"] == st["W_l_exact"]
    assert arb(a0) < (arb(201) / 2) * pt.OmegaP + 6 * pt.gamma0      # a_100 < window top
    low = lg.criterion(a0 * fmpq(499, 500), rX)
    assert low is not None and low["value"] > 2


@pytest.mark.parametrize("l", [20, 50, 100])
def test_stored_near_mode_control(l):
    """A2' (stored): the boxes IN and STRAD contain the l-th PS mode and are not certified;
    OUT excludes it. The re-run is test_near_mode_control_rerun_l50."""
    st = json.loads((CERT / f"m3_PB_l{l}.json").read_text())
    nm = st["A2prime_near_mode"]
    assert nm["A2prime_pass"]
    km = float(Point("PB").km.mid())
    re_ps, im_ps = (float(x) * km for x in nm["PS_mode_over_kappa_minus"])
    for name in ("IN", "STRAD"):
        a1, a2, b1, b2 = nm[name]["box"]
        assert a1 < re_ps < a2 and b1 < im_ps < b2
        assert not nm[name]["certified"] and nm[name]["mu_over_muHO"] < -0.05
    b1, b2 = nm["OUT"]["box"][2:]
    assert not b1 <= im_ps <= b2
    assert st["A2_control"]["b_box_design"]["checker_run"] is False


def test_summary_fit_and_verdict():
    s = json.loads((CERT / "m3_PB_summary.json").read_text())
    assert s["verdict"].startswith("GO")
    assert s["fit"]["C_over_l"]["projected_l0"] <= 100


@SLOW
def test_negative_control_rerun_l100(geo):
    pt, _, rays, tails = geo
    st = json.loads((CERT / "m3_PB_l100.json").read_text())
    nodes = [(Fraction(x), Fraction(y)) for x, y in st["ground_state_nodes"]]
    with precision(m3.PREC):
        out = m3.control_l(pt, 100, rays, tails, nodes, st["float_lambda_min"])
    assert out["A2_pass"]


@SLOW
def test_near_mode_control_rerun_l50(geo):
    """A2' replay at l = 50: the checker, run with the stored near-mode ground state on the stored
    boxes, reproduces every verdict and mu; the boxes containing the mode are not certified.
    (The float design step is not bit-reproducible, so the replay uses the stored nodes, as A4
    does for the windows.)"""
    pt, _, rays, tails = geo
    l = 50
    st = json.loads((CERT / f"m3_PB_l{l}.json").read_text())["A2prime_near_mode"]
    nodes = [(Fraction(x), Fraction(y)) for x, y in st["ground_state_nodes"]]
    with precision(m3.PREC):
        muHO = float((2 * (l + arb(1) / 2) * pt.OmegaP * pt.gamma0).mid())
        lam = st["float_lambda_at_design_over_muHO"] * muHO       # the bisection target only
        for name in ("IN", "STRAD", "OUT"):
            box = tuple(arb(v) for v in st[name]["box"])
            c = m3.certify_l(pt, l, rays, tails, nodes=nodes, lam_float=lam, box=box)
            assert c["A1_pass"] == st[name]["certified"]
            assert c["mu_cert"] == st[name]["mu_cert"]
    assert not st["IN"]["certified"] and not st["STRAD"]["certified"]
