"""S2 part-(i) probe, route M3 (certificates/m3/*.json; criteria in notes/s2-preregistration.md).

Fast (~5 s): geometry certificate A0 recomputed; l = 50, 100 windows replayed with the stored
ground states (A1); stored A2/A3/overlap flags and the pre-registered fit.
P2_SLOW=1 (~3 min): l = 20 replay and the A2 negative control re-run at l = 100.
"""
import json
import os
from fractions import Fraction
from pathlib import Path

import pytest

from p2 import m3
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
