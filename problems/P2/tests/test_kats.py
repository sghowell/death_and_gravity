"""KAT-1 ... KAT-5 and KAT-7 results (results/kat_float_*.json); criteria fixed in p2/kats.py.

Fast: stored outcomes + one live Table I mode (~5 s).  P2_SLOW=1 re-runs KAT-1 (~3 min).
"""
import json
import os
from pathlib import Path

import mpmath as mp
import pytest

from p2 import kats
from p2.qnm_float import Geom, qnm

RES = Path(kats.RESULTS)
SLOW = pytest.mark.skipif(os.environ.get("P2_SLOW") is None, reason="P2_SLOW=1")


def _load(key):
    for f in sorted(RES.glob("kat_float_*.json")):
        d = json.loads(f.read_text())
        if key in d:
            return d[key]
    pytest.skip(f"{key} not stored")


@pytest.mark.parametrize("key", ["KAT-1", "KAT-2", "KAT-3", "KAT-4", "KAT-5", "KAT-7"])
def test_stored_outcome_passes(key):
    """Every KAT passed its pre-registered criterion (p2/kats.py docstring) in the S1 run."""
    d = _load(key)
    assert d["pass"] is True


def test_kat3_printed_digits():
    d = _load("KAT-3")
    assert d["round4"] == [0.992, 0.9967] and d["at_q1"]["PS_dominant"]


def test_kat7_agreement_digits():
    d = _load("KAT-7")
    assert all(float(r["rel_diff"]) <= 1e-10 for r in d["rows"]) and len(d["rows"]) == 8


def test_kat1_live_l1_ps():
    g = Geom.frozen("PA")
    v, err = qnm(g, 1, mp.mpc("3.2426164126", "-0.7958326323"), dps=40, K=6)
    assert abs(mp.re(v) - mp.mpf("3.2426164126")) <= mp.mpf("5e-11")
    assert abs(mp.im(v) + mp.mpf("0.7958326323")) <= mp.mpf("5e-11")
    assert err < mp.mpf("1e-20")


@SLOW
def test_kat1_rerun():
    assert kats.kat1()["pass"]
