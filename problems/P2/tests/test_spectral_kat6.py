"""KAT-6 (validated winding + Krawczyk) and KAT-7 (spectral method) checks.

Fast by default (~5 s): stored KAT-6 certificate fields, spectral double-precision spectrum.
P2_SLOW=1: re-run both KAT-6 boxes (~15 s) and the high-precision spectral refinement (~10 s).
"""
import json
import os
from pathlib import Path

import mpmath as mp
import numpy as np
import pytest

from p2 import winding
from p2.qnm_float import Geom
from p2.spectral import Collocation

SLOW = pytest.mark.skipif(os.environ.get("P2_SLOW") is None, reason="P2_SLOW=1")
CERT = Path(winding.CERT_DIR) / "kat6_winding_krawczyk.json"


def test_kat6_certificate_fields():
    c = json.loads(CERT.read_text())
    assert c["pass"]
    for case in c["cases"]:
        assert case["winding"]["certified"] and case["winding"]["winding"] == 1
        assert case["krawczyk"]["certified_unique_zero"] and case["krawczyk"]["float_value_in_Z"]
    assert c["cases"][0]["printed_within_5e-11"]


@SLOW
def test_kat6_rerun():
    r = winding.run_kat6(write=False)          # replay must not overwrite the committed certificate
    assert r["pass"]


def test_spectral_double_reproduces_PB():
    mp.mp.dps = 40
    g = Geom.frozen("PB")
    ev0 = Collocation(g, 0, 40, prec=128).eig_double()
    ev1 = Collocation(g, 1, 40, prec=128).eig_double()
    for ev, z in ((ev0, -0.909381j), (ev0, -1.855200j), (ev1, -1.911925j), (ev1, 7.214864 - 1.753306j)):
        assert np.min(np.abs(ev - z)) < 2e-6
    assert np.min(np.abs(ev0)) < 1e-8          # the l = 0 zero mode


@SLOW
def test_spectral_refined_l1_ps():
    mp.mp.dps = 40
    g = Geom.frozen("PB")
    c = Collocation(g, 1, 60, prec=200)
    _, wb, _ = c.refine(7.214864 - 1.753306j)
    assert abs(float(wb.real.mid()) - 7.214863916658823) < 1e-13
