"""S3 part (iii): the phase-factored W (zsol), the Cauchy-disc winding pipeline (part3) and its
certificates (certificates/part3/; criteria in notes/s3-preregistration.md).

Fast (~20 s): the sympy identities; zsol against the S1 arbw at points; D4; the l = 0 Krawczyk
enclosure against KAT-6; every stored winding certificate re-summed with its coverage checked and
the expected count; two stored segments of l = 0 recomputed.
P2_SLOW=1: the control rectangles C2(c), C3, and a 5% replay sample for l <= 5.
"""
import json
import os
from pathlib import Path

import pytest
from flint import acb, arb

from p2 import part3, zeq_identity, zsol
from p2.arbutil import precision
from p2.params import Point

CERT = part3.CERT_DIR
SLOW = pytest.mark.skipif(os.environ.get("P2_SLOW") is None, reason="P2_SLOW=1")
LS = sorted(int(p.stem.split("_l")[1]) for p in CERT.glob("winding_PB_l*.json")) if CERT.exists() else []


def test_zeq_identity():
    r = zeq_identity.checks()
    assert r["pass"], r


def test_zsol_matches_arbw():
    out = part3.zsol_vs_arbw("PB", [(0, 2.0, -0.55), (3, 12.0, 0.1), (7, 30.0, -0.55)])
    assert all(x["pass"] for x in out), out


def test_d4():
    assert part3.d4_check("PB")["pass"]


def test_krawczyk_l0_matches_kat6():
    kr = part3.krawczyk_z("PB", 0, ("0", "-0.9093809789948054555249"))
    k6 = json.loads((part3.ROOT / "certificates" / "kat6_winding_krawczyk.json").read_text())["cases"][1]["krawczyk"]
    assert kr["certified_unique_zero"]
    with precision(256):
        assert arb(kr["enclosure_over_kappa_minus"]["im"]).overlaps(arb(k6["enclosure_over_kappa_minus"]["im"]))


@pytest.mark.parametrize("l", LS)
def test_certificate_resum_coverage_expected(l):
    assert part3.verify("PB", [l], sample_frac=0.0, min_sample=0) == []
    st = json.loads((CERT / f"winding_PB_l{l}.json").read_text())
    assert st["certified"] and st["winding"] == (1 if l == 0 else 0)


@pytest.mark.skipif(0 not in LS, reason="no l = 0 certificate")
def test_segment_replay_l0():
    st = json.loads((CERT / "winding_PB_l0.json").read_text())
    pt = Point("PB", 256)
    with precision(part3.PREC):
        for r in (st["pieces"][0], st["pieces"][len(st["pieces"]) // 2]):
            seg = {"edge": r["edge"], "idx": r["idx"], "p0": [part3._fq(x) for x in r["p0"]],
                   "p1": [part3._fq(x) for x in r["p1"]], "c": [part3._fq(x) for x in r["c"]],
                   "rho": part3._fq(r["rho"])}
            rr = part3.cauchy_segment(pt, 0, seg, K=r["K"])
            assert arb(rr["increment"]).overlaps(arb(r["increment"]))


@SLOW
def test_control_rectangles():
    for name, l, box, exp, _ in [part3.C2_RECTS[2]] + part3.C3_RECTS:
        r = part3.rectangle_winding(name, l, box, workers=1)
        assert r["certified"] and r["winding"] == exp, r


@SLOW
@pytest.mark.parametrize("l", [l for l in LS if l <= 5])
def test_replay_sample(l):
    assert part3.verify("PB", [l]) == []
