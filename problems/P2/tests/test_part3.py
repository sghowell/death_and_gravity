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
from flint import acb, arb, fmpq

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
    # (7, 40.0, -0.55) is the point registered at 89a2336b; (7, 30.0, -0.55) was added at 3e7fca05
    out = part3.zsol_vs_arbw("PB", [(0, 2.0, -0.55), (3, 12.0, 0.1), (7, 40.0, -0.55), (7, 30.0, -0.55)])
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


def _box_check(l, rho_kappa, corners):
    """zsol's box enclosure (the source of every sup bound M) over [c +- rho]^2 at the corner
    (a_l, -b) of R'_l must contain S1's arbw enclosure at the given box corners (arbw at the lowest
    precision with relative radius < 1e-10; an independent code)."""
    pt = Point("PB", 256)
    st = json.loads((CERT / f"winding_PB_l{l}.json").read_text())
    a, b = part3._fq(st["region"]["a"]), part3._fq(st["region"]["b"])
    rho = fmpq(round(float(pt.km.mid()) * rho_kappa * 2 ** 70), 2 ** 70)
    with precision(part3.PREC):
        c = acb(arb(a), arb(-b))
        box = acb(c.real + arb(0, 1) * arb(rho), c.imag + arb(0, 1) * arb(rho))
        _, Wb = zsol.wronskian_z(pt, l, c, box, theta=part3.THETA_BOX)
        Wc = zsol.wronskian_z(pt, l, c, theta=part3.THETA_POINT)[0]
        relw = float(max(Wb.real.rad(), Wb.imag.rad()) / abs(Wc).mid())
    for sx, sy in corners:
        Wa, prec, rr = part3._arbw_reference("PB", l, (a + sx * rho, -b + sy * rho))
        assert rr < 1e-10
        with precision(256):
            assert Wb.contains(Wa), (l, rho_kappa, sx, sy)
    return relw


@pytest.mark.skipif(19 not in LS, reason="no l = 19 certificate")
def test_box_enclosure_contains_arbw_l19_corner():
    """Review S4: a tight box (half-width 1e-4 kappa_-) at the highest-frequency corner.  The corner
    values lie 1.04e-4 (relative) from the centre value and the box ball's half-width is 5.8e-3, so
    a deviation bound too small by a factor of about 60 or more fails this test."""
    relw = _box_check(19, 1e-4, [(-1, -1), (1, 1)])
    assert relw < 0.1


@SLOW
@pytest.mark.skipif(19 not in LS, reason="no l = 19 certificate")
def test_box_enclosure_l19_corner_wide():
    """The same at half-width kappa_-/16, all four box corners."""
    _box_check(19, 1 / 16, [(-1, -1), (1, -1), (1, 1), (-1, 1)])


@SLOW
@pytest.mark.skipif(19 not in LS, reason="no l = 19 certificate")
def test_segment_replay_l19():
    """Review S4: two stored l = 19 segments recomputed with their stored K (the bottom-edge segment
    at the corner, the most expensive one, and the middle of the top edge)."""
    st = json.loads((CERT / "winding_PB_l19.json").read_text())
    pcs = st["pieces"]
    corner = max((r for r in pcs if r["edge"] == 0), key=lambda r: r["idx"])
    top = [r for r in pcs if r["edge"] == 2]
    for r in (corner, top[len(top) // 2]):
        out = part3._replay_task(("PB", 19, r))
        assert out["overlap"], out


@SLOW
def test_control_rectangles():
    for name, l, box, exp, _ in [part3.C2_RECTS[2]] + part3.C3_RECTS:
        r = part3.rectangle_winding(name, l, box, workers=1)
        assert r["certified"] and r["winding"] == exp, r


@SLOW
@pytest.mark.parametrize("l", [l for l in LS if l <= 5])
def test_replay_sample(l):
    assert part3.verify("PB", [l]) == []
