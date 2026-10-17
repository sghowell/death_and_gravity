"""S1 item 5: P2-2 thresholds W_l, l <= 50, at P_B and P_A (certificates/lg/*.json).

Fast (~10 s): every certificate is replayed (criterion recomputed in Arb at the stored W_l and
matching point).  The float cross-check compares the rigorous LG error bound with the actual
Jost Wronskian of the independent float solver.
"""
import json
from pathlib import Path

import mpmath as mp
import pytest
from flint import arb, fmpq

from p2 import lg_bounds
from p2.params import Point
from p2.qnm_float import Geom, Solver

CERT = Path(lg_bounds.CERT_DIR)


@pytest.mark.parametrize("name", ["PB", "PA"])
def test_replay_all_l(name):
    assert lg_bounds.verify(name, 50) == []


@pytest.mark.parametrize("name", ["PB", "PA"])
def test_threshold_shape(name):
    """W_l/((l+1/2) Omega_P) > 1 (above the barrier), decreasing towards 1, and below the
    planning value 1.25 (l+1/2) Omega_P + 2 kappa_- for l >= 3."""
    ratios = []
    for l in range(51):
        st = json.loads((CERT / f"lg_{name}_l{l}.json").read_text())
        assert st["criterion"]["certified"]
        r = float(st["W_l_over_l_half_OmegaP"].split()[0].strip("["))
        ratios.append(r)
        assert r > 1
        wl = float(st["W_l_over_kappa_minus"].split()[0].strip("["))
        assert wl > st["sqrt_Vmax_plus_beta2_over_kappa_minus_float"]
        if l >= 3:
            assert wl < st["planning_value_over_kappa_minus"]
    assert ratios[50] < 1.04 and ratios[10] < 1.13


def test_criterion_monotone_in_threshold():
    pt = Point("PB")
    lg = lg_bounds.LG(pt, 4)
    st = json.loads((CERT / "lg_PB_l4.json").read_text())
    a0 = fmpq(*[int(x) for x in st["W_l_exact"].split("/")])
    rX = fmpq(*[int(x) for x in st["matching_point_exact"].split("/")])
    c1 = lg.criterion(a0, rX)["value"]
    c2 = lg.criterion(a0 * fmpq(6, 5), rX)["value"]
    assert c2 < c1 < arb(2)


@pytest.mark.parametrize("name,l", [("PB", 2), ("PA", 10)])
def test_bound_dominates_actual_wronskian(name, l):
    """|W_l e^{i Theta}/(2 i w) - 1| <= criterion - 1 at Re w = W_l, Im w = -beta (float check of
    the whole LG chain: normalisation, phase Theta = int (p - w) dr*, and the bound)."""
    mp.mp.dps = 30
    g = Geom.frozen(name)
    st = json.loads((CERT / f"lg_{name}_l{l}.json").read_text())
    crit = mp.mpf(st["criterion"]["value_upper"].split()[0].strip("["))
    beta = (mp.mpf(1) / 2 + mp.mpf(st["delta"].split("/")[0]) / mp.mpf(st["delta"].split("/")[1])) * g.km
    w = mp.mpc(st["W_l_float"], -beta)
    p = lambda r: mp.sqrt(w * w - g.V(r, l))
    theta = mp.quad(lambda r: (p(r) - w) / g.F(r), [g.rp, g.rP, g.rc])
    W = Solver(g).wronskian(w, l)
    dev = abs(W * mp.exp(1j * theta) / (2j * w) - 1)
    assert dev <= crit - 1
    assert dev < mp.mpf("0.05")          # the bound is conservative (documented in s1-report)
