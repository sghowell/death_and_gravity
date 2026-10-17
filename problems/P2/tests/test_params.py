"""S1 item 1: certified parameters (certificates/params.json) -- replay and known answers.

Fast (~1 s).  Run: .venv/bin/python -m pytest problems/P2/tests -q   (from the repo root)
"""
import json
from pathlib import Path

import pytest
from flint import arb

from p2 import params
from p2.arbutil import precision

CERT = Path(params.CERT_DIR) / "params.json"


@pytest.fixture(scope="module")
def stored():
    return json.loads(CERT.read_text())


@pytest.fixture(scope="module")
def points():
    return {nm: params.Point(nm) for nm in ("PB", "PA")}


def test_replay_identical(stored):
    """Recomputing every ball reproduces the stored certificate exactly."""
    assert params.verify(CERT) == []


@pytest.mark.parametrize("name", ["PB", "PA"])
def test_all_required_inequalities_certified(stored, name):
    ineq = stored["points"][name]["inequalities"]
    required = [k for k in ineq if k != "CGNS18_rho"]
    assert all(ineq[k]["certified"] for k in required), [k for k in required if not ineq[k]["certified"]]
    assert stored["points"][name]["defining_polynomials"]["all_roots_located"]


def test_cgns_context_split(stored):
    """CGNS18's rho < 9/7 holds at P_B and is certified FALSE at P_A (context row of section 3)."""
    assert stored["points"]["PB"]["inequalities"]["CGNS18_rho"]["certified"]
    assert stored["points"]["PA"]["inequalities"]["CGNS18_rho"]["certified_false"]


# S0 float table (FORMULATION section 1.3), printed digits
S0_TABLE = {
    "PB": {"Q": "1.00984316771", "r_n": "-7.96033597", "r_minus": "0.99736659", "r_plus": "1.09434835",
           "r_c": "5.86862103", "kappa_minus": "0.04254199", "kappa_plus": "0.03500736", "kappa_c": "0.09338246",
           "r_P": "1.95873037", "Omega_P_over_kappa_minus": "4.91883525", "gamma0": "0.07160718",
           "gamma0_over_kappa_minus": "1.683212"},
    "PA": {"Q": "1.00681060565", "r_n": "-7.96006333", "r_minus": "0.94927598", "r_plus": "1.14311753",
           "r_c": "5.86766982", "kappa_minus": "0.09426070", "kappa_plus": "0.06379974", "kappa_c": "0.09332621",
           "r_P": "1.97187372", "Omega_P_over_kappa_minus": "2.20956624", "gamma0": "0.07204371",
           "gamma0_over_kappa_minus": "0.764303"},
}


@pytest.mark.parametrize("name", ["PB", "PA"])
def test_s0_float_values_round_correctly(points, name):
    pt = points[name]
    with precision(256):
        qs = pt.quantities()
        for key, printed in S0_TABLE[name].items():
            nd = len(printed.split(".")[1])
            half_ulp = arb(10) ** (-nd) / 2
            assert abs(qs[key] - arb(printed)) < half_ulp, (key, qs[key].str(15), printed)


def test_table_I_eikonal_known_answer(points):
    """KAT-2 (certified part): gamma0/kappa_- at P_A rounds to 0.7643026660 and lies below the
    CCDHJ18 Table I l = 100 value 0.7643094446 (exact decimal)."""
    pt = points["PA"]
    with precision(256):
        g = pt.gamma0 / pt.km
        assert abs(g - arb("0.7643026660")) < arb("5e-11")
        assert g < arb("0.7643094446")


def test_strip_fractions(points):
    with precision(256):
        assert abs(points["PB"].beta / points["PB"].gamma0 - arb("0.327")) < arb("0.0005")
        assert abs(points["PA"].beta / points["PA"].gamma0 - arb("0.667")) < arb("0.0005")


def test_qmax_reduction_exact():
    """Q_max^2 = r - r^4/50 equals r(3 - r)/2 modulo r^3 - 25 r + 25 (exact)."""
    sp = pytest.importorskip("sympy")
    s = sp.Symbol("s")
    assert sp.rem(sp.expand(s - s**4 / 50 - s * (3 - s) / 2), s**3 - 25 * s + 25, s) == 0
