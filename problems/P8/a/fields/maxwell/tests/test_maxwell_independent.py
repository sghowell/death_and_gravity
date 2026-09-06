"""Independent rational contractions and controls against forbidden imports."""

from fractions import Fraction
from pathlib import Path

import pytest
import sympy as sp
from p8a_maxwell import bounds, functional, independent, stress


def test_fraction_covariant_reconstruction_matches_every_primary_component():
    assert independent.replay()["stress_jet_coefficients"] == stress.jet_coefficients()


def test_fraction_radiation_control_keeps_anomaly_and_removes_only_R_squared():
    data = independent.replay()["radiation_control"]
    assert data == independent.serialize(functional.radiation_control())
    assert data["reference_EED_pi_squared_t_fourth_over_hbar"] == "31/2560"
    assert data["absolute_zero_times_t_fourth"] == "4601/1280"
    assert data["absolute_zero_times_t_fourth"] != data["squared_zero_times_t_fourth"]


def test_fraction_envelope_arithmetic_requires_scheme_and_does_not_infer_cosmology():
    assert independent.replay()["envelope_controls"] == independent.serialize(bounds.calibration())
    data = independent.replay()["focusing_incompatibility"]
    assert Fraction(data["K_trace_times_tau_cap"]) < Fraction(data["A1_zeta_zero_gradient_lower"])
    assert data["small_Q2_implies_focusing"] is False


def test_independent_generic_IBP_polynomials_match_symbolic_proper_derivation():
    t = sp.Symbol("t", real=True)
    h = sp.Function("H")(t)
    jets = (h, sp.diff(h, t), sp.diff(h, t, 2), sp.diff(h, t, 3))
    for name, expression in functional.ibp_coefficients(h, t).items():
        poly = sp.Poly(expression, *jets)
        expected = {",".join(map(str, powers)): str(value) for powers, value in poly.terms() if value}
        assert independent.replay()["proper_IBP_coefficients"][name] == expected


def test_independent_engine_is_fraction_only_and_does_not_hide_a_missing_jet():
    source = Path(independent.__file__).read_text()
    assert "import sympy" not in source
    assert "from ." not in source
    with pytest.raises(ValueError, match="fourth proper Hubble jet"):
        independent.derivative(independent.jet(3))


def test_all_source_sign_coarsenings_are_upward_and_negative_controls_nonzero():
    cases = independent.source_controls()
    assert all(case["upward_coarsening_gap"] >= 0 for case in cases.values())
    assert cases["Lambda=1,other_lower=-2"]["raw_source"] == 5
    assert cases["Lambda=-1,other_lower=2"]["raw_source"] == -5
    assert cases["Lambda=-1,other_lower=2"]["coarsened_source"] == 0
    assert cases["Lambda=1,other_lower=-2"]["coarsened_source"] != 1


def test_energy_is_not_EED_even_for_zero_type_D_de_sitter():
    polynomials = independent.covariant_polynomials()
    jets = (1, 0, 0, 0)
    assert independent.evaluate(polynomials["rho"][0], jets) == 186
    assert independent.evaluate(polynomials["EED"][0], jets) == -186
    assert independent.evaluate(polynomials["I_EED"][0], jets) == 0
