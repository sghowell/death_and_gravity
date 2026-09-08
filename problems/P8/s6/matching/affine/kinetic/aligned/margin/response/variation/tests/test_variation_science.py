"""Actual canonical variation and nontrivial initial-state controls."""
import sympy as sp
from p8_vector_variation import modes, preparation, readouts


def assert_zero(values):
    for value in values.values():
        assert all(entry == 0 for entry in value) if isinstance(value, sp.MatrixBase) else value == 0


def test_complete_homogeneous_vector_reduction_and_mass_variations():
    assert_zero(modes.checks())


def test_physical_Cauchy_impedance_and_initial_mixing():
    assert_zero(preparation.checks())


def test_actual_source_has_nonzero_initial_lapse_and_wrong_frozen_covariance():
    assert all(value is True for value in preparation.proof_checks().values())


def test_derivative_free_phase_form_retains_time_dependent_canonical_map():
    data = modes.symplectic_variation()
    assert data["original_variation_has_no_metric_derivatives"] is True
    assert data["two_causal_variation_forms_agree"] == sp.zeros(2)


def test_transverse_leading_physical_data_are_a_conformal_control():
    value = preparation.impedances()["ordinary_transverse_physical_impedance"]
    assert sp.factor(value-modes.scale*sp.sqrt(modes.q)) == 0


def test_initial_impedance_changes_for_a_true_nearby_negative_lapse():
    p = preparation.kernel.geometry.P
    ratio2 = preparation.impedances()["longitudinal_impedance_ratio_squared"]
    assert ratio2.subs(p, sp.Rational(1, 2)) == 1
    assert 0 < ratio2.subs(p, sp.Rational(501, 1000)) < 1
    assert preparation.cauchy_mismatch(sp.Rational(1, 2))["leading_mismatch_vanishes"] is True
    assert preparation.cauchy_mismatch(sp.Rational(501, 1000))["leading_mismatch_vanishes"] is False


def test_initial_force_lower_keeps_both_state_and_local_matching():
    data = preparation.initial_force()
    assert data["state_error_upper"] > 0
    assert data["local_force_coefficients"][2] < 0
    assert data["initial_lapse_force_lower"] > sp.Rational(1, 2*10**14)


def test_actual_energy_pressure_and_mass_second_jet_contact_terms():
    assert_zero(readouts.checks())


def test_contact_matrix_is_not_dropped_from_covariance_response():
    energy = readouts.matrices()["transverse"]["energy_contact_matrix"]
    assert any(value.has(readouts.beta2) for value in energy)
    assert energy != sp.zeros(2)
