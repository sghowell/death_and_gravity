"""Pole and finite-local normalization controls before promotion."""
import sympy as sp
from p8_vector_dimensional import adm, heat, local, mass, metric, replay


def test_dimensionally_continued_local_coefficients_match_known_poles():
    assert all(value == 0 for value in local.checks().values())


def test_dimension_dependent_action_and_canonical_normalization():
    assert all(value == 0 for value in adm.variations().values())


def test_evanescent_polarization_and_normalization_cannot_be_dropped():
    assert local.integrated()[0]["evanescent_finite_contribution"] != 0


def test_covariant_counterterm_variation_precedes_four_dimensional_limit():
    assert all(value == 0 for value in metric.ordinary_matching()["identities"].values())


def test_componentwise_finite_part_alone_is_not_a_covariant_prescription():
    d = metric.ordinary_matching()["ordinary_covariantly_subtracted_local_terms"]
    assert d[1]["radial_Ward_defect"] == 8*local.H*local.Hd
    assert sp.factor(d[1]["energy"]-local.H**2*(6*local.ell-10)) == 0


def test_actual_clock_mass_poles_have_an_independent_covariant_reconstruction():
    assert all(value == 0 for value in mass.checks().values())


def test_frozen_finite_flat_potential_jets_are_not_changed():
    assert all(value == 0 for value in mass.frozen_flat_finite_match().values())


def test_general_dimension_reduces_to_actual_frozen_WKB_and_observable():
    assert all(value == 0 for value in replay.physical_checks().values())


def test_radial_Gamma_residue_finite_term_and_angular_continuation():
    assert all(value == 0 for value in replay.gamma_checks().values())


def test_independent_dimensional_heat_action_reduces_to_frozen_curvature_input():
    assert all(value == 0 for value in heat.checks().values())
