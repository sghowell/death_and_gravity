"""Dimensional geometry, exact contractions and finite-counterterm controls."""
import pytest
import sympy as sp
from p8_vector_quadratic_matching import (
    bimetric,
    contractions,
    evanescence,
    geometry,
    kernel,
    response,
    verify,
)


def zeros(values):
    assert all(value == 0 for value in values.values())


def test_warped_product_geometry_returns_four_dimensions():
    zeros(geometry.checks())


@pytest.mark.parametrize("rank", (1, 2, 3, 4))
def test_exact_dimensional_assignment_counts(rank):
    assert contractions.checks()["exact_index_partition_count_rank_"+str(rank)] == 0


def test_every_continued_second_invariant_matches():
    zeros(contractions.checks())


def test_complete_auxiliary_density_and_bilinear_match():
    zeros(bimetric.checks())


@pytest.mark.parametrize("order", (2, 4))
def test_actual_unreduced_and_compact_dimensional_variations(order):
    zeros(response.checks(order))


@pytest.mark.parametrize("order", (2, 4))
def test_bare_counterterm_finite_sign(order):
    assert response.finite_counterterm_operator(order) == 2*response.operator_jet(order, 1)


def test_dimension_dependent_measure_cannot_be_frozen_early():
    assert response.measure_control(4)["bounce_fixture"] == -sp.Rational(1024, 6561)
    assert response.measure_control(2)["omitted_measure_variation"] != 0


def test_gauss_bonnet_finite_boundary_control():
    value = evanescence.gauss_bonnet_control()
    assert value["exact_curvature_identity"] == value["exact_boundary_identity"] == value["scalar_compact_variation_identity"] == 0
    assert value["nonzero_dimension_derivative_fixture"] == sp.Rational(1, 15)
    assert value["topological_factor_has_simple_zero_at_D_three"] == sp.Rational(1, 90)


def test_native_dimensional_jet_validation_after_warmup():
    kernel.dimension_jet(2, 0)
    with pytest.raises(ValueError):
        kernel.dimension_jet(2, False)
    with pytest.raises(ValueError):
        kernel.dimension_jet(sp.Integer(2), 0)


def test_unsupported_inputs_are_rejected():
    assert verify.controls()["rejected_inputs"] == 189


def test_exact_polynomial_and_unchanged_clock_scope():
    assert all(value is True for value in verify.proof_checks().values())
