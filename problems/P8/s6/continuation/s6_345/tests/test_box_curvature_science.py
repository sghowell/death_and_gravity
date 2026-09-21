"""Whole generic ordered-box matching and strict comparison-scope tests."""

import pytest
import sympy as s
from p8_vacuum_affine_box_curvature_coefficient import (
    audit,
    basis,
    calibration,
    jets,
    moment,
    radiation,
    source,
)

ROWS, GATES = audit.residuals(), audit.gates()


@pytest.mark.parametrize("name", tuple(ROWS))
def test_every_exact_residual(name):
    value = ROWS[name]
    assert all(
        v == 0 for v in (list(value) if isinstance(value, s.MatrixBase) else [value])
    )
    assert not s.sympify(value).has(s.Float)


@pytest.mark.parametrize("name", tuple(GATES))
def test_every_proof_gate(name):
    assert GATES[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[row[0] for row in audit.bad_cases()]
)
def test_all_rejected_inputs_and_scopes(name, call, args):
    assert name
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize("index", range(6))
def test_fixed_complete_word_basis(index):
    assert audit.require_word_index(index) == basis.WORDS[index]


@pytest.mark.parametrize("mass", (s.Rational(1, 2), 2, 3, 10))
def test_exact_mass_coefficient(mass):
    actual = moment.coefficient(mass)
    assert not actual.has(s.Float)
    assert actual == moment.CLOSED.subs(moment.N, mass)


def test_equal_mass_value_not_a_zero_default():
    assert moment.coefficient(1) == -s.Rational(58, 945)


def test_full_flat_basis_and_all_coefficient_projection():
    packet = basis.data()
    assert packet["whole_flat_full_monomial_count"] == 220
    assert packet["whole_flat_Bose_orbit_count"] == 20
    assert packet["whole_flat_basis_rank"] == 6
    assert len(packet["whole_box_flat_lift_coefficients"]) == 6


def test_entire_generic_radiative_identity_not_fixtures():
    packet = radiation.data()
    assert packet["whole_generic_coefficient_count"] == 2430
    assert packet["whole_radial_coefficient_count"] == 1866
    assert packet["checks"]["whole_generic_every_radiative_coefficient"] == s.zeros(
        2430, 1
    )


def test_literal_connections_and_original_convention():
    packet = jets.data()
    assert len(packet["checks"]) == 67
    assert packet["checks"]["whole_original_Galileon_metric_and_connections"] == 0
    assert "connection" in packet["whole_third_jet_identity"]


def test_original_source_and_frontiers_unchanged():
    n, g = source.HEAVY_MASS2, source.CUBIC
    assert source.CONTACT == -g * g * (3 / (n - 2) - 2 / (n - 2) ** 2)
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 201
    assert audit.matching()[:-1] == audit.previous.matching()
    assert len(audit.qualifications()) == 6


def test_actual_module_and_calibration_counts():
    assert [
        len(module.data()["checks"])
        for module in (source, basis, jets, radiation, moment, calibration)
    ] == [283, 228, 67, 4, 15, 156]
    counts = calibration.data()["whole_original_counts"]
    assert (
        counts["entire_vector_polynomials"]
        == counts["nonzero_curvature_polynomials"]
        == 6
    )
    assert counts["frozen_kernel_checks"] == 72


def test_triangle_conversion_and_extra_parent_matching_not_claimed():
    text = audit.observable()["not_established"]
    assert "triangle convention must first be converted" in text
    assert "Independent extra parent chi" in text
    assert "P8 remain open" in text


def test_large_coefficient_bound_not_physical_perturbative_smallness():
    text = moment.data()["whole_original_coefficient_bound"]
    assert "LARGE" in text and "not perturbative smallness" in text
    assert "different" in moment.data()["whole_scheme_boundary"]
