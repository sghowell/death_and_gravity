"""Exact finite triangle coefficient and strict comparison-scope tests."""

import pytest
import sympy as s
from p8_vacuum_affine_triangle_curvature_coefficient import (
    audit,
    calibration,
    dressing,
    jets,
    moment,
    source,
    triangle,
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


def test_equal_mass_coefficient_is_exact_not_zero_default():
    assert moment.coefficient(1) == -s.Rational(1, 630)


@pytest.mark.parametrize("mass", (s.Rational(1, 2), 2, 3, 10))
def test_exact_positive_mass_coefficient(mass):
    actual = moment.coefficient(mass)
    assert not actual.has(s.Float)
    assert actual == moment.CLOSED.subs(moment.N, mass)


def test_generic_entire_simplex_polynomial():
    data = triangle.data()["whole_entire_difference_by_derivative_order"]
    assert data[2] == data[4] == 0
    assert (
        s.expand(
            data[6]
            + 8
            * triangle.X**2
            * triangle.Y**2
            * (1 - triangle.X - triangle.Y) ** 2
            * triangle.Q
        )
        == 0
    )


@pytest.mark.parametrize("split", (0, 1, 2))
def test_each_literal_metric_position_nonzero(split):
    assert triangle.literal_line(split, 3) != 0


def test_dressing_uses_shifted_outer_triangle_and_all_labels():
    packet = dressing.data()
    assert "C_triangle(u;w,q)" in packet["whole_complete_dressed_identity"]
    assert "24" in packet["whole_first_order_normalization"]
    assert "double count" in packet["whole_matching_boundary"]


def test_original_source_and_frontiers_unchanged():
    n, g = source.HEAVY_MASS2, source.CUBIC
    assert source.CONTACT == -g * g * (3 / (n - 2) - 2 / (n - 2) ** 2)
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 200
    assert audit.matching()[:-1] == audit.previous.matching()
    assert len(audit.qualifications()) == 6


def test_actual_module_counts_and_original_kernel_calibrations():
    assert [
        len(module.data()["checks"])
        for module in (source, jets, triangle, moment, dressing, calibration)
    ] == [279, 9, 18, 13, 9, 186]
    counts = calibration.data()["whole_independent_counts"]
    assert counts["vector_order_checks"] == counts["frozen_kernel_checks"] == 54
    assert counts["nonzero_curvature_words"] > 0


def test_known_coefficient_not_full_parent_matching():
    text = audit.observable()["not_established"]
    assert "independent extra parent chi" in text
    assert "box coefficient" in text and "P8 remain open" in text


def test_local_jet_does_not_replace_physical_above_threshold_amplitude():
    text = moment.data()["whole_original_known_coefficient_sign_and_bound"]
    assert "NOT a truncation-error bound" in text
    assert "uniformly" in jets.data()["whole_analytic_jet_proof"]
