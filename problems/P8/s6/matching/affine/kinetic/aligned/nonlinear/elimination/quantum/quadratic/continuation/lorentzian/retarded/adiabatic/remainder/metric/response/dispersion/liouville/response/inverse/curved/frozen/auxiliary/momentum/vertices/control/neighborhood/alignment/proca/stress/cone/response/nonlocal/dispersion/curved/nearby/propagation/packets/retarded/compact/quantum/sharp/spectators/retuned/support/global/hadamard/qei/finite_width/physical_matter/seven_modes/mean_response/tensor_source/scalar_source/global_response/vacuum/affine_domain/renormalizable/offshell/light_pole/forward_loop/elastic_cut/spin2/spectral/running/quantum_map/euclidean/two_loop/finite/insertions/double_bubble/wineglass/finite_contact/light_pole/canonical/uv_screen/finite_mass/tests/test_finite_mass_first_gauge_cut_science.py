"""Independent tensor, composition, rational-bound and negative-control tests."""

from fractions import Fraction
from math import comb

import pytest
import sympy as sp
from p8_vacuum_finite_mass_gauge_cut import (
    audit,
    calibration,
    cuts,
    dirac,
    kinematics,
    series,
)


@pytest.mark.parametrize("name", list(audit.residuals()))
def test_every_named_exact_identity(name):
    value = audit.residuals()[name]
    assert (
        all(v == 0 for v in value) if isinstance(value, sp.MatrixBase) else value == 0
    )


@pytest.mark.parametrize("name", list(audit.gates()))
def test_every_explicit_gate(name):
    assert audit.gates()[name] is True


@pytest.mark.parametrize(
    "name,call,args", cuts.bad_cases(), ids=[x[0] for x in cuts.bad_cases()]
)
def test_every_invalid_exact_input(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize("n", range(13))
def test_composition_counts_with_independent_nested_integer_loops(n):
    actual = series.compositions(n, 4)
    independent = {
        (i, j, k, n - i - j - k)
        for i in range(n + 1)
        for j in range(n + 1 - i)
        for k in range(n + 1 - i - j)
    }
    assert set(actual) == independent
    assert len(actual) == len(independent) == comb(n + 3, 3)


@pytest.mark.parametrize("args", [(True, 4), (1.0, 4), (-1, 4), (1, 0)])
def test_invalid_composition_domain(args):
    with pytest.raises((TypeError, ValueError)):
        series.compositions(*args)


def test_explicit_Dirac_spectrum_and_transverse_norm():
    gamma = dirac.data()["gamma_matrices"]
    v = [sp.Rational(3, 5), sp.Rational(4, 5), 0, 0]
    slash = sum((v[i] * gamma[i] for i in range(4)), sp.zeros(4))
    assert slash.eigenvals() == {-1: 2, 1: 2}
    assert slash.conjugate().T * slash == sp.eye(4)


def test_complete_degree_two_tensor_system_has_only_expected_kernel():
    d = kinematics.data()["degree_two_complete_parity_even_Ward_system"]
    A = d["linear_constraint_matrix"]
    null = A.nullspace()
    assert A.shape == (36, 15)
    assert len(null) == 1
    assert sp.Matrix.hstack(null[0], d["one_dimensional_solution"]).rank() == 1
    # eta*k^2 alone is not transverse to both independent gauge momenta.
    wrong = sp.Matrix([1] + [0] * 14)
    assert any(v != 0 for v in A * wrong)


def test_all_routes_close_and_have_no_more_than_three_external_shifts():
    d = kinematics.data()
    labels = {tuple(r["vertices"]) for r in d["six_cyclic_box_routes"]}
    assert len(labels) == 6
    for route in d["six_cyclic_box_routes"]:
        assert len(route["partial_momenta"]) == 4
        assert route["partial_momenta"][0] == sp.zeros(4, 1)
        end = (
            route["partial_momenta"][-1]
            + d["Euclidean_all_outgoing_external_momenta"][route["vertices"][-1]]
        )
        assert end.applyfunc(sp.factor) == sp.zeros(4, 1)


@pytest.mark.parametrize("x", (Fraction(0), Fraction(1, 4), Fraction(1, 2)))
def test_tail_majorant_fraction_controls(x):
    n3_tail = x**3 * (8 - 7 * x) / (6 * (1 - x) ** 2)
    assert n3_tail <= Fraction(16, 3) * x**3 <= 6 * x**3


@pytest.mark.parametrize("s", (1, Fraction(3, 2), 2, Fraction(5, 2), 3))
def test_cut_bounds_using_independent_Fraction_arithmetic(s):
    mass = 10**200
    eps = Fraction(10**6, mass)
    expected_error = 18 * (eps + eps**2 / 8)
    D = Fraction(s) ** 2 - (4 - Fraction(s)) ** 2
    result = cuts.point(s, mass)
    assert result["leading_boundary_difference_over_dA_z2_div_9pi"] == sp.Rational(D)
    assert result[
        "absolute_boundary_difference_error_over_dA_z2_div_9pi_upper"
    ] == sp.Rational(expected_error)
    assert result["strict_nonzero_boundary_difference_at_this_order"] == (
        abs(D) > expected_error
    )


def test_crossing_sign_not_added_in_same_direction():
    a = cuts.point(Fraction(3, 2))
    b = cuts.point(Fraction(5, 2))
    key = "leading_boundary_difference_over_dA_z2_div_9pi"
    assert a[key] == -4
    assert b[key] == 4
    assert (
        a["relative_boundary_difference_error_upper"]
        == b["relative_boundary_difference_error_upper"]
    )


def test_center_has_no_relative_division_by_zero():
    d = cuts.point(2)
    assert d["relative_boundary_difference_error_upper"] is None
    assert d["strict_nonzero_boundary_difference_at_this_order"] is False


def test_valid_inconclusive_points_are_not_rejected():
    assert (
        cuts.point(mass=24)["strict_nonzero_boundary_difference_at_this_order"] is False
    )
    assert (
        cuts.point(2 + Fraction(1, 10**200))[
            "strict_nonzero_boundary_difference_at_this_order"
        ]
        is False
    )


def test_actual_strict_interior_bound_is_finite_mass_not_field_series():
    d = calibration.data()
    assert (
        0
        < d["first_two_gauge_cut_relative_error_upper_at_interior"]
        < sp.Rational(5, 10**194)
    )
    assert d["fermion_mass"] == 10**200
    assert d["interior"]["strict_nonzero_boundary_difference_at_this_order"] is True


def test_counts_and_nonclosure():
    assert len(audit.residuals()) == 190
    assert audit.scalar_entries() == 606
    assert len(audit.gates()) == 34
    assert len(audit.controls()) == 9
    assert audit.rejected_inputs() == 31
    assert audit.controls()["original_V_G_B_and_P8_not_closed"] is True
