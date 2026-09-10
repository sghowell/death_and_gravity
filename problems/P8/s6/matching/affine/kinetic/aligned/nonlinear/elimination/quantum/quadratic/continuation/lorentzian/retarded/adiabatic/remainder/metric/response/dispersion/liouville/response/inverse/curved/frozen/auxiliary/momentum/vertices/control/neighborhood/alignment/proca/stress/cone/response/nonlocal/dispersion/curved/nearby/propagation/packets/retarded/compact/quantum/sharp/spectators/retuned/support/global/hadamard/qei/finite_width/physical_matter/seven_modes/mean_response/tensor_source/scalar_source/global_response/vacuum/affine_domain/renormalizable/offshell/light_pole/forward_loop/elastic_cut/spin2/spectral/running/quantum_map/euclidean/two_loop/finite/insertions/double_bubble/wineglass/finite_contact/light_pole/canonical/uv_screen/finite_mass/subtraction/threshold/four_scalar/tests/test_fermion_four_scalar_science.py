"""Independent exact interval, tensor, routing and ownership regressions."""

from fractions import Fraction
from math import comb

import pytest
import sympy as sp
from p8_vacuum_fermion_four_scalar import (
    audit,
    calibration,
    kinematics,
    low_degree,
    reference,
    series,
)
from p8_vacuum_finite_mass_gauge_cut import series as inherited_series


@pytest.mark.parametrize("name", list(audit.residuals()))
def test_every_named_identity(name):
    value = audit.residuals()[name]
    assert (
        all(v == 0 for v in value) if isinstance(value, sp.MatrixBase) else value == 0
    )


@pytest.mark.parametrize("name", list(audit.gates()))
def test_every_explicit_gate(name):
    assert audit.gates()[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[x[0] for x in audit.bad_cases()]
)
def test_invalid_inputs(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize("mass", (36, 10**6, 10**200))
@pytest.mark.parametrize("Y", (Fraction(1, 4), Fraction(1, 10**206)))
def test_tail_upper_using_independent_fraction(mass, Y):
    expected = Fraction(300000000, 144) * Y * Y / mass**4
    assert calibration.tail_bound(mass, Y) == sp.Rational(expected)


@pytest.mark.parametrize("degree", range(9))
def test_composition_count_using_independent_nested_loops(degree):
    points = {
        (i, j, k, degree - i - j - k)
        for i in range(degree + 1)
        for j in range(degree - i + 1)
        for k in range(degree - i - j + 1)
    }
    assert len(points) == comb(degree + 3, 3)
    assert set(inherited_series.compositions(degree, 4)) == points


@pytest.mark.parametrize("x", (Fraction(0), Fraction(1, 4), Fraction(1, 2)))
def test_all_degree_majorant_against_fraction_partial_sums(x):
    majorant = x**4 * (9 - 8 * x) / (6 * (1 - x) ** 2)
    partial = sum(Fraction((n + 3) * (n + 1), 6 * n) * x**n for n in range(4, 101))
    assert partial <= majorant <= Fraction(10, 3) * x**4 <= 4 * x**4


@pytest.mark.parametrize(
    "T,E",
    (
        (1, 0),
        (1, Fraction(1, 4)),
        (1, 1),
        (1, 2),
        (Fraction(1, 10**600), Fraction(1, 10**1204)),
    ),
)
@pytest.mark.parametrize(
    "lo,hi",
    (
        (0, 0),
        (0, Fraction(1, 2)),
        (Fraction(1, 4), Fraction(1, 2)),
        (Fraction(1, 2), Fraction(3, 4)),
    ),
)
def test_canonical_band_all_corners_and_interior_samples(T, E, lo, hi):
    T, E, lo, hi = map(Fraction, (T, E, lo, hi))
    values = [
        (T + error) / (1 - slope) ** 2
        for error in (-E, Fraction(0), E)
        for slope in (lo, (lo + hi) / 2, hi)
    ]
    d = calibration.canonical_band(T, E, lo, hi)
    assert d["lower"] == sp.Rational(min(values))
    assert d["upper"] == sp.Rational(max(values))
    assert d["positive_selected_coefficient"] == (min(values) > 0)
    assert d["strict_selected_increase_over_tree"] == (min(values) > T)


@pytest.mark.parametrize(
    "mass,momentum", ((1, 0), (2, 3), (Fraction(3, 2), Fraction(5, 7)))
)
def test_Dirac_trace_on_axis_from_independent_eigenvalue_pairs(mass, momentum):
    d = low_degree.data()
    m = d["fermion_mass"]
    trace = d["zero_momentum_Dirac_trace"]
    variables = sorted(
        (s for s in trace.free_symbols if str(s).startswith("real_loop_q")), key=str
    )
    sub = {m: sp.Rational(mass), **{s: sp.Integer(0) for s in variables}}
    sub[variables[0]] = sp.Rational(momentum)
    expected = (
        2 * (sp.Rational(mass) - sp.I * sp.Rational(momentum)) ** 4
        + 2 * (sp.Rational(mass) + sp.I * sp.Rational(momentum)) ** 4
    )
    assert sp.simplify(trace.subs(sub) - expected) == 0


def test_complete_quadratic_space_rejects_a_single_external_square():
    d = low_degree.data()
    matrix = d["complete_quadratic_S4_constraint_matrix"]
    assert matrix.shape == (42, 10)
    basis = d["quadratic_S4_invariant_basis"]
    assert sp.Matrix.hstack(*matrix.nullspace(), *basis).rank() == 2
    assert any(v != 0 for v in matrix * sp.Matrix([1] + [0] * 9))


def test_all_routes_and_independent_sign_permutations():
    d = kinematics.data()
    assert len(d["six_complete_forward_routes"]) == 6
    assert len(d["all_twenty_four_external_permutation_checks"]) == 24
    for route in d["six_complete_forward_routes"]:
        assert len(route["partial_momenta"]) == 4
        assert route["partial_momenta"][0] == sp.zeros(4, 1)
    assert d["complex_s_disc_radius"] == 5
    assert d["conservative_Neumann_route_norm_upper"] == 18
    assert d["minimum_fermion_mass"] == 36


def test_derivative_operator_is_kept_but_constant_on_shell():
    d = low_degree.data()
    c = d["two_derivative_action_coefficient"]
    assert sp.simplify(d["on_shell_degree_two_four_vertex"] + 8 * c) == 0
    assert d["Euclidean_degree_two_four_vertex"].has(sp.Symbol("sum_external_squares"))
    assert "constant on the equal-mass shell" in d["scope"]


def test_local_threshold_is_counted_once():
    A0, AF, v4, k = sp.symbols("A0 AF v4 k")
    correct = (A0 - v4) / k**2 + (AF + v4) / k**2
    wrong = (A0 - v4) / k**2 + AF / k**2
    assert sp.factor(correct - (A0 + AF) / k**2) == 0
    assert sp.factor(wrong - correct) == -v4 / k**2
    assert "not" in reference.data()["scope"]


def test_large_valid_error_is_inconclusive_not_false_positive():
    d = calibration.canonical_band(1, 2, 0, Fraction(1, 2))
    assert d["lower"] == -4
    assert d["upper"] == 12
    assert d["positive_selected_coefficient"] is False
    assert d["strict_selected_increase_over_tree"] is False


def test_actual_selected_band_and_full_model_nonclosure():
    d = calibration.data()
    T = d["same_reference_tree_second_coefficient"]
    band = d["selected_exact_normalized_coefficient_interval"]
    assert T < band["lower"] < band["upper"] < T * (1 + sp.Rational(1, 10**204))
    assert d["rational_box_second_coefficient_relative_upper"] < sp.Rational(1, 10**604)
    assert d["formal_one_loop_fermion_increment_interval"][0] > 0
    assert "not the full new-model amplitude" in d["scope"]
    assert "not a bound on later primitive loop orders" in series.data()["scope"]


def test_exact_counts_and_scope():
    assert len(audit.residuals()) == 127
    assert audit.scalar_count() == 254
    assert len(audit.gates()) == 29
    assert len(audit.controls()) == 9
    assert audit.rejected_inputs() == 43
    assert audit.controls()["original_P8_not_closed"] is True
