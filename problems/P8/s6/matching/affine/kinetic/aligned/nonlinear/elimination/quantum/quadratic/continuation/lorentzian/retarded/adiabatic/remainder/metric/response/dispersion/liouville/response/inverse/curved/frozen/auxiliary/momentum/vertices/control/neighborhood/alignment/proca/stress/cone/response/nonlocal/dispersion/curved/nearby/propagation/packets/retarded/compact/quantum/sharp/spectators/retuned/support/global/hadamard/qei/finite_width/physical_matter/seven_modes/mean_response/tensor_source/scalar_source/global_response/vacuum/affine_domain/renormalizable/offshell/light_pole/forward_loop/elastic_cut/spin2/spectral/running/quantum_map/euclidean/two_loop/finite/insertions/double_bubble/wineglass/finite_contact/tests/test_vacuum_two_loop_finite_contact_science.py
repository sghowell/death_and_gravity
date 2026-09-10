"""Fixed-potential, radial, insertion and actual-error checks."""

from fractions import Fraction

import pytest
import sympy as sp
from p8_vacuum_two_loop_finite_contact import (
    audit,
    calibration,
    insertion,
    potential,
    radial,
    subtraction,
)


@pytest.mark.parametrize(
    "name,value", list(audit.residuals().items()), ids=list(audit.residuals())
)
def test_exact_Hessian_integral_insertion_or_reference_identity(name, value):
    assert value == 0, name


@pytest.mark.parametrize(
    "name,value", list(audit.gates().items()), ids=list(audit.gates())
)
def test_fixed_condition_analytic_bound_or_original_scope(name, value):
    assert bool(value), name


@pytest.mark.parametrize(
    "name,call,args",
    calibration.bad_cases(),
    ids=[row[0] for row in calibration.bad_cases()],
)
def test_reject_unsupported_radial_input(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize(
    "name,value",
    [(k, v) for k, v in audit.controls().items() if k != "rejected_inputs"],
    ids=[k for k in audit.controls() if k != "rejected_inputs"],
)
def test_negative_control_or_scope_boundary(name, value):
    assert bool(value), name


@pytest.mark.parametrize("y", (0, Fraction(1, 2), 1, 4, 10**200, 10**400))
def test_exact_all_radius_negative_contact_kernel(y):
    d = calibration.point(y)
    assert d["radial_invariant"] == sp.Rational(y)
    assert 0 < d["actual_positive_F"] < d["actual_positive_A"]
    assert d["negative_finite_contact_kernel"] < 0
    density = d["finite_contact_density_without_loop_measure"]
    assert (
        density
        == 6 * sp.Rational(y) * d["negative_finite_contact_kernel"] / (y + 1) ** 2
    )
    assert density == 0 if y == 0 else density < 0


@pytest.mark.parametrize("index", (1, 2))
def test_each_radial_primitive_independently_anchored(index):
    d = radial.data()
    prefix = "one" if index == 1 else "two"
    P = d[prefix + "_heavy_factor_primitive"]
    f = d[prefix + "_heavy_factor_integrand"]
    y = d["y"]
    assert sp.factor(sp.diff(P, y) - f) == 0
    assert sp.limit(P, y, sp.oo) == 0
    assert sp.simplify(-P.subs(y, 0) - d["J" + str(index)]) == 0


def test_fixed_contact_is_concrete_not_an_adjustable_forward_parameter():
    d = calibration.data()
    assert not d["actual_fixed_finite_potential_contact"].free_symbols
    assert d["actual_fixed_finite_potential_contact"].has(
        sp.log(d["actual_heavy_mass_squared"])
    )
    assert "already fixed" in potential.data()["scope"]
    assert "zero second momentum derivative" in insertion.data()["scope"]


def test_actual_complete_insertion_and_known_error_budget():
    d = calibration.data()
    L, tree = d["actual_quartic"], d["actual_tree_b2"]
    assert d["actual_contact_absolute_upper"] == sp.Rational(50, 3) * L**2
    assert (
        d["actual_one_loop_quartic_derivative_b2_upper"] == sp.Rational(4812, 144) * L
    )
    error = d["actual_finite_contact_insertion_b2_absolute_upper"]
    assert error == sp.Rational(10025, 18) * L**3
    assert 0 < error < sp.Rational(7, 10**612)
    assert error / tree < sp.Rational(2, 10**12)
    two = d["raw_graph_and_finite_contact_two_loop_upper"]
    assert two == error + d["previous_four_raw_graph_groups_upper"]
    assert two / tree < sp.Rational(1, 10**7)
    assert d["known_one_and_two_loop_contributions_upper"] / tree < sp.Rational(
        1, 10**6
    )
    assert d["positive_tree_minus_known_corrections_lower"] > 0


def test_negative_contact_bound_does_not_remove_reference_counterterms():
    d = subtraction.data()
    ct = d["local_reference_continuation_counterterms"]
    assert ct["delta_polynomial_quartic"] != 0
    assert ct["delta_cubic_squared"] != 0
    assert ct["delta_heavy_mass_squared"] == 0
    assert all(v == 0 for v in d["checks"].values())


def test_partial_derivative_keeps_g_M_and_both_triangle_integrals():
    d = insertion.data()
    derivative = d["fixed_g_and_M_quartic_variation"]
    TA = sp.Symbol("first_integrated_heavy_triangle")
    TB = sp.Symbol("second_integrated_heavy_triangle")
    assert sp.diff(derivative, TA) == -sp.Rational(1, 2)
    assert sp.diff(derivative, TB) == -sp.Rational(1, 2)
    assert d["diagnostic_nonzero_local_contact_loop_coefficient"] != 0


def test_equal_inexact_inputs_stay_rejected_after_exact_point():
    calibration.point(1)
    for bad in (True, 1.0, sp.Float(1)):
        with pytest.raises((TypeError, ValueError)):
            calibration.point(bad)


def test_known_contributions_are_not_full_pole_LSZ_or_original_completion():
    assert "completed two-loop pole/residue/LSZ" in calibration.data()["scope"]
    assert audit.gates()["two_loop_source_aware_matching_not_computed"]
    assert audit.controls()["full_V_G_B_and_original_P8_not_closed"]


def test_exact_counts():
    assert len(audit.residuals()) == 40
    assert len(audit.gates()) == 34
    assert len(audit.controls()) == 9
    assert audit.rejected_inputs() == 14
