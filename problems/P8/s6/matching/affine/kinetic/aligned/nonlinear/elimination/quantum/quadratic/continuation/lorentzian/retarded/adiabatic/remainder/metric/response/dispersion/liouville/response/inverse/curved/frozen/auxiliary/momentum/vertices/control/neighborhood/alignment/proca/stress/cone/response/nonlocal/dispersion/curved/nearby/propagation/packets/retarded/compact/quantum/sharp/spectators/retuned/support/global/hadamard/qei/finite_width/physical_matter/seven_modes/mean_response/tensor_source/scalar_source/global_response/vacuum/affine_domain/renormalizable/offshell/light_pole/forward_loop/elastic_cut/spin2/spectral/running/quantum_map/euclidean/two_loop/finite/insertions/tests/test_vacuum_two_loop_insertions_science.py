"""Actual kernel, independent routing, subtraction and integrated bound tests."""

from fractions import Fraction

import pytest
import sympy as sp
from p8_vacuum_light_pole import kernel as original_kernel
from p8_vacuum_two_loop_insertions import (
    audit,
    calibration,
    kernel,
    radial,
    routing,
    subtraction,
)


@pytest.mark.parametrize(
    "name,value", list(audit.residuals().items()), ids=list(audit.residuals())
)
def test_exact_kernel_routing_subtraction_or_bound_identity(name, value):
    assert value == 0, name


@pytest.mark.parametrize(
    "name,value", list(audit.gates().items()), ids=list(audit.gates())
)
def test_analytic_domain_subtraction_or_scope_gate(name, value):
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
def test_negative_control_or_original_scope(name, value):
    assert bool(value), name


@pytest.mark.parametrize("channel", ("s", "t", "u"))
def test_actual_bilinear_routing_and_both_internal_light_lines(channel):
    d = routing.data()
    r = d["channels"][channel]
    P = r["pair_momentum"]
    assert sp.simplify(P.dot(P) + r["channel_invariant"]) == 0
    assert all(
        sp.expand(value) == 0
        for value in r["second_light_line"] - r["first_light_line"] - P
    )
    assert len(r["first_vertex_heavy_exchange_shifts"]) == 2
    assert d["centered_light_real_mass_gap"] == sp.Rational(1, 4)


@pytest.mark.parametrize("y", (0, Fraction(1, 2), 1, 4, 10**200, 10**400))
def test_exact_radial_enclosures_without_loop_cutoff(y):
    d = calibration.point(y)
    p = calibration.data()
    M, g = p["actual_heavy_mass_squared"], p["actual_cubic_squared"]
    assert d["radial_invariant"] == sp.Rational(y)
    expected = 1 / M if y == 0 else sp.log(1 + sp.Rational(y) / M) / y
    assert d["decaying_remainder_logarithmic_upper"] == g * expected / 72
    assert d["decaying_remainder_uniform_upper"] == g / (144 * M)
    assert d["combined_light_denominator_real_lower"] == y + sp.Rational(1, 4)


@pytest.mark.parametrize(
    "u", (sp.Rational(-1, 2), -sp.Rational(1, 4) + sp.I, 0, 2 + 3 * sp.I)
)
def test_elementary_complex_fraction_strip_bound(u):
    assert sp.expand_complex(abs(1 + u) ** 2 - abs(u) ** 2) == 1 + 2 * sp.re(u)
    assert abs(u) ** 2 <= abs(1 + u) ** 2


def test_kernel_is_the_actual_parent_not_a_new_mass_surrogate():
    old = original_kernel.data()
    d = kernel.data()
    assert d["heavy_mass_squared"] == old["heavy_mass_squared"]
    assert d["b"] == old["normalized_parameter_weight"]
    assert d["prefactor"] == old["self_energy_prefactor"]
    assert d["fixed_asymptotic_multiplier"] == d["prefactor"] * sp.Integral(
        d["b"], (d["x"], 0, 1)
    )
    assert "inherited outer subtraction" in d["scope"]


def test_integrated_radial_majorant_is_anchored_at_both_endpoints():
    d = radial.data()
    primitive, y = d["anchored_radial_primitive"], d["y"]
    assert sp.factor(sp.diff(primitive, y) - d["complete_radial_integrand"]) == 0
    assert sp.limit(primitive, y, sp.oo) == 0
    assert sp.simplify(-primitive.subs(y, 0) - d["full_zero_to_infinity_integral"]) == 0
    assert d["unit_disc_Cauchy_coefficient_factor"] == 1


def test_exact_family_weight_both_insertions_and_all_sixty_four_refinements():
    d = subtraction.data()
    assert d["grouped_tadpole_family_Wick_weight"] == sp.Rational(3, 2)
    assert (
        d["outer_choice_count"] * d["inner_choice_count"]
        == d["grouped_refinement_count"]
        == 64
    )
    assert len(d["total_outer_local_counterterms_for_this_group"]) == 3
    assert all(v == 0 for v in d["checks"].values())


def test_actual_error_keeps_nonzero_constant_term_and_previous_finite_group():
    d = calibration.data()
    error = d["complete_grouped_insertion_b2_absolute_upper"]
    const = d["constant_insertion_group_b2_upper"]
    rest = d["decaying_insertion_group_b2_upper"]
    assert error == const + rest
    assert 0 < const < sp.Rational(1, 10**813)
    assert 0 < error < sp.Rational(1, 10**614)
    assert error / d["actual_tree_b2"] < sp.Rational(2, 10**15)
    combined = d["two_integrated_groups_combined_upper"]
    assert combined == error + d["previous_UV_finite_group_b2_upper"]
    assert combined < sp.Rational(3, 10**607)
    assert combined / d["actual_tree_b2"] < sp.Rational(1, 10**7)
    assert d["raw_refinements_represented_in_both_groups"] == 152
    assert d["remaining_subtraction_dependent_raw_refinements"] == 40


def test_equal_inexact_inputs_remain_rejected_after_exact_api_call():
    calibration.point(1)
    for bad in (True, 1.0, sp.Float(1)):
        with pytest.raises((TypeError, ValueError)):
            calibration.point(bad)


def test_partial_group_is_not_a_complete_two_loop_or_P8_result():
    assert (
        "complete two-loop pole/LSZ normalization remain open"
        in calibration.data()["scope"]
    )
    assert audit.gates()["forty_other_subtraction_dependent_refinements_still_open"]
    assert audit.controls()["full_two_loop_and_original_P8_not_closed"]


def test_exact_counts():
    assert len(audit.residuals()) == 51
    assert len(audit.gates()) == 34
    assert len(audit.controls()) == 9
    assert audit.rejected_inputs() == 14
