from fractions import Fraction as Q

import pytest
import sympy as sp
from p8_trimetric_cones import background, bounds, independent, rolling, tensor


@pytest.mark.parametrize("module", (background, tensor, rolling, bounds))
def test_exact_core_identities(module):
    assert all(sp.simplify(value) == 0 for value in module.checks().values())


def test_independent_polynomial_and_full_rolling_checks():
    d = independent.checks()
    assert len(d["coefficientwise_identities"]) == 4
    assert d["literal_source_lapse_scale_jet_directions"] == 2
    assert len(d["full_rolling_Fraction_fixtures"]) == 4
    assert len(d["independent_physical_clock_pullbacks"]) == 2
    for fixture in independent.clock_fixtures():
        primary = tensor.derive()
        point = {primary[key]: sp.Rational(fixture[key].numerator, fixture[key].denominator)
                 for key in ("A", "N", "G")}
        for key in ("G_T_h", "F_T_h", "common_speed_squared"):
            assert primary[key].subs(point) == sp.Rational(fixture[key].numerator, fixture[key].denominator)


def test_positive_weighted_cones_cannot_all_be_at_most_one_with_rolling_matter():
    # Denominator-free independent finite fixtures of the written inequality.
    for pe, pf in ((Q(1), Q(2)), (Q(3, 7), Q(5, 11))):
        for ce, cf in ((Q(1), Q(1)), (Q(1, 2), Q(3, 4))):
            assert pe*(ce-1)+pf*(cf-1) <= 0
            assert pe*(ce-1)+pf*(cf-1) != Q(1, 100)


def test_mixed_sign_auxiliary_point_has_no_full_solution_claim():
    c = background.controls()
    assert c["mixed_sign_auxiliary_lapse"] == c["mixed_sign_auxiliary_space"] == 0
    assert (c["mixed_sign_auxiliary_c_e"], c["mixed_sign_auxiliary_c_f"]) == (sp.Rational(1, 2), 1)
    assert c["mixed_sign_actual_null_stress"] == 2


def test_mixed_sign_full_flat_countercontrol_prevents_positive_link_inference():
    c = background.controls()
    assert c["mixed_sign_full_flat_e"] == c["mixed_sign_full_flat_v"] == c["mixed_sign_full_flat_u"] == 0
    assert c["mixed_sign_full_flat_relative_spring"] == 4


def test_auxiliary_response_keeps_external_contact_separate_from_propagation():
    d = tensor.derive()
    pg, pf, j, ge, gf = (d[key] for key in ("P_g", "P_f", "j", "gamma_e", "gamma_v"))
    point = {pg: 2, pf: 3, j: 5, ge: 7, gf: 11}
    assert d["sourced_auxiliary_map"].subs(point) == sp.Rational(57, 5)
    assert d["source_free_auxiliary_map"].subs(point) == sp.Rational(47, 5)
    assert d["source_contact"].subs(point) == 5


def test_source_free_relative_mode_does_not_change_physical_common_tensor():
    d = tensor.derive()
    sub = {d["P_g"]: 1, d["P_f"]: 1, d["gamma_e"]: d["delta"], d["gamma_v"]: -d["delta"]}
    assert sp.simplify(d["source_free_auxiliary_map"].subs(sub)) == 0
    assert sp.diff(d["exchange_even_odd_density"], d["gamma_dot"], d["delta_dot"]) == 0
    assert sp.diff(d["exchange_even_odd_density"], d["gamma"], d["delta"]) == 0


def test_actual_clock_and_kinetic_coefficients_at_a_real_rolling_point():
    d = tensor.derive()
    point = {d["A"]: 2, d["N"]: sp.Rational(2, 7), d["G"]: 1}
    assert d["G_T_h"].subs(point) == sp.Rational(1, 14)
    assert d["F_T_h"].subs(point) == sp.Rational(7, 2)
    assert d["common_speed_squared"].subs(point) == 49
    assert d["relative_algebraic_mass_h_squared"].subs({**point, d["q"]: 1}) == 28


def test_actual_solution_not_just_an_auxiliary_algebraic_point():
    d = rolling.derive()
    for key in ("A", "N", "rho", "H_r_squared", "D_h_psi"):
        assert d[key].is_positive
    initial = rolling.initial_data()
    assert initial["A"] == 2 and initial["N"] == sp.Rational(2, 7)
    assert initial["physical_speed"] == 7
    assert sp.simplify(d["H_h"]).is_positive


def test_source_strength_does_not_suppress_the_gap_at_fixed_canonical_stress():
    a, q, nh = sp.symbols("A q n_h", positive=True)
    result = bounds.speed_budget(a, q, nh, 1)
    assert result["speed_gap"] == a*nh/(4*q)
    assert "epsilon" not in str(result)


def test_remainder_budget_is_a_conditional_threshold_not_an_error_estimate():
    c = bounds.controls()
    assert c["actual_example_F_minus_G"] == sp.Rational(24, 7)
    assert c["protected_kinetic_margin"] == sp.Rational(3, 56)
    assert c["protected_cone_margin"] == sp.Rational(95, 28)
    assert c["nonrolling_has_no_strict_gap"] == 0


@pytest.mark.parametrize("bad", (True, 0.1, float("inf"), sp.Float("0.1"), sp.oo, sp.nan, sp.I))
def test_certified_numeric_budget_rejects_inexact_or_nonfinite_inputs(bad):
    with pytest.raises(TypeError):
        bounds.speed_budget(1, 1, bad, 1)


@pytest.mark.parametrize("args", ((0, 1, 1, 1), (1, -1, 1, 1), (1, 1, -1, 1), (1, 1, 1, 0), (1, 1, 1, 1, -1)))
def test_certified_budget_rejects_false_domains(args):
    with pytest.raises(ValueError):
        bounds.speed_budget(*args)
