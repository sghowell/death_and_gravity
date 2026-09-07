from fractions import Fraction as Q

import pytest
import sympy as sp
from p8_composite_vacuum import branch, cayley, independent, matching, vacuum


@pytest.mark.parametrize("module", (cayley, vacuum, branch, matching))
def test_exact_main_identities(module):
    assert all(sp.simplify(value) == 0 for value in module.checks().values())


def test_conserved_physical_source_projector():
    assert all(value == 0 for value in vacuum.source_projector_checks().values())


def test_separate_fraction_determinant_algebra():
    out = independent.checks()
    assert len(out["coefficientwise_identities"]) == 4
    assert out["Fraction_physical_matching_fixtures"]["physical_mass_squared_over_m_squared"] == "1/4"


def test_cayley_map_on_a_regular_exact_lorentzian_diagonal_chart():
    physical = sp.diag(1, -2, -3, -4)
    relative = sp.diag(sp.Rational(1, 3), sp.Rational(-1, 4), sp.Rational(1, 5), 0)
    out = cayley.matrix_map(physical, relative)
    assert out["g"]*out["S"]**2 == out["f"]
    assert out["g"]*(sp.eye(4)+out["S"])**2 == physical
    assert (sp.eye(4)-out["S"])*(sp.eye(4)+out["S"]).inv() == relative
    assert all(value > 0 for value in out["S"].diagonal())


def test_exchange_asymmetry_and_chart_controls_are_not_silently_accepted():
    controls = cayley.controls()
    assert controls["asymmetric_beta0_linear_relative_tadpole"] == -cayley.M2*cayley.m2/16
    assert controls["singular_Cayley_plus_example"] == 0
    with pytest.raises(ValueError):
        cayley.matrix_map(sp.eye(3), sp.eye(3))


def test_constant_vacuum_shift_and_wrong_mass_controls():
    d = vacuum.derive()
    c = vacuum.controls()
    assert d["shifted_beta_values"] == (-sp.Rational(3, 8), -sp.Rational(3, 8), sp.Rational(5, 8), -sp.Rational(3, 8), -sp.Rational(3, 8))
    assert c["lost_vacuum_shift_mass_error"] == 3*d["m2"]/4
    assert d["physical_massive_source_response"] == 0
    assert d["relative_auxiliary_probe_residue"] > 0
    assert d["physical_FP_mass_squared"] > 0


def test_same_potential_is_not_supplied_by_constant_vacuum_choice():
    c = vacuum.controls()
    eps = next(symbol for symbol in c["constant_vacuum_potential_mismatch_cleared"].free_symbols if symbol.name == "epsilon")
    assert sp.Poly(c["constant_vacuum_potential_mismatch_cleared"], eps).all_coeffs() == [36, 72, 39, -2, 3]
    for value in (sp.Rational(1, 100000), sp.Rational(1, 10000)):
        assert c["pinned_bounce_initial_potential_over_M2m2"].subs(eps, value)+sp.Rational(3, 8) >= sp.Rational(3749, 10000)


def test_differential_boundary_invertibility_is_not_algebraic_mass_positivity():
    c = branch.controls()
    assert c["positive_mass_Dirichlet_resonant_equation"] == 0
    assert c["positive_mass_Dirichlet_left"] == c["positive_mass_Dirichlet_right"] == 0
    assert c["positive_mass_Dirichlet_nonzero_solution"] == 1
    assert c["branch_switch_loses_linear_inverse"] == 0


def test_heavy_state_source_and_loop_controls():
    c = branch.controls()
    assert c["nonzero_heavy_initial_state"] == 1
    assert c["broken_symmetry_linear_source"] == 1
    assert c["even_heavy_quadratic_coefficient_can_produce_light_loop_dependence"] != 0


def test_no_tree_heavy_graph_but_nonzero_loop_graphs_are_allowed():
    assert branch.closed_heavy_graph_loop_bound([2]) == 1
    assert branch.closed_heavy_graph_loop_bound([2, 2, 2]) == 1
    assert branch.closed_heavy_graph_loop_bound([2, 4]) == 2
    for invalid in ([], [1, 1], [True, 2], [3, 3]):
        with pytest.raises(ValueError):
            branch.closed_heavy_graph_loop_bound(invalid)


def test_conditional_functional_residual_and_ball_arithmetic():
    assert branch.residual_bound(2, Q(1, 4), Q(3, 100)) == Q(2, 25)
    assert branch.ball_gate(2, Q(1, 4), Q(3, 100), Q(1, 10))
    assert not branch.ball_gate(2, Q(1, 4), Q(3, 100), Q(1, 20))
    for invalid in (True, 0.1, float("inf"), float("nan"), sp.Float("0.1")):
        with pytest.raises(TypeError):
            branch.residual_bound(1, 0, invalid)
    for args in ((0, 0, 0), (1, 1, 0), (1, Q(-1, 4), 0), (1, 0, -1)):
        with pytest.raises(ValueError):
            branch.residual_bound(*args)


def test_original_M1_operator_and_physical_metric_budgets():
    d = matching.derive()
    assert d["target_F2X_at_centre"] == -sp.Rational(1, 2)
    assert d["target_A3_at_centre"] == 1
    assert d["endpoint_H_error_threshold"] == sp.Rational(8, 5)
    assert d["exact_null_residual_lower_magnitude"] == sp.Rational(801, 100)
    approximate = d["approximate_null_residual_lower_magnitude"].subs(
        {d["epsilon_Hdot"]: sp.Rational(1, 10), d["epsilon_chidot"]: sp.Rational(1, 100)})
    assert approximate == sp.Rational(78081, 10000)


def test_static_endpoint_control_does_not_claim_old_M1_rolls_in_vacuum():
    d = matching.derive()
    assert sp.simplify(d["tau"]*d["target_H_right"]) == d["endpoint_H_error_threshold"]
    assert matching.controls()["dropping_free_chi_understates_null_residual_budget"] == Q(1, 100)


def test_actual_null_residual_budget_has_an_exact_domain_guard():
    assert matching.null_residual_budget(0, 0) == Q(801, 100)
    assert matching.null_residual_budget(Q(1, 10), Q(1, 100)) == Q(78081, 10000)
    for errors in ((4, 0), (-1, 0), (0, Q(1, 10)), (0, Q(-1, 100))):
        with pytest.raises(ValueError):
            matching.null_residual_budget(*errors)
    for inexact in (True, 0.01, sp.Float("0.01"), float("inf")):
        with pytest.raises(TypeError):
            matching.null_residual_budget(0, inexact)


def test_unguarded_chi_error_would_claim_a_false_square_lower_bound():
    # At error 1/5, an allowed actual chidot=0 falsifies the unguarded
    # claim chidot² >= (1/10-1/5)². It is outside the public bound domain.
    epsilon = Q(1, 5)
    actual_chidot = Q(0)
    assert abs(actual_chidot-Q(1, 10)) <= epsilon
    assert actual_chidot**2 < (Q(1, 10)-epsilon)**2
