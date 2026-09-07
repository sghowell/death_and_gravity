import pytest
import sympy as sp
from p8_hr_tree import components, controls, edge, graph, independent


def test_literal_edge_action_and_both_reciprocal_weights():
    assert all(sp.simplify(value) == 0 for value in edge.checks().values())
    assert all(sp.simplify(value) == 0 for value in edge.vertex_checks().values())
    assert all(sp.simplify(value) != 0 for value in edge.controls().values())


@pytest.mark.parametrize("size", (1, 2, 3, 5))
def test_conditional_component_algebra_includes_nonnegative_Einstein_weights(size):
    d = components.derive(size)
    assert d["A"].is_positive is True
    assert all(sp.simplify(value) == 0 for value in components.checks(size).values())


def test_independent_polynomial_and_jet_replay():
    report = independent.checks()
    assert len(report["coefficientwise_identities"]) == 6
    assert len(report["literal_edge_stress_fixtures"]) == 5
    assert len(report["Gaussian_incidence_fixtures"]) == 4
    assert len(report["component_first_jet_fixtures"]) == 5


def test_directed_tree_fluxes_are_zero_only_for_zero_divergence():
    assert len(graph.checks()) == 5
    assert graph.solve_fluxes(3, [(0, 1), (2, 1)], [1, -1, 0]) == (1, 0)


@pytest.mark.parametrize("count,edges", [(3, [(0, 1), (1, 2), (2, 0)]),
                                        (3, [(0, 1)]), (2, [(0, 1), (1, 0)]),
                                        (2, [(0, 0)]), (2, [(0, 2)]), (True, []),
                                        (2, [(False, 1)]), (0, [])])
def test_cycle_disconnected_and_invalid_graphs_are_not_silently_admitted(count, edges):
    with pytest.raises((TypeError, ValueError)):
        graph.validate_tree(count, edges)


def test_algebraic_zero_edge_is_deleted_even_at_an_isolated_or_double_root():
    edges = [(1, 0), (2, 0)]
    assert graph.active_component(3, edges, [0, 2]) == (0, 2)
    assert graph.active_component(3, edges, [sp.Rational(1, 1000), -2]) == (0, 1, 2)
    assert graph.active_component(3, edges, [0, 0]) == (0,)


def test_auxiliary_zero_EH_vertex_is_allowed_but_target_zero_EH_is_not():
    d = components.reconstruct([1, 0, 2], [1, 2, 3], [1, 3, 2], [1, 0, 2], 1, [0, 1, 2])
    assert d["A"] == 19
    assert abs(d["lambda"]) <= d["component_log_rate_bound"]
    with pytest.raises(ValueError):
        components.reconstruct([0, 1], [1, 1], [1, 1], [0, 0], 0, [0])


@pytest.mark.parametrize("bad", (True, 0.1, sp.Float("0.1"), sp.oo, sp.I))
def test_certificate_arithmetic_rejects_float_bool_and_nonfinite_inputs(bad):
    with pytest.raises((TypeError, ValueError)):
        components.reconstruct([1], [1], [1], [0], bad, [0])


def test_all_actual_solution_controls_satisfy_full_background_and_scalar_equations():
    assert all(sp.simplify(value) == 0 for value in controls.checks().values())
    assert controls.zero_target_exception()["H_root_prime_at_zero"] == 2
    assert controls.nec_violating_control()["null_at_bounce"] == -4


def test_actual_mixed_algebraic_tree_rejects_including_the_inactive_leaf():
    d = controls.mixed_algebraic_tree()
    assert d["correct_component_null_residual"] == 0
    assert d["false_all_vertex_null_residual"] == 2/(3*d["T"]**2)
    assert graph.active_component(3, d["edges"], d["P_values"]) == d["component"]


def test_tuned_edge_and_separate_stiff_fields_allow_nontrivial_regular_expansion():
    d = controls.stiff_tree()
    assert d["G"] == (1, 0, 2, 3)
    assert all(sp.simplify(h-1/(3*d["T"])) == 0 for h in d["H"])
    assert all(value == 0 for value in d["interaction_rho"]+d["interaction_pressure"])


def test_strict_actual_vertex_CD_endpoint_bound():
    assert components.cd_endpoint_test(sp.Rational(1, 2), 0)["necessary_normalized_endpoint_error"] == sp.Rational(8, 5)
    assert components.cd_endpoint_test(sp.Rational(1, 2), sp.Rational(8, 5))["status"] == "INCONCLUSIVE"
