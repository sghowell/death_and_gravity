"""Homogeneous action, actual source, continuous response and chart lift."""
import sympy as sp
from p8_margin_response import bounds, system


def assert_zero(values):
    for value in values.values():
        assert all(entry == 0 for entry in value) if isinstance(value, sp.MatrixBase) else value == 0


def test_literal_homogeneous_action_and_physical_metric_source():
    assert_zero(system.checks())


def test_continuous_derivative_and_chart_reconstructions():
    assert_zero(bounds.checks())


def test_all_continuous_response_and_geometry_margins():
    assert all(value is True for value in bounds.proof_checks().values())


def test_lapse_is_not_an_independent_zero_initial_datum():
    data = system.reduction()
    initial = data["lapse"].subs({system.old.v: 0, system.Y: 0})
    assert sp.factor(initial+system.Fn/(2*system.Jnew)) == 0


def test_free_matter_momentum_is_preserved_jointly_not_separately():
    assert system.reduction()["preserved_physical_matter_momentum"] == 0
    assert sp.diff(system.reduction()["matter_dot"], system.Fn) != 0


def test_zero_stress_has_zero_response_budget():
    data = bounds.response_bounds(0, 0)
    assert all(value == 0 for name, value in data.items() if name.endswith("_upper"))


def test_actual_vector_state_is_inside_the_rounded_controlled_example():
    data = bounds.vector_example(10**12, 1000)["response"]
    assert data["eta_value"] < sp.Rational(1, 10**14)
    assert data["eta_time_derivative"] < sp.Rational(1, 10**14)
    assert data["full_self_consistent_quantum_solution_or_nonlinear_residual_claim"] is False


def test_exact_chart_Hubble_bound_is_not_a_full_semiclassical_solution():
    d = bounds.response_bounds(sp.Rational(1, 10**14), sp.Rational(1, 10**14))
    assert d["physical_fractional_scale_correction_upper"] == sp.Rational(9, 2*10**12)
    assert d["exact_chart_lift_Hubble_difference_over_inverse_tau_upper"] < sp.Rational(1, 10**9)
