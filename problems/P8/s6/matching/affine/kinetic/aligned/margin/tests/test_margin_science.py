"""Literal action, both physical principal charts and finite-jet controls."""
import sympy as sp
from p8_aligned_margin import action, bounds, dynamics


def assert_zero(values):
    for value in values.values():
        assert all(entry == 0 for entry in value) if isinstance(value, sp.MatrixBase) else value == 0


def test_lower_scalar_preserves_background_but_changes_second_variation():
    assert_zero(action.checks())
    assert action.deformation()["clock_x_second"].is_positive is True


def test_temporal_and_shift_constraints_and_regular_Hamiltonian():
    data = dynamics.scalar()
    assert_zero({name: data[name] for name in ("shift_constraint", "temporal_constraint_unchanged",
                "unitary_kinetic_change", "regular_lapse_equation", "regular_Hamiltonian")})
    assert data["no_Theta_inverse"] is True


def test_regular_gamma_principal_keeps_moving_time_boundary():
    data = dynamics.gamma()
    assert_zero({name: data[name] for name in ("finite_q_kinetic_identity",
                "beta_leading_symmetric", "gamma_has_no_q_squared_principal",
                "principal_time_boundary_retained")})


def test_both_charts_have_same_positive_subluminal_clock_and_luminal_matter():
    assert_zero(dynamics.cone_checks())


def test_continuous_global_and_compact_polynomial_bounds():
    assert_zero(bounds.checks())
    assert all(value is True for value in bounds.proof_checks().values())


def test_center_threshold_is_a_velocity_chart_not_a_phase_pole():
    data = dynamics.center()
    eps, q = action.epsilon, dynamics.old.q
    K = data["finite_q_kinetic"]
    assert sp.factor(K.det()-q*(sp.Rational(1199, 800)+4*eps)/(2*(q/4-sp.Rational(3, 2)-4*eps))) == 0
    assert data["velocity_chart_threshold"] == 6+16*eps
    assert data["center_gradient_identity"] == sp.zeros(2)


def test_negative_deformation_is_a_superluminal_control_not_an_allowed_model():
    clock_speed = dynamics.center()["principal_squared_speeds"][1]
    assert clock_speed.subs(action.epsilon, 0) == 1
    assert clock_speed.subs(action.epsilon, -sp.Rational(1, 10**6)) > 1


def test_omitted_time_boundary_would_erase_the_center_clock_gradient():
    old = dynamics.old
    without = -dynamics.gamma()["gamma_leading"]/2
    without = without.subs({old.H: 0, old.theta: 0, old.lam: -sp.Rational(1, 2)})
    actual = dynamics.center()["independently_differentiated_principal_gradient"]
    assert without[0, 0] == 0
    assert actual[0, 0] == 6


def test_example_margin_is_explicit_and_does_not_assert_quantum_stability():
    data = bounds.physical_bounds(sp.Rational(1, 10**6), 10**12, 1000)
    assert data["compact_clock_cone_margin_lower"] > sp.Rational(1, 10**7)
    assert data["isolated_potential_jet_is_dominated"] is True
    assert data["full_quantum_cone_or_cutoff_claim"] is False
