"""Exact selected-state cancellation, continuous window and front checks."""
import pytest
import sympy as sp
from p8_clock_tadpole import fronts, profiles, state, verify, window


def test_all_actual_profile_and_front_identities():
    verify.affine.certify_residuals(verify.residuals())
    assert len(verify.residuals()) == 30


def test_fixed_selected_profiles_cancel_all_three_first_variations():
    data = profiles.covariant()
    assert data["energy"] == -data["rho"]
    assert data["spatial_pressure"] == -data["pressure"]
    assert sp.simplify(data["clock_source"]-sp.diff(data["rho"], profiles.u)
                       -3*data["H"]*(data["rho"]+data["pressure"])) == 0
    another = sp.Symbol("another_state_energy", real=True)
    assert data["energy"]+another != 0


def test_actual_point_chart_lapse_square_at_the_bounce():
    data = profiles.point_chart()
    assert sp.factor(data["coefficients"]["N_square_coefficient"].subs(data["h"], 1)
                     -data["rho"]+sp.Rational(7, 8)*data["pressure"]) == 0
    assert all(value == 0 for value in data["checks"].values())


def test_retained_vector_principal_speed_uses_the_full_canonical_rate():
    data = fronts.principal()
    assert data["modes"]["transverse"]["squared_principal_speed"] == 1
    assert data["modes"]["longitudinal"]["squared_principal_speed"] == fronts.modes.bm/fronts.modes.am
    assert all(value == 0 for value in data["checks"].values())
    assert all(fronts.proof_checks().values())
    assert fronts.speed(sp.Rational(49, 100))["longitudinal_squared_principal_speed"] < 1
    assert fronts.speed(sp.Rational(1, 2))["longitudinal_squared_principal_speed"] == 1
    assert fronts.speed(sp.Rational(51, 100))["longitudinal_squared_principal_speed"] > 1


def test_window_is_flat_on_the_whole_original_tube_and_zero_near_x_zero():
    assert all(window.proof_checks().values())
    for value in (-1, -sp.Rational(9, 10), -sp.Rational(11, 10), -sp.Rational(7, 8)):
        assert window.cutoff(value) == 1
    for value in (0, -sp.Rational(1, 2), -sp.Rational(5, 4)):
        assert window.cutoff(value) == 0
    assert window.cutoff(-1+sp.sqrt(sp.Rational(3, 128))) == sp.Rational(1, 2)


def test_continuous_window_derivative_bounds_and_exact_chain_rule():
    assert [window.derivative_bound(j) for j in range(6)] == [
        1, 8704, 22874112, 92896886784, 509829585567744, 3527551905599324160]
    assert all(value == 0 for value in window.chain_rule_checks().values())


def test_full_mixed_profile_jet_budget_and_literal_lapse_coefficient():
    data = profiles.profile_bounds(10**24, 1000)
    mixed = data["all_clock_norm_mixed_derivative_bounds_total_order_at_most_five"]
    assert len(mixed) == 21
    assert all(0 < value < sp.Rational(1, 10**18) for value in mixed.values())
    assert max(mixed.values()) < sp.Rational(1456, 10**22)
    assert sp.Rational(3, 8)-data["lapse_square_coefficient_absolute_upper"] > sp.Rational(1, 3)
    assert all(row["P_xx"] == 0 for row in data["fixed_profile_derivative_bounds_on_clock_tube"].values())
    assert data["profile_held_fixed_under_subsequent_state_and_metric_variations"]
    assert not data["full_quantum_or_UV_completion_claim"]


def test_actual_frozen_initial_state_is_not_replaced_by_the_reference():
    for data in state.initial_data_examples().values():
        assert data["active_higher_derivative_orders"] == [6]
        assert data["included_correction_terms"][6]["weight"] == sp.Rational(1, 2)
        assert data["all_order_frequency"] > 0
        assert data["all_order_frequency"] != data["frozen_frequency"]


def test_explicit_profile_definition_has_the_physical_scale_normalization():
    data = state.definition(10**24, 1000)
    assert data["finite_integral_prefactor"] == sp.Rational(1, 10**48)
    assert data["initial_u"] == -sp.Rational(1, 2)
    assert set(data["local_finite_profile_terms"]) == {"energy", "pressure"}


@pytest.mark.parametrize("order", (True, 1.0, sp.Float(1), sp.Integer(1)))
def test_native_derivative_order_validation_after_cache_warmup(order):
    window.derivative_bound(1)
    with pytest.raises(ValueError, match="native"):
        window.derivative_bound(order)


def test_all_continuous_proof_checks_and_input_rejections():
    assert all(verify.proof_checks().values())
    assert verify.controls()["rejected_inputs"] == 110
