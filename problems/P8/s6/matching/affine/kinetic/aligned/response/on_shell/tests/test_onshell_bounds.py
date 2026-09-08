"""Continuous compact bounds, exact scaling and physical-response transfer."""
import pytest
import sympy as sp
from p8_aligned_onshell import bounds as b
from p8_aligned_onshell import phase


def test_continuous_proofs_and_rational_growth_majorant():
    assert all(value is True for value in b.proof_checks().values())
    data = b.continuous()
    assert data["growth_exponent"] == sp.Rational(4771877971625, 20115881984)
    assert data["growth_bound"] == 3**238
    assert data["normalization_C"] > data["B_N"]+data["B_K"]


@pytest.mark.parametrize("point", (phase.LEFT, -sp.Rational(1, 4), 0, sp.Rational(1, 4), phase.RIGHT))
def test_independent_dense_rational_bounds(point):
    data, limits = phase.system(), b.continuous()
    for key, bound_key, orders in (("A", "A_time_jet_bounds", 3),
                                   ("lapse_row", "lapse_observer_time_jet_bounds", 4),
                                   ("trace_row", "trace_observer_time_jet_bounds", 4)):
        for order in range(orders):
            actual = data[key].diff(phase.u, order).subs(phase.u, point)
            assert actual.norm(sp.oo) <= limits[bound_key][order]


def test_explicit_nonzero_initial_scaling_satisfies_both_jet_bounds():
    epsilon = sp.Rational(1, 1000)
    family, limits = b.prepared_family(epsilon), b.continuous()
    assert family["initial_scale"] > 0
    assert family["initial_scale"]*limits["growth_bound"]*limits["normalization_C"] == epsilon
    assert family["initial_state"].norm(sp.oo) == family["initial_scale"]
    assert 0 < family["lapse_jet_upper"] < epsilon
    assert 0 < family["trace_jet_upper"] < epsilon
    assert family["nonzero_source_third_cosine_squared_coefficient"] != 0
    assert family["nonlinear_solution_or_higher_order_remainder_claim"] is False


def test_actual_on_shell_source_transfers_to_the_prepared_response():
    result = b.transfer(sp.Rational(1, 1000), sp.Rational(1, 10**9))
    assert result["family"]["source_jet_E"] == sp.Rational(136139, 8000000)
    response = result["response"]
    assert response["k_com_squared"] == sp.Rational(1, 4)
    assert response["source_jet_E"] == result["family"]["source_jet_E"]
    assert response["spatial_readout_error_against_zeta_sqrt_q_times_Sprime_plus_2rhoS"] == sp.Rational(408417, 4000000000000000)
    assert response["temporal_readout_error_against_S"] == sp.Rational(136139, 16000000000000)
    assert result["source_class_realized_by_nonzero_linear_light_solution"] is True
    assert result["higher_order_nonlinear_or_quantum_error_claim"] is False


@pytest.mark.parametrize("value", (0, -1, 0.001, True, sp.oo, sp.Symbol("epsilon", positive=True)))
def test_nontrivial_family_requires_exact_positive_amplitude(value):
    with pytest.raises((TypeError, ValueError)):
        b.prepared_family(value)
