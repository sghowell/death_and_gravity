"""All-order recurrence, Cauchy preparation and constrained-state controls."""
import pytest
import sympy as sp
from p8_vector_hadamard import cutoffs, preparation, proca, series, transfer


def assert_zero(values):
    for value in values.values():
        assert all(entry == 0 for entry in value) if isinstance(value, sp.MatrixBase) else value == 0


def test_recurrence_preserves_frozen_coefficients_and_cancels_sixth_order_residual():
    assert_zero(series.low_order_checks())


def test_all_order_recurrence_valuation_and_smooth_cutoff_controls():
    assert_zero(cutoffs.structural_checks())
    assert_zero(cutoffs.checks())
    assert cutoffs.turn_on(sp.Integer(1)) == 0
    assert cutoffs.turn_on(sp.Integer(2)) == 1
    assert cutoffs.turn_on(sp.Rational(3, 2)) == sp.Rational(1, 2)


def test_positive_Cauchy_frequency_and_bounded_initial_mixing():
    assert all(value is True for value in cutoffs.proof_checks().values())


def test_recursive_cutoffs_are_derived_from_both_value_and_slope():
    for kind in ("transverse", "longitudinal"):
        assert cutoffs.threshold(kind, 3, 1000) == 2000
        bounds = series.coefficient_bounds(kind, 4)
        assert bounds["coefficient_reconstruction"] == 0
        assert bounds["slope_reconstruction"] == 0
        threshold = cutoffs.threshold(kind, 4, 1000)
        assert threshold >= 2*cutoffs.threshold(kind, 3, 1000)
        assert threshold >= 16*(1+max(bounds["coefficient_upper"], bounds["frequency_coefficient_slope_upper"]))


def test_warm_caches_preserve_native_integer_validation():
    for function, order in ((series.coefficient, 1), (series.coefficient_bounds, 1)):
        function("transverse", order)
        for invalid in (True, 1.0, sp.Integer(1), sp.Float(1)):
            with pytest.raises(ValueError, match="native integer"):
                function("transverse", invalid)
    cutoffs.threshold("transverse", 4, 1000)
    for invalid in (4.0, sp.Integer(4), sp.Float(4)):
        with pytest.raises(ValueError, match="order>=3"):
            cutoffs.threshold("transverse", invalid, 1000)


def test_low_momentum_data_remain_exactly_the_frozen_preparation():
    for kind in ("transverse", "longitudinal"):
        for nu in (1000, 2000):
            data = preparation.initial_data(kind, 1000, nu)
            assert data["active_higher_derivative_orders"] == []
            assert data["all_order_frequency"] == data["frozen_frequency"]
            assert data["all_order_frequency_slope"] == data["frozen_frequency_slope"]


def test_locally_finite_Cauchy_preparation_at_half_and_full_cutoff():
    for kind in ("transverse", "longitudinal"):
        for nu, weight in ((3000, sp.Rational(1, 2)), (4000, 1)):
            data = preparation.initial_data(kind, 1000, nu)
            assert data["active_higher_derivative_orders"] == [6]
            assert data["included_correction_terms"][6]["weight"] == weight
            assert data["first_inactive_threshold"] > nu
            assert data["all_order_frequency"] > 0
            assert data["all_order_frequency"] != data["frozen_frequency"]
            assert sp.factor(data["all_order_frequency_slope"]-2*data["all_order_frequency"]*data["all_order_half_log_rate"]) == 0


def test_Minkowski_seed_keeps_three_positive_constrained_vector_polarizations():
    assert_zero(proca.minkowski_seed())
    assert_zero(proca.reconstruction_checks())


def test_state_change_and_complex_regulator_integrability_constants():
    assert all(value is True for value in transfer.proof_checks().values())
    assert all(value is True for value in transfer.scale_checks().values())
    assert len(transfer.scale_checks()) == 5


def test_physical_transfer_keeps_the_finite_matching_prescription():
    d = transfer.physical_bounds(10**12, 1000)
    assert d["state_change_energy_or_pressure_over_reference_density"] == sp.Rational(1, 5*10**21)
    assert d["state_change_density_derivative_over_reference_density_rate"] == sp.Rational(1, 10**18)
    assert d["full_quantum_solution_or_VGB_claim"] is False
