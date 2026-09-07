"""Exact full density/pressure equation and open low-branch controls."""

import pytest
import sympy as sp
from p8a_maxwell_thermal import dynamics


def test_all_exact_see_and_clock_identities():
    assert all(sp.simplify(value) == 0 for value in dynamics.identities().values())


def test_named_initial_constraint_is_actual_not_a_frozen_radiation_substitution():
    data = dynamics.branch_point(2, dynamics.DELTA_MAX)
    assert data["a_fourth_with_a0_one"] == 1
    assert data["normal_Hubble_jets_scaled_by_tau"][1] != -8
    assert data["kappa_EED_times_tau_squared"] > 0


def test_endpoint_is_finite_positive_but_outside_the_regular_domain():
    data = dynamics.endpoint_data(dynamics.DELTA_MAX)
    assert 0 < data["a_endpoint_fourth"] < 1
    assert data["endpoint_x_strict_lower"] == sp.Rational(1, 16)
    assert data["endpoint_x_strict_upper"] == sp.Rational(1, 4)
    assert data["endpoint_is_excluded_from_smooth_metric_domain"]
    assert not data["fundamental_EFT_validity_at_endpoint_proved"]


def test_positive_scale_factor_on_an_algebraic_high_branch_is_not_admitted():
    # This y has lambda*y² between1/2 and1: a⁴ is positive, but the branch is wrong.
    lam = dynamics.dimensionless_coupling(dynamics.DELTA_MAX)
    assert sp.Rational(1, 2) < lam*20000**2 < 1
    with pytest.raises(ValueError, match="low-curvature"):
        dynamics.branch_point(20000, dynamics.DELTA_MAX)


def test_h4_omission_changes_the_initial_source_charge():
    assert dynamics.calibration()["Q_in_one_over_kappa_tau_squared_units"] != 12


@pytest.mark.parametrize("bad", [True, 0, -1, 0.01, sp.oo, sp.nan, sp.Symbol("delta"), sp.Rational(1, 10**7)])
def test_actual_physical_branch_requires_positive_exact_small_delta(bad):
    with pytest.raises((TypeError, ValueError)):
        dynamics.dimensionless_coupling(bad)
