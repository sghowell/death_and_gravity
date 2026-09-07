"""Generic finite beta, robust caps and exact source budget tests."""

import pytest
import sympy as sp
from p8a_maxwell_cosmology import calibration, geometry
from p8a_maxwell_focusing import envelope


def test_expanded_affine_constants_are_exact_and_subcritical():
    data = calibration.affine_cost()
    assert data["C0"] == sp.Rational(1440182116554809, 112500000)
    assert data["Cbeta"] == sp.Rational(3345309847373, 10500000)
    assert data["C_at_abs_beta_one"] == data["C0"]+data["Cbeta"]
    assert data["C0_upper_margin"] > 0 and data["Cbeta_upper_margin"] > 0


@pytest.mark.parametrize("beta", [0, 1, -1, sp.Rational(7, 3), -10**8])
def test_every_finite_beta_cost_matches_immutable_parent(beta):
    actual = envelope.coefficients(geometry.EXPANDED_CAPS, geometry.FUTURE_CAPS,
                                   geometry.RATIO, beta_m=beta)["total_cost"]
    assert calibration.cost_at_beta(beta) == actual


def test_joint_gate_is_not_a_fixed_beta_or_delta_alone_assumption():
    assert all(sp.simplify(value) == 0 for value in calibration.identities().values())
    data = calibration.theorem_gate(sp.Rational(1, 10**8), 5)
    assert data["strict_focusing_margin_lower"] == sp.Rational(47079, 350000)
    assert data["margin_above_one_eighth"] == sp.Rational(3329, 350000)
    # delta=1e-8 with |beta|=100 is not inside this uniform joint gate.
    with pytest.raises(ValueError):
        calibration.theorem_gate(sp.Rational(101, 10**8), 5)


def test_excess_source_budget_is_rejected():
    with pytest.raises(ValueError):
        calibration.theorem_gate(0, 6)


def test_tight_reference_regression_does_not_replace_the_robust_constants():
    data = calibration.calibration()
    tight = data["tight_reference_regression"]
    assert tight["C0"] == sp.Rational(1583767184219, 126000)
    assert tight["Cbeta"] == sp.Rational(364312, 875)
    assert tight["C0"] < data["expanded_affine_cost"]["C0"]
