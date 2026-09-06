import pytest
import sympy as sp
from p8_composite_light import background, bounds, canonical, independent


@pytest.mark.parametrize("check", [canonical.physical_checks, canonical.generic_action_checks,
                                  canonical.operator_checks, canonical.routh_checks, background.checks])
def test_exact_physical_action_operator_charge_and_full_ODE_checks(check):
    assert check()
    assert all(sp.simplify(value) == 0 for value in check().values())


def test_independent_fraction_coefficients_and_bounds():
    result = independent.checks()
    assert len(result["coefficientwise_identities"]) == 4
    assert all(value == "0" for value in result["coefficientwise_identities"].values())


def test_actual_asymmetric_hierarchy_bounds_are_strict():
    d = bounds.build()
    assert all(value > 0 for value in d["strict_positive_margins"].values())
    assert d["strict_mass_led_ratio_lower"] == sp.Rational(50, 9)


def test_symmetric_zero_mass_and_positive_Routh_are_distinct():
    d = background.initial_jets(1)
    assert d["mass_squared"] == 0
    assert d["Routh_frequency_squared"].subs(background.A, sp.Rational(4, 25)) == sp.Rational(153, 100)
    assert d["heavy_diagonal_at_zero_k"].subs(background.A, sp.Rational(4, 25)) == -sp.Rational(43579, 52500)


def test_asymmetric_locked_cone_does_not_certify_a_heavy_gap():
    d = bounds.build()
    cone = d["conditional_locked_speed_squared"]
    assert 0 < cone["lower"] < cone["upper"] < 1
    assert d["asymmetric_mass_squared_over_m_squared"]["upper"] < sp.Rational(1, 20)
    assert d["small_A_omega_squared_lower"] > sp.Rational(1, 3)


@pytest.mark.parametrize("value", [1.0, 2.0, True, 3, sp.Rational(1, 2)])
def test_initial_state_contract_is_not_silently_generalized(value):
    with pytest.raises(ValueError):
        background.initial_jets(value)


@pytest.mark.parametrize("value", [0.1, sp.Float("0.1"), sp.oo, sp.nan])
def test_approximate_or_nonfinite_enclosure_inputs_rejected(value):
    with pytest.raises((TypeError, ValueError, sp.PolynomialError)):
        bounds.linear_root_enclosure(value)


def test_nonlinear_root_expression_is_outside_affine_bound_contract():
    symbolic = sp.Symbol("uncontracted_coefficient", real=True)
    with pytest.raises(TypeError):
        bounds.linear_root_enclosure(symbolic)


def test_omitting_moving_eigenvector_changes_heavy_forcing():
    t = sp.Symbol("T", real=True)
    assert canonical.B(t, t, 1, 0, 0) == 2
    assert canonical.B(t, t, 0, 0, 0) == 0


def test_omitting_heavy_homogeneous_data_changes_light_equation():
    t = sp.Symbol("T", real=True)
    hom = sp.sin(t)
    assert canonical.B_adjoint(hom, t, 1, 0, 0) == -2*sp.cos(t)


def test_nonconstant_mass_residual_retains_mass_derivatives():
    t = sp.Symbol("T", positive=True)
    light = t
    h0 = canonical.mass_led_pair(light, t, 1, 0, 0, t**2)
    assert h0 == -2/t**2
    # With V_HH=m_alg², the remaining exact heavy residual is h0''.
    assert sp.diff(h0, t, 2) == -12/t**4


def test_zero_algebraic_mass_is_not_used_in_inverse_mass_representative():
    t = sp.Symbol("T", real=True)
    with pytest.raises(ValueError, match="nonzero mass_squared"):
        canonical.mass_led_pair(t, t, 1, 0, 0, 0)
