"""Exact full-operator, chart, source and phase omission controls."""

import sympy as sp
from p8_variable_response import connection, operator, source


def test_complete_chart_residuals():
    assert all(sp.simplify(value) == 0 for value in operator.checks().values())


def test_exact_hypergeometric_and_power_connection_residuals():
    assert all(sp.simplify(value) == 0 for value in connection.checks().values())


def test_actual_physical_source_residuals():
    assert all(sp.simplify(value) == 0 for value in source.checks().values())


def test_moving_weight_operator_is_not_a_locked_scalar():
    d = operator.derive()
    center = {d["u"]: 0, d["c"]: 2}
    assert sp.simplify(d["C"].subs(center)) == sp.Rational(48, 5)
    assert sp.simplify(d["E"].subs(center)) == 0
    # omega=0 is not permission to delete its derivative in the adjoint.
    assert d["omega"].subs(center) == 0
    assert sp.diff(d["omega"], d["u"]).subs(center) == -sp.Rational(24, 5)


def test_indicial_powers_are_complex_not_constant_heavy_data():
    mu = connection.MU
    for root in (sp.Rational(1, 2)+sp.I*mu, sp.Rational(1, 2)-sp.I*mu):
        assert sp.simplify(root*(root-1)+10) == 0
        assert sp.simplify(root*(root-1)) == -10  # omission residual


def test_source_off_g_jet_keeps_the_relative_spring():
    d, obs = operator.derive(), source.g_observability()
    for sign in (-1, 1):
        value = sp.factor(obs["relative_coefficient"].subs({d["u"]: sp.Rational(sign, 100), d["c"]: 2}))
        assert value.is_positive is True
    assert obs["matrix"][2, 2] != 0
    assert obs["matrix"][3, 3] != 0


def test_raw_proper_clock_is_not_silently_u_clock():
    delta = sp.Rational(1, 10**12)
    unit = operator.physical_at(sp.Rational(1, 100), delta)
    scaled = operator.physical_at(sp.Rational(1, 100), delta, tau=7)
    assert scaled == sp.diag(1, sp.Rational(1, 7), 1, sp.Rational(1, 7))*unit


def test_reflection_is_not_the_only_phase_dependent_entry():
    controls = connection.controls()
    assert controls["nonzero_reflection_squared"].is_positive is True
    assert sp.simplify(controls["transmission_not_small"]) == 0
    assert controls["B_zero_phase_derivative"].is_nonzero is True
    assert controls["canonical_Cauchy_symplectic_claim_from_power_determinant"] is False


def test_matching_sign_and_fractional_weight_cannot_be_deleted():
    eta = sp.Rational(1, 8000)
    plus, minus = operator.boundary_map(eta, 1), operator.boundary_map(eta, -1)
    assert plus[2, 2] != 1
    assert minus[3, 3] == -plus[3, 3]
    assert plus[3, 2] != 0


def test_source_couples_to_both_modes_at_the_center_limit():
    d = source.coefficients()
    at = {d["u"]: 0, d["c"]: 2}
    assert sp.simplify(d["j_heavy"].subs(at)/d["j_light"].subs(at)) == -2
    assert sp.simplify(d["physical_J_over_sigma"].subs(at)) == sp.Rational(1, 2)
