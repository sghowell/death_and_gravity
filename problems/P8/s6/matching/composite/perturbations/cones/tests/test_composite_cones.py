import pytest
import sympy as sp
from p8_composite_cones import background, cone, independent


@pytest.mark.parametrize("name", list(cone.checks()))
def test_literal_cone_and_pressure_identity(name):
    assert cone.checks()[name] == 0


@pytest.mark.parametrize("name", list(background.checks()))
def test_full_lapse_background_identity(name):
    assert background.checks()[name] == 0


def test_independent_coefficientwise_replay():
    assert len(independent.checks()["coefficientwise_identities"]) == 6


def test_positive_pressure_control_if_vector_degeneracy_is_allowed():
    d = background.coincident_pressure_control()
    assert d["coincident_lapse_identity"] == d["initial_H"] == 0
    assert d["initial_H_prime"] == sp.Rational(1, 6)
    assert d["initial_y_prime"] == -sp.sqrt(sp.Rational(7, 3))
    assert d["initial_null"] == 1
    assert d["initial_Ng"] == d["initial_Nf"] == sp.Rational(1, 2)


def test_omission_and_outside_domain_controls():
    controls = cone.controls()
    assert controls["g_can_exceed_matter_cone"] == sp.Rational(5, 4)
    assert controls["f_can_exceed_matter_cone"] == sp.Rational(7, 9)
    assert controls["single_metric_beta_zero_does_not_force_c_equal_y"] == 0
    assert controls["omitted_pressure_falsely_keeps_relative_stiffness"] != 0
    assert controls["squared_speed_not_velocity_average"] != 0


def test_zero_and_negative_Xi_are_not_regular_positive_all_k_vectors():
    from p8_composite_modes import model as m
    from p8_composite_modes import vector

    spring = vector.spring()
    assert spring["harmonic"].subs(spring["C"], 0) == 0
    d = vector.derive()
    fixture = {m.G: 1, m.F: 1, m.a: 1, m.Ng: 1, m.y: 1, m.c: 1,
               d["Xi"]: -1, d["k"]: 2}
    assert d["K"].subs(fixture) == -sp.Rational(1, 32)
