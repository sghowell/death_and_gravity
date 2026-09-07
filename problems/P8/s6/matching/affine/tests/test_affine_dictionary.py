"""Original witness, general inverse and nonunit source conventions."""
from fractions import Fraction

import pytest
import sympy as sp
from p8_affine import dictionary as d


@pytest.mark.parametrize("name", ["f", "alpha1", "alpha2", "alpha3", "alpha4", "alpha5",
                                  "Delta_one", "p_square", "inverse_ODE", "clock_normalization"])
def test_original_cd_principal_dictionary(name):
    assert d.principal_identities()[name] == 0


@pytest.mark.parametrize("name", ["target_f", "target_A1", "target_A3", "two_displayed_denominators"])
def test_general_inverse(name):
    assert d.generic_inverse_identities()[name] == 0


@pytest.mark.parametrize("time,expected", [(0, (sp.Rational(-1), sp.Rational(11, 4), 1)),
                                          (1, (sp.Rational(-1, 8), sp.Rational(39, 256), sp.Rational(1, 64)))])
def test_distinct_clock_points(time, expected):
    values = d.target()
    assert tuple(sp.simplify(values[key].subs({d.u: time, d.x: -1}))
                 for key in ("alpha3", "alpha4", "alpha5")) == expected


@pytest.mark.parametrize("physical_X", [sp.Rational(9, 10), 1, sp.Rational(11, 10), Fraction(2401, 2500)])
def test_closed_domain_calibrations(physical_X):
    result = d.tube_point(physical_X=physical_X)
    assert sp.Rational(9, 40) <= result["dimensionless_p_squared"] <= sp.Rational(11, 40)
    assert result["additional_quotient_factor"] >= sp.Rational(4, 5)


def test_nonunit_coefficient_and_action_units():
    result = d.tube_point(physical_X=Fraction(2401, 2500), mass_squared=3, time_scale=2)
    assert result["p_physical"] == sp.Rational(147, 100)
    assert result["Delta_physical"] == 3
    assert result["cubic_and_effective_braiding_scale"] == sp.Rational(3, 2)
    assert result["parent_scalar_coefficient_scale"] == sp.Rational(3, 4)
    assert result["four_volume_normalized_action_scale"] == 12


@pytest.mark.parametrize("value", [True, False, sp.true, sp.false, 1.0, sp.Float(1), "1", sp.I, sp.oo, sp.nan])
def test_inexact_or_nonphysical_inputs_rejected(value):
    with pytest.raises((TypeError, ValueError)):
        d.tube_point(physical_X=value)


@pytest.mark.parametrize("value", [0, -1, sp.Rational(1, 2), sp.Rational(6, 5)])
def test_outside_actual_tube_rejected(value):
    with pytest.raises(ValueError):
        d.tube_point(physical_X=value)


@pytest.mark.parametrize("kwargs", [{"mass_squared": 0}, {"mass_squared": -1}, {"time_scale": 0}, {"time_scale": -1}])
def test_nonpositive_units_rejected(kwargs):
    with pytest.raises(ValueError):
        d.tube_point(**kwargs)


def test_scalar_input_is_original_nontrivial_witness():
    assert d.target()["braiding"] == 0
    assert d.target()["scalar_F"].has(d.x, d.u)
    assert sp.diff(d.target()["scalar_F"], d.x, 2) != 0
