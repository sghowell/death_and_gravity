from fractions import Fraction

import sympy as sp
from p8.signs import even_rational_positive
from p8_s5 import nonlinear_d as m
from p8_s5 import scales


def test_all_time_weighted_lapse_bound():
    J = m.functions()["J"]
    even_rational_positive(m.d*J-sp.Rational(1, 2), m.u)
    assert sp.factor(2/m.d-J) == 4*m.u**4*(m.u**2+3)/m.d**5


def test_gamma_chart_finite_q_kinetic_bound():
    f = m.functions()
    even_rational_positive(-f["Lambda"]-sp.Rational(1, 2), m.u, upper=Fraction(1, 16))
    q = scales.derive()["q"]
    assert (8/(q-8)).subs(q, 64) == sp.Rational(1, 7)
    assert scales.derive()["residuals"]["relative_error"] == 0


def test_curvature_reference_does_not_vanish_at_bounce():
    inside, outside = scales.derive()["E_curvature_squared_over_E_ref_squared"]
    assert inside == 2
    x = next(iter(outside.free_symbols))
    assert sp.cancel(outside-2-4*(x-1)/(x+1)) == 0
    assert sp.cancel(6-outside-8/(x+1)) == 0
    assert scales.derive()["curvature_bounds_residuals"]["inner_equals_two"] == 0


def test_scale_restoration_and_background_variation():
    output = scales.derive()
    assert output["residuals"]["physical_H_definition"] == 0
    assert output["residuals"]["A3_variation_squared_bound"] == 36/m.d**2
    # Four-volume times either F, F2 R or Ai Li gives M^2 tau^2.
    M, tau = sp.symbols("M tau", positive=True)
    assert sp.cancel(tau**4*M**2/tau**2-output["action_prefactor"]) == 0


def test_dimensionful_family_solves_covariant_background_and_principal_dictionary():
    assert all(value == 0 for value in scales.scaled_background_checks().values())
