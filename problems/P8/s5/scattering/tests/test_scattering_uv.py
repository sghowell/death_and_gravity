import sympy as sp
from p8_s5 import nonlinear_d as m
from p8_scattering import uv_audit


def test_specific_polynomial_continuation_has_no_constant_clock_flat_vacuum():
    out = uv_audit.derive()
    assert out['minus_F_strictly_positive']['numerator']['roots'] == 0
    F = m.functions()['F']
    assert F.subs({m.u: 0, m.X: 0}) == -19
    assert sp.limit(m.u**2*F.subs(m.X, 0), m.u, sp.oo) == -10


def test_complete_geometry_is_not_an_asymptotic_high_frequency_band_for_fixed_modes():
    out = uv_audit.derive()
    assert out['conformal_time_length_over_tau'] == 'pi'
    assert out['past_ratio_limit'] == out['future_ratio_limit'] == '0'
