"""Independent action, time-scale, momentum and final-germ readout checks."""
from fractions import Fraction

import sympy as sp
from p8_own_retarded import action


def test_literal_action_and_original_g_readout_identities():
    values = action.identities()
    assert len(values) == 10
    assert set(values.values()) == {0}


def test_variable_kinetic_drift_is_retained():
    data = action.normalization()
    k, h, s, q = (data[name] for name in ("k", "hidden", "s", "q"))
    frozen_kinetic = k*sp.diff(h, action.X, 2)+s*(h-q)
    assert sp.expand(data["own_Euler"]-frozen_kinetic) == sp.diff(k, action.X)*sp.diff(h, action.X)


def test_physical_time_and_action_measure_are_not_dropped():
    d = action.normalization()
    assert d["scale"] == action.TAU*sp.sqrt(action.DELTA)
    assert d["physical_Euler_conversion"] == action.M2/(action.TAU**2*action.DELTA)
    assert d["dictionary"][d["s"]] == 2*action.TAU**2*action.DELTA*d["beta1"]*d["b"]/action.M2


def test_independent_sharp_green_constants_and_response_floors():
    d = action.independent_volterra_constants()
    assert d["positive_flux"] == Fraction(283, 608)
    assert d["R_x"] == Fraction(1415, 12768)
    assert d["R_over_duration"] == Fraction(7495, 38304)
    assert d["R_margin"] == Fraction(313, 38304)
    assert d["R_x_margin"] == Fraction(691, 63840)
    assert d["Q0_lower_per_amplitude"] == Fraction(63, 1024)
    assert d["Qx0_lower_per_amplitude"] == Fraction(21, 80)


def test_original_metric_final_euler_difference_and_nonclaims():
    d = action.readout()
    assert sp.factor(d["normalized_Euler_magnitude_lower"]-3*d["eta"]) > 0
    assert d["normalized_original_g_Euler_at_zero_germ"] == -128*d["hidden0"]/action.C
    assert not d["is_a_conserved_matter_source_or_on_shell_physical_g_claim"]
    # A nonzero readout means keeping prescribed g may require driving;
    # it is not a solution of the undriven physical-g equation.


def test_constant_coefficient_calibration_does_not_set_variable_drift_to_zero():
    # A non-load-bearing exact calibration: k=4,s=64 gives frequency4,
    # not the coupled relative frequency sqrt(80).
    duration = sp.Symbol("duration", real=True)
    kernel = sp.sin(4*duration)/16
    assert sp.diff(kernel, duration).subs(duration, 0) == sp.Rational(1, 4)
    assert sp.diff(kernel, duration, 2)+16*kernel == 0
    assert action.independent_volterra_constants()["R_x"] < sp.Rational(1, 8)
