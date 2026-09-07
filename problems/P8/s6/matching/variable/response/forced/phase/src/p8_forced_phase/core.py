"""Exact continuous amplitude and error arithmetic, not mode sampling."""

from functools import cache

import sympy as sp
from p8_variable_forced import operator

R_MAX = sp.Rational(1, 1000)
DEFICIT_MAX = sp.Rational(1, 100)
MU = sp.sqrt(39)/2


def rational(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, sp.Rational)):
        raise TypeError(f"{name} must be an exact integer or SymPy rational")
    return sp.Rational(value)


def calibration(r=R_MAX, deficit=DEFICIT_MAX):
    r, deficit = rational(r, "r"), rational(deficit, "deficit/r")
    if not 0 < r <= R_MAX or not 0 <= deficit <= sp.Rational(1, 100):
        raise ValueError("Require 0<r<=1/1000 and 0<=deficit/r<=1/100")
    t = r*r
    moment = sp.Rational(5, 28)-deficit-4*t
    amplitude = sp.Rational(32, 65)*(1-9*t)*moment
    extra = sp.Rational(202, 15)*(10+60*t)
    return {"r": r, "deficit_over_r": deficit,
            "moment_over_r_three_halves": moment,
            "amplitude_over_r_squared_lower": amplitude,
            "additional_feedback_error_over_r_fourth": extra,
            "simplified_response_error_over_r_fourth": 1000+extra,
            "rounded_amplitude_over_r_squared_lower": sp.Rational(2, 25),
            "rounded_response_error_over_r_fourth": sp.Integer(1200),
            "phase_separation_over_r_squared_lower": sp.Rational(4, 25)-2400*t,
            "claimed_separation_over_r_squared": sp.Rational(3, 20),
            "smooth_transition_width_over_r": sp.Rational(1, 800),
            "smooth_pulse_deficit_over_r_upper": sp.Rational(1, 200)}


def separation_bound(r=R_MAX, amplitude=1, momentum_squared=0):
    """Conditional lower bound; numerical inputs do not certify a profile."""
    c = calibration(r)
    amp, kk = rational(amplitude, "amplitude"), rational(momentum_squared, "K")
    if amp <= 0 or not 0 <= kk <= 4:
        raise ValueError("Require a positive scaling amplitude and 0<=K<=4")
    return amp*c["claimed_separation_over_r_squared"]*c["r"]**2


@cache
def checks():
    t, q = sp.symbols("t q", nonnegative=True)
    d = operator.derive()
    u, delta = d["u"], d["delta"]
    j0 = -2/sp.sqrt(5)
    map0 = -4/sp.sqrt(5)
    jr = d["jH"].subs({u: sp.sqrt(t), delta: 0})/j0
    gr = d["bg"].subs({u: sp.sqrt(t), delta: 0})/map0
    exponent = sp.Rational(3, 2)+sp.I*MU
    x = sp.Symbol("x", positive=True)
    residuals = {
        "physical_source_squared_ratio": sp.factor(jr**2-5*(1+t)**6/((1+t)**12+4)),
        "source_lower_polynomial": sp.expand(5*q-q*q-4-(q-1)*(4-q)),
        "physical_map_squared_ratio": sp.factor(gr**2-5/((1+t)**6*((1+t)**12+4))),
        "source_upper_cubic": sp.expand((1+t)**3-1-3*t-3*t*t-t**3),
        "radial_moment_primitive": sp.simplify((sp.diff(x**exponent/exponent, x)-sp.sqrt(x)*sp.exp(sp.I*MU*sp.log(x))).rewrite(sp.exp)),
        "moment_denominator_modulus": sp.expand((sp.Rational(3, 2)+sp.I*MU)*(sp.Rational(3, 2)-sp.I*MU)-12),
        "physical_constant_product": sp.simplify(j0*map0-sp.Rational(8, 5)),
        "double_clock_scale": sp.Integer(4)**2*2-32,
    }
    c = calibration()
    margins = {
        "source_near_constant_domain": 4-(1+R_MAX**2)**6,
        "cubic_envelope": 4-(3+3*R_MAX**2+R_MAX**4),
        "inverse_ninth_positive": 1-9*R_MAX**2,
        "moment_survives_pulse_smoothing": c["moment_over_r_three_halves"],
        "nonzero_physical_phase_amplitude": c["amplitude_over_r_squared_lower"]-sp.Rational(2, 25),
        "full_feedback_and_pole_error": 1200-c["simplified_response_error_over_r_fourth"],
        "strict_limsup_liminf_separation": c["phase_separation_over_r_squared_lower"]-sp.Rational(3, 20),
        "smooth_pulse_deficit": sp.Rational(1, 100)-c["smooth_pulse_deficit_over_r_upper"],
        "radial_two_power_bound_squared": sp.Rational(9, 64)-sp.Rational(1, 8),
        "moment_root_twelve_bound_squared": sp.Rational(49, 4)-12,
        "frequency_upper_bound_squared": sp.Rational(169, 16)-MU**2,
    }
    if any(value != 0 for value in residuals.values()):
        raise ValueError(f"Exact residuals failed: { {key: value for key, value in residuals.items() if value != 0} }")
    if any(value.is_positive is not True for value in margins.values()):
        raise ValueError("A continuous phase-separation margin failed")
    return {"residuals": residuals, "margins": margins}
