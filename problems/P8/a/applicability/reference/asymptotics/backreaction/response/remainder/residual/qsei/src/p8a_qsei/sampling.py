"""Exact spectral moments and proper-time all-H2-sampler bounds."""

from functools import cache

import sympy as sp
from p8a_remainder.mode_bounds import exact_nonnegative

from . import mode_bounds

DELTA_MAX = sp.Rational(1, 10**14)


def spectral_identities():
    alpha, k, u = sp.symbols("alpha k u", positive=True)
    a, b = sp.symbols("a b", real=True)
    return {"positive_frequency_UV_integral": sp.integrate((alpha+k)**-4, (alpha, 0, sp.oo))-1/(3*k**3),
            "full_UV_moment": sp.integrate(1/(3*k**2), (k, 1, sp.oo))-sp.Rational(1, 3),
            "IR_mode_measure": sp.integrate(k, (k, 0, 1))-sp.Rational(1, 2),
            "flat_k_squared_moment": sp.integrate(k**3, (k, 0, u))-u**4/4,
            "flat_constant_moment": sp.integrate(k, (k, 0, u))-u**2/2,
            "flat_cross_moment": sp.integrate(2*k**2, (k, 0, u))-2*u**3/3,
            "flat_positive_square_completion": sp.expand(a*a+2*b*b+sp.Rational(8, 3)*a*b
                                                       -(a+sp.Rational(4, 3)*b)**2-sp.Rational(2, 9)*b*b)}


def clock_identities():
    s = sp.Symbol("s", positive=True)
    a, psi = sp.Function("a", positive=True)(s), sp.Function("psi")(s)
    dx = lambda value: a*sp.diff(value, s)
    f = a**-sp.Rational(3, 2)*psi
    h, hubble = sp.diff(a, s), sp.diff(a, s)/a
    hubble_dot = sp.diff(hubble, s)
    y, delta = sp.symbols("y delta", positive=True)
    plateau_a = (4*s*s-delta)**sp.Rational(1, 4)
    actual_h = sp.diff(plateau_a, s)/plateau_a
    z = delta/(4*s*s)
    return {"sampler_zero_derivative_clock": sp.simplify(f/sp.sqrt(a)-psi/a**2),
            "sampler_first_derivative_clock": sp.simplify(dx(f)/sp.sqrt(a)
                -(sp.diff(psi, s)-sp.Rational(3, 2)*hubble*psi)/a),
            "sampler_second_derivative_clock": sp.simplify(dx(dx(f))/sp.sqrt(a)
                -sp.diff(psi, s, 2)+2*hubble*sp.diff(psi, s)
                -(sp.Rational(3, 4)*hubble**2-sp.Rational(3, 2)*hubble_dot)*psi),
            "weighted_sampler_derivative_clock": sp.simplify(dx(h*f)/sp.sqrt(a)
                -hubble*sp.diff(psi, s)-(hubble_dot-hubble**2/2)*psi),
            "exact_plateau_Hubble": sp.factor(actual_h-1/(2*s*(1-z))),
            "exact_plateau_Hubble_derivative": sp.factor(sp.diff(actual_h, s)+(1+z)/(2*s*s*(1-z)**2)),
            "proper_time_plateau_endpoints": sp.expand((y*y/2).subs(y, 3)-(y*y/2).subs(y, 2)-sp.Rational(5, 2))}


@cache
def calibration():
    data = mode_bounds.calibration()
    if DELTA_MAX > data["delta_bar"]:
        raise ValueError("QSEI amplitude cap exceeds the prepared-state domain")
    hubble, derivative = sp.Rational(8, 31), sp.Rational(132, 961)
    first_poincare, second_poincare = sp.Rational(5, 6), sp.Rational(25, 36)
    f0 = second_poincare/2
    f1 = first_poincare+sp.Rational(3, 2)*hubble*second_poincare
    f2 = 1+2*hubble*first_poincare+(sp.Rational(3, 4)*hubble**2
                                   +sp.Rational(3, 2)*derivative)*second_poincare
    weighted_first = hubble*first_poincare+(derivative+hubble**2/2)*second_poincare
    # 2 sqrt(2) >= 8/3 and sqrt(2)<3/2 prove this norm bound for the
    # exact positive flat-kernel quadratic form, even when its cross term is negative.
    auxiliary_root = f2+sp.Rational(3, 2)*weighted_first
    e0, e1, e2 = data["ultraviolet_error_derivatives_per_delta"]
    combo = e0*f2+2*e1*f1+e2*f0
    # W<=3 and the UV moment=1/3; 2/sqrt(pi)<2. The IR full-Parseval
    # bound contributes exactly 2*E_IR*||F|| to the same normalized root.
    error_root = 2*combo+2*data["infrared_error_per_delta"]*f0
    worst_coefficient = (auxiliary_root+DELTA_MAX*error_root)**2
    if worst_coefficient >= 5:
        raise ValueError("The explicit all-sampler coefficient C=5 is not proved")
    return {"dimensionless_proper_sampling_span": sp.Rational(5, 2),
            "Hubble_cap": hubble, "absolute_Hubble_derivative_cap": derivative,
            "inverse_scale_squared_cap": sp.Rational(1, 2),
            "sampler_norms_per_proper_second_derivative": [f0, f1, f2],
            "weighted_first_sampler_norm": weighted_first,
            "auxiliary_flat_kernel_root_coefficient": auxiliary_root,
            "mode_error_root_coefficient_per_delta": error_root,
            "maximum_amplitude": DELTA_MAX,
            "maximum_QSEI_coefficient": worst_coefficient,
            "rounded_QSEI_coefficient": sp.Integer(5),
            "strict_rounding_margin": 5-worst_coefficient}


def coefficient(delta):
    delta = exact_nonnegative(delta)
    if (DELTA_MAX-delta).is_nonnegative is not True:
        raise ValueError("The absolute all-sampler theorem requires delta<=10^-14")
    data = calibration()
    return (data["auxiliary_flat_kernel_root_coefficient"]
            +delta*data["mode_error_root_coefficient_per_delta"])**2


def absolute_bound(delta, proper_second_derivative_norm_squared, *, hbar=1):
    """Magnitude of the negative lower bound, in the named A.8 scheme.

For real H2_0 samplers supported in 2<y<3 and any Hadamard target state:
integral E_target*h^2 dt >= -absolute_bound(...). No SEE is assumed here.
"""
    norm = exact_nonnegative(proper_second_derivative_norm_squared)
    hbar = exact_nonnegative(hbar)
    return hbar*coefficient(delta)*norm/(16*sp.pi**2)
