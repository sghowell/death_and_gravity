"""Uniform spectral quadratic remainder and finite slope estimates."""

from functools import cache

import sympy as sp
from p8_vacuum_global_one_loop_insertions import halfplane

N = sp.Integer(6)


def bound(mass, yukawa_upper, cubic_squared, heavy_mass_squared, Q_lower):
    m, Y, g, M, Q = map(
        halfplane.exact,
        (mass, yukawa_upper, cubic_squared, heavy_mass_squared, Q_lower),
    )
    if m < 2 or Y < 0 or g < 0 or M < 1 or Q <= 0:
        raise ValueError("Require m>=2, Y>=0, g>=0, M>=1 and Q>0")
    T = 4 * m * m
    C = 2 * N * Y / Q
    slope = 4 * C * g / (Q * T)
    remainder = 4 * C * g / (3 * Q * T * T)
    return {
        "fermion_mass": m,
        "Y_upper": Y,
        "g": g,
        "M": M,
        "Q_lower": Q,
        "spectral_threshold": T,
        "finite_positive_outer_slope_upper": slope,
        "uniform_outer_OS_remainder_coefficient_upper": remainder,
        "strict_family_slope_and_spacelike_remainder_positive": bool(Y > 0 and g > 0),
        "scope": "Only this first covariance-insertion quadratic family. The remainder bound is uniform on Re(s)<=4; it is not a complete two-loop pole error or exact heavy spectrum.",
    }


@cache
def data():
    x, u, T, C, g, Q, m, Y = sp.symbols("x u T C g Q m Y", positive=True)
    remainder = 4 * C * g / (3 * Q * T * T)
    slope = 4 * C * g / (Q * T)
    checks = {
        "first_derivative_parameter_moment": sp.integrate(1 - x, (x, 0, 1))
        - sp.Rational(1, 2),
        "second_derivative_parameter_moment": sp.integrate((1 - x) ** 2, (x, 0, 1))
        - sp.Rational(1, 3),
        "spectral_first_moment": sp.integrate(u ** (-2), (u, T, sp.oo)) - 1 / T,
        "spectral_second_moment": sp.integrate(u ** (-3), (u, T, sp.oo))
        - 1 / (2 * T * T),
        "outer_Taylor_half_not_doubled": (g / Q)
        * sp.Rational(1, 2)
        * (16 * C / 3)
        * (1 / (2 * T * T))
        - remainder,
        "finite_outer_slope_bound_factor": (g / Q) * (4 * C) * (1 / T) - slope,
        "mass_coupling_remainder_dictionary": sp.factor(
            remainder.subs({C: 2 * N * Y / Q, T: 4 * m * m}, simultaneous=True)
            - N * Y * g / (6 * Q * Q * m**4)
        ),
        "slope_remainder_ratio": sp.factor(slope / remainder - 3 * T),
    }
    return {
        "uniform_remainder_upper": remainder,
        "finite_slope_upper": slope,
        "large_spectral_mass_bounds": "w(u)<=4C/u, u-4>=u/2. The outer derivative majorants are 1/[2(u-4)] and 1/[3(u-4)^2].",
        "analytic_remainder": "|Pi_spec,R(s)|<=B_spec |s-1|^2 on Re(s)<=4, B_spec=4Cg/(3QT^2).",
        "spacelike_bound": "0<=Pi_spec,R(-t)<=(t+1) r_spec, 0<=r_spec<=4Cg/(QT). Strict positivity for nonzero Y and g.",
        "scope": "The whole proper inner subtraction and complete outer OS subtraction precede each bound. No unregulated tadpole or outer mass constant is assigned a finite value.",
        "checks": checks,
    }
