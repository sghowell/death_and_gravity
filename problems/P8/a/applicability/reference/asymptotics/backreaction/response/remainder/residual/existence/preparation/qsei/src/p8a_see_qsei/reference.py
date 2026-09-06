"""Positive reference credit from the NEW exact SEE, not old-metric transfer."""

from functools import cache

import sympy as sp
from p8a_preparation import bounds as prior

from .scattering import nonnegative


def energy_dimensionless(scale, hubble, potential, delta):
    """Actual reference EED in hbar/(pi^2 A^4 eta_star^8) units.

Requires the unforced SEE with the original ordinary radiation constant.
It is not a formula for arbitrary target states on the fixed geometry.
    """
    a, h, u, delta = map(sp.sympify, (scale, hubble, potential, delta))
    return (a**2*(u+h**2)-1)/(960*delta*a**4)


@cache
def calibration():
    data = prior.calibration()
    delta, length = data["delta"], data["length"]
    distance = data["actual_gate"]["fixed_point_distance"]
    cap = sp.Rational(72, 10**9)
    b = data["geometry"]["u_cap"]
    loss = (9*length+sp.Rational(9, 2)*length**2+(3*b+sp.Rational(3, 4))*length**3)*cap
    numerator = delta/54
    eed = numerator/(960*delta*81)
    margins = {"distance": cap-distance, "numerator_loss": delta/54-loss}
    if any(value <= 0 for value in margins.values()):
        raise ValueError("the actual exact-solution reference credit failed")
    return {"actual_fixed_point_distance": distance, "distance_cap": cap,
            "old_numerator_lower": delta/27, "numerator_loss_upper": loss,
            "new_numerator_lower": numerator, "positive_EED_dimensionless": eed,
            "strict_margins": margins}


def positive_credit(*, normalization=1, eta_star=1, hbar=1):
    normalization = nonnegative(normalization, positive=True)
    eta_star = nonnegative(eta_star, positive=True)
    hbar = nonnegative(hbar)
    return hbar*calibration()["positive_EED_dimensionless"]/(sp.pi**2*normalization**4*eta_star**8)


def identities():
    x = sp.Symbol("x", real=True)
    a = sp.Function("a", positive=True)(x)
    delta, y = sp.symbols("delta y", positive=True)
    h, u = sp.diff(a, x)/a, -sp.diff(a, x, 2)/a
    proper_h = h/a
    geometric = -3*(sp.diff(proper_h, x)/a+proper_h**2)/(2880*delta)
    radiation = 3/(2880*delta*a**4)
    f = 1-delta/y**4
    aa, hh, uu = y*f**sp.Rational(1, 4), 1/(y*f**sp.Rational(3, 4)), 2*delta/(y**6*f**sp.Rational(3, 2))
    return {
        "actual_SEE_EED_minus_fixed_radiation": sp.simplify(
            geometric-radiation-energy_dimensionless(a, h, u, delta)),
        "old_geometry_numerator": sp.simplify(aa**2*(uu+hh**2)-1-3*(delta/y**4)/f),
        "old_geometry_EED": sp.simplify(energy_dimensionless(aa, hh, uu, delta)-1/(320*y**8*f**2)),
    }
