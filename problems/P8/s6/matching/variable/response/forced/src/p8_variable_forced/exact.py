"""The r-variable, zero-spatial-momentum-inclusive physical domain."""

import sympy as sp

R_MAX = sp.Rational(1, 100)
DELTA_MAX = sp.Rational(1, 10**9)
MU = sp.sqrt(39)/2


def number(value, name="value"):
    if isinstance(value, (bool, float)):
        raise TypeError(f"{name} must be exact, finite and real")
    value = sp.sympify(value)
    if (value.has(sp.Float, sp.oo, sp.zoo, sp.nan) or value.is_number is not True
            or value.is_real is not True or value.is_finite is not True):
        raise ValueError(f"{name} must be an exact finite real number")
    return value


def radius(value):
    value = number(value, "r")
    if value.is_positive is not True or (R_MAX-value).is_nonnegative is not True:
        raise ValueError("Require 0<r<=1/100")
    return value


def parameters(delta, r=R_MAX, momentum_squared=1):
    delta, r, kk = number(delta, "delta"), radius(r), number(momentum_squared, "K")
    if (delta.is_positive is not True or (DELTA_MAX-delta).is_nonnegative is not True
            or (r**2/100-delta).is_nonnegative is not True
            or kk.is_nonnegative is not True or (4-kk).is_nonnegative is not True):
        raise ValueError("Require 0<delta<=min(1e-9,r^2/100) and 0<=K<=4")
    return delta, r, kk


def nonnegative(value, name="value"):
    value = number(value, name)
    if value.is_nonnegative is not True:
        raise ValueError(f"{name} must be nonnegative")
    return value
