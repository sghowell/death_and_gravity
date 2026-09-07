"""Exact domains; complex analytic radii are proof devices, not new parents."""

import sympy as sp

RADIUS = sp.Rational(1, 20)
DELTA_RADIUS = RADIUS**2
SLICE = sp.Rational(1, 100)
DELTA_MAX = sp.Rational(1, 10**9)


def number(value, name="value"):
    if isinstance(value, (bool, float)):
        raise TypeError(f"{name} must be exact, finite and real")
    value = sp.sympify(value)
    if (value.has(sp.Float, sp.oo, sp.zoo, sp.nan) or value.is_number is not True
            or value.is_real is not True or value.is_finite is not True):
        raise ValueError(f"{name} must be an exact finite real number")
    return value


def parameters(delta, momentum_squared=1, *, transfer=False):
    delta = number(delta, "delta")
    momentum_squared = number(momentum_squared, "K")
    lower = 1 if transfer else 0
    if (delta.is_positive is not True or (DELTA_MAX-delta).is_nonnegative is not True
            or (momentum_squared-lower).is_nonnegative is not True
            or (4-momentum_squared).is_nonnegative is not True):
        raise ValueError(f"Require 0<delta<=1e-9 and {lower}<=K<=4")
    return delta, momentum_squared


def nonnegative(value, name="value"):
    value = number(value, name)
    if value.is_nonnegative is not True:
        raise ValueError(f"{name} must be nonnegative")
    return value


def integer(value, name="index"):
    value = nonnegative(value, name)
    if value.is_integer is not True:
        raise ValueError(f"{name} must be an integer")
    return int(value)
