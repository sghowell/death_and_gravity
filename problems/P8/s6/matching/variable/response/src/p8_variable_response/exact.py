"""Strict exact domains for the finite-parameter response theorem."""

from fractions import Fraction

import sympy as sp

DELTA_MAX = sp.Rational(1, 10**9)
SLICE = sp.Rational(1, 100)
MU = sp.sqrt(39)/2


def number(value, name="value"):
    if isinstance(value, (bool, float)) or not isinstance(value, (int, Fraction, sp.Basic)):
        raise TypeError(f"{name} must be an exact real number, not a bool or float")
    value = sp.sympify(value)
    if value.has(sp.Float, sp.oo, sp.zoo, sp.nan) or value.is_number is not True:
        raise ValueError(f"{name} must be finite, exact and numerical")
    if value.is_real is not True or value.is_finite is not True:
        raise ValueError(f"{name} must be finite and real")
    return value


def nonnegative(value, name="value"):
    value = number(value, name)
    if value.is_nonnegative is not True:
        raise ValueError(f"{name} must be nonnegative")
    return value


def parameters(delta, kbar=1):
    delta, kbar = number(delta, "delta"), number(kbar, "kbar")
    if delta.is_positive is not True or (DELTA_MAX-delta).is_nonnegative is not True:
        raise ValueError("Require 0 < delta <= 10^-9; delta=0 is not a regular parent")
    if (kbar-1).is_nonnegative is not True or (2-kbar).is_nonnegative is not True:
        raise ValueError("The fixed momentum band is 1 <= kbar <= 2")
    return delta, kbar


def side(value):
    if isinstance(value, bool) or value not in (-1, 1):
        raise ValueError("The side must be the exact integer -1 or +1")
    if not isinstance(value, (int, sp.Integer)):
        raise TypeError("The side must be the exact integer -1 or +1")
    return int(value)
