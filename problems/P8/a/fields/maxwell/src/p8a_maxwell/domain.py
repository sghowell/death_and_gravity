"""Strict exact inputs for numerical enclosures; geometry remains a hypothesis."""

from fractions import Fraction

import sympy as sp


def rational(value):
    if isinstance(value, (bool, float, sp.Float)):
        raise TypeError("use exact rational inputs, not booleans or binary floats")
    if isinstance(value, Fraction):
        value = sp.Rational(value.numerator, value.denominator)
    result = sp.sympify(value)
    if not isinstance(result, sp.Rational):
        raise TypeError("a finite exact rational input is required")
    return result


def nonnegative(value, *, positive=False):
    if not isinstance(positive, bool):
        raise TypeError("positive must be boolean")
    result = rational(value)
    if result < 0 or (positive and result == 0):
        raise ValueError("the input is outside the stated positive domain")
    return result


def support_interval(domain_left, domain_right, left, right):
    """Check enclosing compact conformal intervals, not a field or sampler."""
    dl, dr, left, right = map(rational, (domain_left, domain_right, left, right))
    if not dl < left < right < dr:
        raise ValueError("sampler support must be compact inside the conformal strip")
    return True
