"""Small exact-rational public-input boundary; symbolic derivations are separate."""

from fractions import Fraction

import sympy as sp


def rational(value):
    """Accept finite rationals only, never silently rationalize a binary float."""
    if isinstance(value, (bool, float, sp.Float)):
        raise TypeError("Expected an exact rational, not bool or float")
    if isinstance(value, Fraction):
        return sp.Rational(value.numerator, value.denominator)
    if isinstance(value, (int, sp.Rational)):
        return sp.Rational(value)
    if isinstance(value, str):
        try:
            result = sp.Rational(value)
        except (TypeError, ValueError, ZeroDivisionError) as exc:
            raise ValueError("Expected a finite rational string") from exc
        if result.is_finite is True:
            return result
    raise TypeError("Expected a finite exact rational")


def coefficients(beta):
    if not isinstance(beta, (tuple, list)) or len(beta) != 5:
        raise ValueError("A link has exactly five constant coefficients beta0..4")
    return tuple(rational(value) for value in beta)


def positive(value, name):
    result = rational(value)
    if result <= 0:
        raise ValueError(f"{name} must be strictly positive")
    return result


def nonnegative(value, name):
    result = rational(value)
    if result < 0:
        raise ValueError(f"{name} must be nonnegative")
    return result
