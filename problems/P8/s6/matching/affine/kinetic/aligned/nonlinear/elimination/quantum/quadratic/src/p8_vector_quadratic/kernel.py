"""Accepted curved quadratic mass-insertion pole, separated by derivative order."""
from functools import cache, lru_cache

import sympy as sp

from . import bimetric, invariants, tensors


@lru_cache(maxsize=None, typed=True)
def pole(order):
    if type(order) is not int or order not in (0, 2, 4):
        raise ValueError("Require native covariant derivative order 0,2 or 4")
    if order == 4:
        return bimetric.fourier_bilinear()
    return invariants.evaluate(order)/(1 if order == 0 else 48*tensors.mass**2)


@cache
def independent_second():
    # Garcia-Recio--Salcedo (114), derivative index FIRST in Y_abc.
    # Converted from their epsilon=(d-4)/2 to ours (4-d)/2.
    Y, Ric, Riem, R = invariants.Y, invariants.Ric, invariants.Riem, invariants.R
    terms = ((-sp.Rational(1, 12), (Y("mn", "m"), Y("aa", "n"))),
             (-sp.Rational(1, 48), (Y("nn", "m"), Y("aa", "m"))),
             (sp.Rational(1, 24), (Y("na", "m"), Y("na", "m"))),
             (-sp.Rational(1, 4), (Y("na", "m"), Y("ma", "n"))),
             (-sp.Rational(1, 12), (Y("mm"), Y("na"), Ric("na"))),
             (sp.Rational(1, 48), (Y("mm"), Y("nn"), R)),
             (-sp.Rational(1, 24), (Y("mn"), Y("mn"), R)),
             (sp.Rational(1, 6), (Y("mn"), Y("ab"), Riem("manb"))))
    return sp.expand(sum(coefficient*tensors.contraction(factors) for coefficient, factors in terms)/tensors.mass**2)
