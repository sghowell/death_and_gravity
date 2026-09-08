"""Chosen covariant pole continuation and its first dimensional Taylor jet."""
from functools import lru_cache

import sympy as sp

from . import bimetric, contractions, geometry


@lru_cache(maxsize=None, typed=True)
def pole(order):
    if type(order) is not int or order not in (2, 4):
        raise ValueError("Require native derivative order 2 or 4; flat potential is not re-added")
    return contractions.second_pole() if order == 2 else bimetric.fourier_bilinear()


@lru_cache(maxsize=None, typed=True)
def dimension_jet(order, jet):
    if type(order) is not int or order not in (2, 4) or type(jet) is not int or jet not in (0, 1):
        raise ValueError("Require derivative order 2/4 and dimensional Taylor jet 0/1")
    return sp.expand(sp.diff(pole(order), geometry.dimension, jet).subs(geometry.dimension, 3))
