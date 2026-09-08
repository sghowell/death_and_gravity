"""Literal Ruf--Steinwachs (115) basis; fourth order is a WITHHELD control.

Use kernel.pole for accepted results. The literal fourth-order transcription
fails the independent curved scalar-mass test and is never used there.
"""
from functools import lru_cache

import sympy as sp

from . import tensors


def Y(indices, derivatives=""):
    return ("Y", derivatives, indices)


def Ric(indices):
    return ("Ric", "", indices)


def Riem(indices):
    return ("Riem", "", indices)


R = ("R", "", "")

ZERO = ((-sp.Rational(1, 8), (Y("mn"), Y("mn"))),
        (-sp.Rational(1, 16), (Y("mm"), Y("nn"))))

# Delta=-nabla^a nabla_a. Its minus signs are already incorporated
# below, while every ordered covariant derivative is retained.
SECOND = (
    (-2, (R, Y("mn"), Y("mn"))),
    (12, (Ric("mn"), Y("mr"), Y("nr"))),
    (1, (R, Y("mm"), Y("nn"))),
    (-4, (Ric("mn"), Y("mn"), Y("rr"))),
    (-4, (Riem("mrns"), Y("mn"), Y("rs"))),
    (12, (Y("mn"), Y("mr", "nr"))),
    (4, (Y("tt"), Y("nr", "rn"))),
    (-2, (Y("mn"), Y("mn", "aa"))),
    (1, (Y("tt"), Y("rr", "aa"))))

FOURTH = (
    (-2, (Ric("mn"), Ric("mn"), Y("aa"), Y("bb"))),
    (-1, (R, R, Y("aa"), Y("bb"))),
    (8, (Ric("mn"), R, Y("aa"), Y("mn"))),
    (-2, (R, R, Y("mn"), Y("mn"))),
    (-16, (Ric("an"), Riem("nmrs"), Y("ar"), Y("ms"))),
    (-16, (Ric("mr"), Ric("mn"), Y("aa"), Y("nr"))),
    (8, (Ric("mn"), R, Y("mr"), Y("nr"))),
    (32, (Ric("mn"), Riem("mrns"), Y("aa"), Y("rs"))),
    (24, (Ric("mr"), Ric("ns"), Y("mn"), Y("rs"))),
    (-8, (Ric("mn"), Ric("rs"), Y("mn"), Y("rs"))),
    (-8, (R, Riem("mrns"), Y("mn"), Y("rs"))),
    (-16, (Riem("manb"), Riem("rasb"), Y("mn"), Y("rs"))),
    (-8, (Ric("mr"), Ric("mn"), Y("ns"), Y("rs"))),
    (-4, (Ric("mn"), Ric("mn"), Y("rs"), Y("rs"))),
    (-16, (Riem("mrns"), Y("mn"), Y("rs", "aa"))),
    (32, (Riem("msna"), Y("mn"), Y("rs", "ar"))),
    (32, (Riem("nrsa"), Y("mn"), Y("rs", "ma"))),
    (8, (Ric("nr"), Y("nr"), Y("ss", "aa"))),
    (16, (Ric("ns"), Y("nr"), Y("sr", "aa"))),
    (16, (Riem("mrns"), Y("aa"), Y("mn", "sr"))),
    (16, (Ric("sn"), Y("rm"), Y("sr", "mn"))),
    (-32, (Ric("ns"), Y("rm"), Y("ns", "mr"))),
    (2, (R, Y("aa", "m"), Y("bb", "m"))),
    (16, (Ric("nr"), Y("nr", "m"), Y("aa", "m"))),
    (24, (R, Y("aa", "m"), Y("mn", "n"))),
    (-16, (Ric("sn"), Y("sr", "m"), Y("rm", "n"))),
    (-4, (R, Y("mn"), Y("aa", "nm"))),
    (-12, (R, Y("aa"), Y("mn", "nm"))),
    (-16, (Ric("sn"), Y("sr"), Y("rm", "nm"))),
    (16, (Riem("mnrs"), Y("aa", "m"), Y("nr", "s"))),
    (-8, (Ric("mr"), Y("nr"), Y("aa", "nm"))),
    (16, (R, Y("mn"), Y("mr", "nr"))),
    (-40, (Ric("mn"), Y("aa"), Y("mr", "nr"))),
    (32, (Riem("mrns"), Y("mn", "a"), Y("ar", "s"))),
    (-80, (Ric("sn"), Y("rm"), Y("rm", "ns"))),
    (4, (Ric("mn"), Y("aa", "m"), Y("bb", "n"))),
    (8, (R, Y("mr", "n"), Y("mn", "r"))),
    (24, (Ric("sm"), Y("sr", "n"), Y("mn", "r"))),
    (40, (Ric("ns"), Y("nr"), Y("sm", "rm"))),
    (24, (R, Y("mn"), Y("mn", "aa"))),
    (16, (Ric("mn"), Y("rr"), Y("mn", "aa"))),
    (-16, (Ric("nr"), Y("sm", "m"), Y("nr", "s"))),
    (16, (Ric("mn"), Y("aa", "m"), Y("nr", "r"))),
    (-16, (Y("rr"), Y("mn", "aanm"))),
    (-4, (Y("mn"), Y("mn", "aabb"))),
    (-16, (Y("mn"), Y("rs", "nmsr"))),
    (8, (Y("mn"), Y("mr", "naar"))),
    (-2, (Y("rr"), Y("ss", "aabb"))))


@lru_cache(maxsize=None, typed=True)
def evaluate(order):
    if type(order) is not int or order not in (0, 2, 4):
        raise ValueError("Require native covariant derivative order 0,2 or 4")
    terms = {0: ZERO, 2: SECOND, 4: FOURTH}[order]
    return sp.expand(sum(coefficient*tensors.contraction(factors) for coefficient, factors in terms))
