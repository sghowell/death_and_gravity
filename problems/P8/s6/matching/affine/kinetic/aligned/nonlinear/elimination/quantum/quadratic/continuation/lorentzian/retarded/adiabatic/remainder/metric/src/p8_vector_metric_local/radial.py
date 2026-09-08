"""Full metric Laurent coefficients and the common-dimensional local limit."""
from functools import cache

import sympy as sp
from p8_vector_dimensional import local

from . import counterterms, jets, variation


def laurent(output, order):
    output, order = jets.output(output), jets.order(order)
    return _laurent(output, order)


@cache
def _laurent(output, order):
    angular = local.radial_polynomial(variation.combined(output, order), order)
    factor = sp.Rational((-1)**(2-order), sp.factorial(2-order))
    pole = sp.expand(factor*angular.subs(jets.D, 3))
    evanescent = sp.expand(-2*factor*sp.diff(angular, jets.D).subs(jets.D, 3))
    return {"pole": pole, "radial_evanescent": evanescent,
            "finite_MSbar_mu_m": sp.expand(sp.harmonic(2-order)*pole+evanescent)}


def matched(output, order):
    output, order = jets.output(output), jets.order(order)
    return _matched(output, order)


@cache
def _matched(output, order):
    raw = counterterms.mass**(4-2*order)*jets.actual(laurent(output, order)["finite_MSbar_mu_m"])
    return jets.linear_clean(raw+2*counterterms.operator(output, order, 1))


@cache
def pole_checks():
    out = {}
    for output in ("N", "Z"):
        for order in (0, 1, 2):
            raw = counterterms.mass**(4-2*order)*jets.actual(laurent(output, order)["pole"])
            out[output+"_full_metric_pole_"+str(2*order)] = jets.linear_clean(raw-counterterms.operator(output, order))
    return out
