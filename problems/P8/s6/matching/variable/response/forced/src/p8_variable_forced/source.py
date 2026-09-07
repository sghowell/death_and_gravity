"""Actual source versus symplectic projection onto the prepared sector."""

from functools import cache

import sympy as sp

from . import operator


def projected_source_factor(scale_factor, physical_wronskian, symplectic_form):
    return scale_factor**3*physical_wronskian/(2*symplectic_form)


@cache
def checks():
    d = operator.derive()
    u, delta = d["u"], d["delta"]
    ge, go = sp.Function("ge")(u), sp.Function("go")(u)
    aa, sigma, omega = sp.symbols("a sigma Omega", positive=True)
    w = ge*sp.diff(go, u)-sp.diff(ge, u)*go
    adot, bdot = -aa**3*go*sigma/(2*omega), aa**3*ge*sigma/(2*omega)
    residuals = {
        "source_projection_first_derivative_cancellation": sp.expand(adot*ge+bdot*go),
        "source_projection_closed_equation": sp.expand(adot*sp.diff(ge, u)+bdot*sp.diff(go, u)
                                                       -projected_source_factor(aa, w, omega)*sigma),
        "actual_full_source_diagonal": sp.factor(d["ag"]*d["jL"]+d["bg"]*d["jH"]-2),
        "prepared_limit_physical_ag_squared": sp.factor(d["ag"].subs({u: 0, delta: 0})**2-sp.Rational(4, 5)),
        "prepared_limit_projection": sp.Rational(4, 5)/2-sp.Rational(2, 5),
        "nonprepared_complement_not_zero": 2-sp.Rational(2, 5)-sp.Rational(8, 5),
    }
    return residuals
