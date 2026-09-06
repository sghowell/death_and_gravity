"""Exact full-equation decomposition, retaining every rolling/state remainder.

The constant-coefficient block is a preconditioner, not the Frechet derivative
of the rolling prepared state's full semiclassical equation.
"""

import sympy as sp
from p8a_preparation.actual_map import local_coefficient


def anomaly_polynomial(hubble, potential):
    h, u = map(sp.sympify, (hubble, potential))
    return u**2/4+h**2*(u+h**2)/30


def full_trace_in_s(scale, wick, delta, time):
    """The actual unforced trace residual, multiplied by a^2.

This is a symbolic identity interface, not a numerical state constructor.
Wick means the dimensionless S of A.11, with its named finite prescription.
    """
    a, s, delta = map(sp.sympify, (scale, wick, delta))
    h, u = sp.diff(a, time)/a, -sp.diff(a, time, 2)/a
    return (sp.diff(s, time, 2)-2*h*sp.diff(s, time)+2*(u+h**2)*s
            -a**2*u/(60*delta)+anomaly_polynomial(h, u))


def full_remainder(baseline_residual, rolling_einstein_history,
                   anomaly_history, auxiliary_product_difference,
                   nonlinear_wick_derivative, delta):
    """Return G_f in (L_f-c I^2/2)X=G_f, conditional on exact inputs.

rolling_einstein_history is integral[(a^2-af^2)u-(ab^2-af^2)ub].
anomaly_history is integral(P[h,u]-P[hb,ub]). The auxiliary product is
a^2*h*q-ab^2*hb*qb. nonlinear_wick_derivative is derivative of
R[u]-R[ub]+(d(a)-df)u-(d(ab)-df)ub, including the actual shared history.
No one of these inputs is assumed small or omitted by this interface.
    """
    base, einstein, anomaly, auxiliary, nonlinear, delta = map(sp.sympify, (
        baseline_residual, rolling_einstein_history, anomaly_history,
        auxiliary_product_difference, nonlinear_wick_derivative, delta))
    return base+einstein/(60*delta)-anomaly+2*auxiliary-nonlinear


def identities():
    x = sp.Symbol("x", real=True)
    delta, af = sp.symbols("delta af", positive=True)
    a, ss = sp.Function("a", positive=True)(x), sp.Function("S")(x)
    h, u = sp.diff(a, x)/a, -sp.diff(a, x, 2)/a
    q = ss/a**2
    trace_q = (sp.diff(q, x, 2)+2*h*sp.diff(q, x)-u/(60*delta)
               +anomaly_polynomial(h, u)/a**2)
    d, db, df, ubp, xx, rp, rpb, hu, hbub = sp.symbols(
        "d db df ubprime X Rprime Rbprime hu hbub", real=True)
    iau, i2x, ip, aux, residual = sp.symbols("Iau I2X IP aux residual", real=True)
    # A.11's actual, integrated full map before resumming the fixed block.
    original = (residual+iau/(60*delta)-ip+2*aux-(rp-rpb)
                -d*xx-(d-db)*ubp+hu/2-hbub/2)
    nonlinear = rp-rpb+(d-df)*xx+(d-db)*ubp-hu/2+hbub/2
    retained = full_remainder(residual, iau-af**2*i2x, ip, aux, nonlinear, delta)
    w = sp.Function("W")(x)
    t = sp.Symbol("t", real=True)
    i2 = sp.Integral((x-t)*sp.diff(w, x).subs(x, t), (t, 0, x))
    return {
        "full_unforced_trace_in_S": sp.simplify(a**2*trace_q-full_trace_in_s(a, ss, delta, x)),
        "named_local_coefficient_not_changed": sp.simplify(
            local_coefficient(af)+sp.Rational(19, 60)+sp.log(af/2)/2),
        "full_integrated_map_without_dropped_terms": sp.expand(
            original+df*xx-af**2*i2x/(60*delta)-retained),
        "double_primitive_derivative": sp.simplify(sp.diff(i2, x)-w+w.subs(x, 0)),
        "double_primitive_initial_value": sp.simplify(i2.subs(x, 0).doit()),
    }
