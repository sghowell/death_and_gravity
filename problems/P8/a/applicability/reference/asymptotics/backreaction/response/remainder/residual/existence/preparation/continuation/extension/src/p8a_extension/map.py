"""The same forced A11 map, with the A13 Einstein/log block resummed."""

import sympy as sp
from p8a_continuation.decomposition import full_remainder
from p8a_preparation.actual_map import weighted_source


def full_forced_remainder(baseline_residual, rolling_einstein_history,
                          anomaly_history, auxiliary_product_difference,
                          nonlinear_wick_derivative, source_difference, delta):
    """Retain the complete conserved-source difference in addition to A13."""
    return full_remainder(baseline_residual, rolling_einstein_history,
                          anomaly_history, auxiliary_product_difference,
                          nonlinear_wick_derivative, delta)+sp.sympify(source_difference)


def source_primitive(coefficient, hubble, initial_ratio, history):
    """Exact source part of P=a^2 q'; history=I[c(1+u/h^2)]."""
    c, h, initial, integral = map(sp.sympify, (
        coefficient, hubble, initial_ratio, history))
    return 8*(c/h-initial-integral)


def identities():
    x = sp.Symbol("x", real=True)
    a = sp.Function("a", positive=True)(x)
    c = sp.Function("c")(x)
    h, u = sp.diff(a, x)/a, -sp.diff(a, x, 2)/a
    source = weighted_source(a, c, x)
    s = sp.Symbol("s", real=True)
    history = sp.Integral((c*(1+u/h**2)).subs(x, s), (s, 0, x))
    primitive = source_primitive(c, h, c.subs(x, 0)/h.subs(x, 0), history)
    base, rolling, anomaly, aux, nonlinear, source_delta, delta = sp.symbols(
        "base rolling anomaly aux nonlinear source_delta delta", real=True)
    return {
        "source_conservation_on_each_unknown_metric": sp.simplify(
            sp.diff(source["density"], x)+3*h*(source["density"]+source["pressure"])),
        "source_primitive_with_both_endpoints": sp.simplify(sp.diff(primitive, x)-8*sp.diff(c, x)/h),
        "source_primitive_common_initial_value": sp.simplify(primitive.subs(x, 0).doit()),
        "full_A13_remainder_plus_unchanged_source": sp.simplify(
            full_forced_remainder(base, rolling, anomaly, aux, nonlinear, source_delta, delta)
            -(base+rolling/(60*delta)-anomaly+2*aux-nonlinear+source_delta)),
    }
