"""Actual coefficient derivative envelopes and polynomial tail majorants."""
from functools import cache

import sympy as sp
from p8_vector_clock_matching import continuation
from p8_vector_state import comparison, wkb
from p8_vector_subtraction import tail

p, r, b2, b4, aa, bb, t, si, wa, wb, c = sp.symbols(
    "P2 P4 B2 B4 residual_A residual_B inverse_frequency_squared inverse_S weight_A weight_B c1")
VARIABLES = (p, r, b2, b4, aa, bb, t, si, wa, wb, c)
RESIDUAL_CONSTANT = comparison.RESIDUAL_CONSTANT
RESIDUAL_DERIVATIVE_CONSTANT = sp.Integer(100_000_000)
TAIL_DERIVATIVE_CONSTANT = sp.Integer(100_000_000)


@cache
def coefficients(kind):
    return {name: sp.factor(value.subs(continuation.local.dimension, 3))
            for name, value in continuation.coefficients(kind).items()
            if name in ("P2", "P4", "B2", "B4", "A", "B")}


def polynomial_majorant(expression):
    if (not isinstance(expression, sp.Expr) or expression.has(sp.Float)
            or expression.has(sp.oo, -sp.oo, sp.zoo, sp.nan)):
        raise ValueError("Require an exact finite symbolic polynomial")
    poly = sp.Poly(expression, *VARIABLES)
    if any(not coefficient.is_Rational for _, coefficient in poly.terms()):
        raise ValueError("Require a fixed rational polynomial")
    return sum(abs(coefficient)*sp.prod(variable**power for variable, power in zip(VARIABLES, powers, strict=True))
               for powers, coefficient in poly.terms())


def value_and_slope(expression, envelopes):
    majorant = polynomial_majorant(expression)
    values = {variable: data[0] for variable, data in envelopes.items()}
    for value, slope in envelopes.values():
        if not isinstance(value, sp.Rational) or not isinstance(slope, sp.Rational) or value < 0 or slope < 0:
            raise ValueError("Require nonnegative exact rational value and derivative envelopes")
    if set(envelopes) != set(VARIABLES):
        raise ValueError("Require the full independent coefficient envelope dictionary")
    value = majorant.subs(values)
    slope = sum(sp.diff(majorant, variable).subs(values)*data[1] for variable, data in envelopes.items())
    return {"value_upper": value, "derivative_upper": slope}


@cache
def abstract_polynomials():
    residual = (-aa*(p+r*t)*si-2*p*r-r**2*t+bb*si
                +sp.Rational(3, 4)*(b2+t*b4)**2*si**2)
    R1 = (2*p*r-p**3+r*(r-p**2)*t)*si
    R2 = (p**2-r+p*r*t)*si
    R3 = (b4-b2*(2*p+(p**2+2*r)*t+2*p*r*t**2+r**2*t**3))*si**2
    tail = wb*R1+wa*c**2*R2+wa*c*R3+wa*(b2+b4*t)**2*si**3/4
    return {"residual_bracket": residual, "reference_tail_bracket": tail}


@cache
def envelopes(kind):
    data = coefficients(kind)
    names = {p: "P2", r: "P4", b2: "B2", b4: "B4", aa: "A", bb: "B"}
    out = {variable: (wkb.box_bound(data[name])["absolute_upper"],
                      wkb.box_bound(wkb.D(data[name]))["absolute_upper"])
           for variable, name in names.items()}
    tmax = 1/wkb.MASS_TIME_MIN**2
    sdot = out[b2][0]*tmax+out[b4][0]*tmax**2
    out.update({t: (tmax, 4*tmax), si: (sp.Integer(2), 4*sdot),
                wa: (sp.Integer(1), sp.Integer(2)),
                wb: (sp.Rational(3, 2), sp.Integer(2)),
                c: (sp.Integer(2), sp.Integer(5))})
    return out


@cache
def bounds(kind):
    data = envelopes(kind)
    polynomials = abstract_polynomials()
    residual = value_and_slope(polynomials["residual_bracket"], data)
    tail = value_and_slope(polynomials["reference_tail_bracket"], data)
    residual_derivative = residual["derivative_upper"]+8*residual["value_upper"]
    tail_derivative = tail["derivative_upper"]+16*tail["value_upper"]
    return {"coefficient_value_derivative_envelopes": {str(key): values for key, values in data.items()},
            "residual_bracket_upper": residual["value_upper"],
            "residual_time_derivative_over_inverse_frequency_fourth_upper": residual_derivative,
            "residual_time_derivative_integer_upper": sp.ceiling(residual_derivative),
            "reference_tail_bracket_upper": tail["value_upper"],
            "reference_tail_density_derivative_bracket_upper": tail_derivative,
            "reference_tail_density_derivative_integer_upper": sp.ceiling(tail_derivative)}


@cache
def checks():
    out = {}
    for kind in ("transverse", "longitudinal"):
        data = coefficients(kind)
        for name, value in data.items():
            out[kind+"_coefficient_derivative_reconstruction_"+name] = wkb.box_bound(wkb.D(value))["reconstruction"]
        out[kind+"_inverse_S_derivative_numerator"] = sp.factor(
            wkb.D(1+data["P2"]*wkb.t+data["P4"]*wkb.t**2)
            -data["B2"]*wkb.t-data["B4"]*wkb.t**2)
    return out


@cache
def replay_checks():
    out = {}
    for kind in ("transverse", "longitudinal"):
        data = coefficients(kind)
        mapping = {p: data["P2"], r: data["P4"], b2: data["B2"], b4: data["B4"],
                   aa: data["A"], bb: data["B"], t: wkb.t,
                   si: 1/(1+data["P2"]*wkb.t+data["P4"]*wkb.t**2)}
        value = wkb.t**2*abstract_polynomials()["residual_bracket"].subs(mapping, simultaneous=True)
        out[kind+"_residual_polynomial_replays_frozen_exact_reference"] = sp.factor(value-wkb.frequency(kind)["residual"])
    mapping = {tail.p: p, tail.r: r, tail.b2: b2, tail.b4: b4,
               tail.c: c, tail.wa: wa, tail.wb: wb, tail.t: t}
    value = abstract_polynomials()["reference_tail_bracket"].subs(si, 1/(1+p*t+r*t**2))
    out["tail_polynomial_replays_frozen_exact_subtraction"] = sp.factor(value-tail.algebra()["remainder"].subs(mapping, simultaneous=True))
    return out
