"""Exact physical-g loading with an explicit H0^2 source cost.

The delta=0 coefficients are used only on the punctured interval J.
Endpoint norms refer to the full Frobenius frame, not raw light columns.
"""

from functools import cache
from math import comb

import sympy as sp


def exact(value, name="value"):
    if isinstance(value, (bool, float)):
        raise TypeError(f"{name} must be exact, finite and real")
    value = sp.sympify(value)
    if (value.has(sp.Float, sp.oo, sp.zoo, sp.nan) or value.is_number is not True
            or value.is_real is not True or value.is_finite is not True):
        raise ValueError(f"{name} must be exact, finite and real")
    return value


def source_cost(target_norm=1):
    norm = exact(target_norm, "Frobenius target norm")
    if norm.is_nonnegative is not True:
        raise ValueError("The target norm must be nonnegative")
    d = calibration()
    return {key: norm*d[key] for key in
            ("source_sup", "source_L1", "source_L2", "source_second_derivative_L2")}


@cache
def derive():
    u, momentum, coordinate = sp.symbols("u K x", real=True)
    ell, radius = sp.Rational(1, 100), sp.Rational(1, 50)
    d = 1+u**2
    kg, kf, gg, gf = d**6/8, 1/(2*d**6), d**2/8, 1/(2*d**2)
    spring = 4*(1-u**2)/(d**8*(d**4-1))
    inverse = d**8*(d**4-1)/(4*(1-u**2))
    jg, jf, friction = spring/kg, spring/kf, sp.diff(kg, u)/kg
    phase = sp.Matrix([[0, 1, 0, 0],
                       [-gg*momentum/kg-jg, -friction, jg, 0],
                       [0, 0, 0, 1],
                       [jf, 0, -gf*momentum/kf-jf, friction]])
    scale = sp.diag(1, ell, 1, ell)
    scaled = (ell*scale*phase*scale.inv()).applyfunc(sp.factor)
    t2, t1, t0 = inverse*kf, inverse*sp.diff(kf, u), 1+inverse*gf*momentum
    coefficients = [
        kg*sp.diff(t0, u, 2)+sp.diff(kg, u)*sp.diff(t0, u)+gg*momentum*t0+gf*momentum,
        kg*(sp.diff(t1, u, 2)+2*sp.diff(t0, u))
        +sp.diff(kg, u)*(sp.diff(t1, u)+t0)+gg*momentum*t1+sp.diff(kf, u),
        kg*(sp.diff(t2, u, 2)+2*sp.diff(t1, u)+t0)
        +sp.diff(kg, u)*(sp.diff(t2, u)+t1)+gg*momentum*t2+kf,
        kg*(2*sp.diff(t2, u)+t1)+sp.diff(kg, u)*t2,
        kg*t2,
    ]
    coefficients = tuple(sp.factor(4*value/d**6) for value in coefficients)
    return {"u": u, "K": momentum, "x": coordinate, "ell": ell, "radius": radius,
            "left": -radius, "right": -ell, "d": d, "a": d**2,
            "K_g": kg, "K_f": kf, "G_g": gg, "G_f": gf,
            "spring": spring, "inverse_spring": inverse,
            "phase": phase, "scaled_phase": scaled, "phase_scale": scale,
            "physical_input": sp.Matrix([0, 2, 0, 0]),
            "g_from_f": (t0, t1, t2), "source_coefficients": coefficients}


@cache
def endpoint_jets():
    """Rows mapping (g,ell*g',f,ell*f') at the right end to f_x^(j)."""
    d = derive()
    u, ell, matrix = d["u"], d["ell"], d["scaled_phase"]
    derivatives = [matrix.applyfunc(lambda entry, j=j: sp.factor(
        ell**j*sp.diff(entry, u, j).subs(u, d["right"]))) for j in range(5)]
    states = [sp.eye(4)]
    for degree in range(5):
        states.append(sum((comb(degree, j)*derivatives[j]*states[degree-j]
                           for j in range(degree+1)), sp.zeros(4)))
    return tuple(row[2, :].applyfunc(sp.cancel) for row in states)


def falling(degree, order):
    return sp.factorial(degree)/sp.factorial(degree-order)


@cache
def hermite():
    """Degree11 Bernstein rows: zero left jets0..5, exact right jets0..5."""
    jets = endpoint_jets()
    rows = [sp.zeros(1, 4) for _ in range(12)]
    for j in range(6):
        rows[11-j] = sum(((-1)**k*comb(j, k)*jets[k]/falling(11, k)
                          for k in range(j+1)), sp.zeros(1, 4)).applyfunc(sp.cancel)
    derivative_rows = []
    for order in range(7):
        derivative_rows.append(tuple((falling(11, order)*sum(
            ((-1)**(order-k)*comb(order, k)*rows[index+k] for k in range(order+1)),
            sp.zeros(1, 4))).applyfunc(sp.cancel) for index in range(12-order)))
    return {"degree": 11, "rows": tuple(rows), "derivative_rows": tuple(derivative_rows)}


def hermite_polynomial(endpoint, momentum_squared=1):
    """Exact f(x), with x=(u-left)/ell, for an exact scaled physical endpoint."""
    momentum = exact(momentum_squared, "K")
    if (momentum-1).is_nonnegative is not True or (4-momentum).is_nonnegative is not True:
        raise ValueError("Require 1<=K<=4")
    if not isinstance(endpoint, (tuple, list)) or len(endpoint) != 4:
        raise ValueError("Require four scaled physical endpoint components")
    values = sp.Matrix([exact(value, "endpoint component") for value in endpoint])
    d, h = derive(), hermite()
    coordinate = d["x"]
    return sp.expand(sum(sp.binomial(11, index)*coordinate**index*(1-coordinate)**(11-index)
                         *(row.subs(d["K"], momentum)*values)[0]
                         for index, row in enumerate(h["rows"])))


def _momentum_row_bound(row):
    momentum = derive()["K"]
    return sum(sum(abs(value)*4**power[0] for power, value in
                   sp.Poly(sp.cancel(entry), momentum).terms()) for entry in row)


def _rational_bound(expression):
    """Continuous |u|<=1/50, |K|<=4 coefficient-l1 / positive-denominator bound."""
    d = derive()
    u, momentum, radius = d["u"], d["K"], d["radius"]
    numerator, denominator = sp.cancel(expression).as_numer_denom()
    upper = sum(abs(value)*radius**power[0]*4**power[1]
                for power, value in sp.Poly(numerator, u, momentum).terms())
    content, factors = sp.factor_list(denominator, u)
    lower = abs(content)
    for factor, exponent in factors:
        if sp.expand(factor-(1+u**2)) == 0:
            factor_lower = sp.Integer(1)
        elif sp.degree(factor, u) == 1:
            polynomial = sp.Poly(factor, u)
            factor_lower = abs(polynomial.nth(0))-abs(polynomial.nth(1))*radius
        else:
            raise ValueError("Unproved rational denominator factor")
        if factor_lower.is_positive is not True:
            raise ValueError("A rational denominator enclosure contains zero")
        lower *= factor_lower**exponent
    return sp.cancel(upper/lower)


@cache
def coefficient_bounds():
    d, h = derive(), hermite()
    source_raw = tuple(tuple(_rational_bound(sp.diff(value, d["u"], order))/d["ell"]**index
                             for index, value in enumerate(d["source_coefficients"]))
                       for order in range(3))
    hermite_raw = tuple(max(_momentum_row_bound(row) for row in rows)
                        for rows in h["derivative_rows"])
    return {"source_raw": source_raw,
            "source_integer": tuple(tuple(sp.ceiling(value) for value in row) for row in source_raw),
            "hermite_raw": hermite_raw,
            "hermite_integer": tuple(sp.ceiling(value) for value in hermite_raw)}


@cache
def calibration():
    d, b = derive(), coefficient_bounds()
    ell, coefficients, f_jets = d["ell"], b["source_integer"], b["hermite_integer"]
    source = sum(coefficients[0][j]*f_jets[j] for j in range(5))
    source_uu = sum(comb(2, order)*coefficients[2-order][j]*ell**(-order)*f_jets[j+order]
                    for j in range(5) for order in range(3))
    # Entrywise majorants of diag(1,ell,1,ell) C0(b) W_-(ell).
    # |ag|<1, |bg|<2, |bf|<1/2 and all their first derivatives<1.
    # sqrt(ell)=1/10 and mu<13/4. The full wave Psi norm is<2.
    frame_rows = ((1, 0, sp.Rational(1, 5), 0),
                  (ell, ell, sp.Rational(101, 1000), sp.Rational(13, 20)),
                  (1, 0, sp.Rational(1, 20), 0),
                  (ell, ell, sp.Rational(13, 500), sp.Rational(13, 80)))
    frame_square = sum(value**2 for row in frame_rows for value in row)
    g_even, g_odd, volume = sp.Rational(101, 100), sp.Rational(21, 1000), sp.Rational(101, 100)
    dual = (volume*g_odd/2, volume*g_even/2)
    return {"ell": ell, "source_coefficient_bounds": coefficients, "hermite_derivative_bounds": f_jets,
            "source_sup_per_scaled_endpoint_linf": source,
            "source_uu_sup_per_scaled_endpoint_linf": source_uu,
            "physical_frame_entry_majorants": frame_rows, "physical_frame_Frobenius_squared": frame_square,
            "full_wave_norm": sp.Integer(2), "scaled_endpoint_norm_factor": sp.Integer(4),
            "source_sup": sp.Integer(2*10**9), "source_L1": sp.Integer(2*10**7),
            "source_L2": sp.Integer(2*10**8), "source_second_derivative_L2": sp.Integer(4*10**14),
            "analytic_gramian_lower": sp.Rational(1, 4*10**16),
            "light_dual_kernel_suprema": dual, "regular_light_L2_lower": sp.Integer(19),
            "regular_even_L2_lower": sp.Integer(900),
            "source_normalization": "sigma=tau² Pi/M²; S_probe=(1/2) integral dT a³ Pi gamma_g"}


@cache
def checks():
    d, h, jets, c = derive(), hermite(), endpoint_jets(), calibration()
    u = d["u"]
    field = sp.Function("f")(u)
    g = sum(value*sp.diff(field, u, order) for order, value in enumerate(d["g_from_f"]))
    f_equation = sp.diff(d["K_f"]*sp.diff(field, u), u)+d["G_f"]*d["K"]*field+d["spring"]*(field-g)
    direct = 4/d["a"]**3*(sp.diff(d["K_g"]*sp.diff(g, u), u)
                                 +d["G_g"]*d["K"]*g+d["spring"]*(g-field))
    assembled = sum(value*sp.diff(field, u, order) for order, value in enumerate(d["source_coefficients"]))
    residuals = {"inverse_spring": sp.factor(d["spring"]*d["inverse_spring"]-1),
                 "physical_source_B2": sp.factor(d["a"]**3/(4*d["K_g"])-2),
                 "literal_unsourced_f": sp.factor(f_equation),
                 "literal_g_source_operator": sp.factor(direct-assembled)}
    for order in range(6):
        left = falling(11, order)*sum(((-1)**(order-k)*comb(order, k)*h["rows"][k]
                                      for k in range(order+1)), sp.zeros(1, 4))
        right = falling(11, order)*sum(((-1)**(order-k)*comb(order, k)*h["rows"][11-order+k]
                                       for k in range(order+1)), sp.zeros(1, 4))
        for index in range(4):
            residuals[f"Hermite_left_{order}_{index}"] = sp.cancel(left[index])
            residuals[f"Hermite_right_{order}_{index}"] = sp.cancel(right[index]-jets[order][index])
    for column in range(4):
        value = sum(coefficient.subs(u, d["right"])*d["ell"]**(-j)*jets[j][column]
                    for j, coefficient in enumerate(d["source_coefficients"]))
        derivative = sum(sp.diff(coefficient, u).subs(u, d["right"])*d["ell"]**(-j)*jets[j][column]
                         +coefficient.subs(u, d["right"])*d["ell"]**(-j-1)*jets[j+1][column]
                         for j, coefficient in enumerate(d["source_coefficients"]))
        residuals[f"target_source_value_zero_{column}"] = sp.cancel(value)
        residuals[f"target_source_derivative_zero_{column}"] = sp.cancel(derivative)
    primitive = sp.Rational(2501, 2500)
    ratio = sp.Rational(2, 5)
    le, lo = 1+sp.Rational(1, 40)*ratio**2, sp.Rational(1, 50)+sp.Rational(1, 800)*ratio**3
    qe, qo, pole = sp.Rational(3, 400)*ratio**2, sp.Rational(1, 100)*ratio, sp.Rational(1, 300)
    margins = {"source_sup": c["source_sup"]-4*c["source_sup_per_scaled_endpoint_linf"],
               "source_L1": c["source_L1"]-4*c["ell"]*c["source_sup_per_scaled_endpoint_linf"],
               "source_L2": c["source_L2"]-4*sp.sqrt(c["ell"])*c["source_sup_per_scaled_endpoint_linf"],
               "source_second_derivative_L2": c["source_second_derivative_L2"]
               -4*sp.sqrt(c["ell"])*c["source_uu_sup_per_scaled_endpoint_linf"],
               "scaled_physical_frame_norm": 4-c["physical_frame_Frobenius_squared"],
               "physical_bg_squared": 4-sp.Rational(16, 5),
               "physical_bf_squared": sp.Rational(1, 4)-primitive**18/5,
               "physical_ag_squared": 1-4*primitive**6/5,
               "physical_first_map_derivatives": 1-36*c["ell"]*2,
               "prepared_D_on_J": pole-2*(primitive**4-1),
               "prepared_volume_on_J": sp.Rational(101, 100)-primitive**6,
               "prepared_g_even_on_J": sp.Rational(101, 100)-le-2*pole*qe,
               "prepared_g_odd_on_J": sp.Rational(21, 1000)-lo-2*pole*qo,
               "regular_light_L2_lower": 1-c["regular_light_L2_lower"]**2*c["ell"]
               *sum(value**2 for value in c["light_dual_kernel_suprema"]),
               "regular_even_L2_lower": 1-c["regular_even_L2_lower"]
               *sp.sqrt(c["ell"])*c["light_dual_kernel_suprema"][0]}
    if any(value != 0 for value in residuals.values()):
        raise ValueError("An exact Hermite/source identity failed")
    if any(value.is_positive is not True for value in margins.values()):
        raise ValueError("A continuous source-cost margin failed")
    return {"residuals": residuals, "strict_margins": margins}


@cache
def parent_checks():
    from p8_variable_response import operator, source

    d, old, probe = derive(), operator.derive(), source.coefficients()
    substitutions = {old["c"]: 2, old["u"]: d["u"]}
    values = {"literal_K_g": old["K1"], "literal_K_f": old["K2"],
              "literal_spring": old["mass"]*old["K_relative"]}
    expected = {"literal_K_g": d["K_g"], "literal_K_f": d["K_f"], "literal_spring": d["spring"]}
    residuals = {key: sp.factor(value.subs(substitutions)-expected[key]) for key, value in values.items()}
    raw = sp.Matrix([0, probe["j_light"], 0, probe["j_heavy"]])
    physical = (operator.physical_map()*raw).subs(substitutions)-d["physical_input"]
    residuals.update({f"physical_source_map_{j}": sp.simplify(value) for j, value in enumerate(physical)})
    if any(value != 0 for value in residuals.values()):
        raise ValueError("The literal pinned physical action/source dictionary changed")
    return residuals


@cache
def independent_checks():
    from . import construction_independent as independent

    d, h, jets, calibration_primary = derive(), hermite(), endpoint_jets(), calibration()
    def polynomial(values):
        return sum(sp.Rational(value)*d["u"]**i*d["K"]**j for (i, j), value in values.items())

    comparisons = 0
    for index, value in enumerate(independent.source_coefficients()):
        for order in range(3):
            other = value.diff(order)
            expression = polynomial(other.numerator)/((1+d["u"]**2)**other.d_power*(1-d["u"]**2)**other.e_power)
            if sp.factor(sp.diff(d["source_coefficients"][index], d["u"], order)-expression) != 0:
                raise ValueError("Independent differential-operator composition failed")
            comparisons += 1
    for actual, expected in zip(jets, independent.endpoint_jets()):
        for value, other in zip(actual, expected):
            if sp.cancel(value-polynomial(other)) != 0:
                raise ValueError("Independent Taylor endpoint derivative failed")
            comparisons += 1
    for rows, expected in zip(h["derivative_rows"], independent.hermite()["derivative_rows"]):
        for row, other in zip(rows, expected):
            for value, poly in zip(row, other):
                if sp.cancel(value-polynomial(poly)) != 0:
                    raise ValueError("Independent power-basis Hermite solve/Bernstein conversion failed")
                comparisons += 1
    independent_calibration = independent.calibration()
    for key in ("source_coefficient_bounds", "hermite_derivative_bounds",
                "source_sup_per_scaled_endpoint_linf", "source_uu_sup_per_scaled_endpoint_linf"):
        if calibration_primary[key] != independent_calibration[key]:
            raise ValueError("Independent continuous coefficient/source bound failed")
    return {"exact_scalar_function_or_polynomial_comparisons": comparisons,
            "independent_source_coefficient_bounds": independent_calibration["source_coefficient_bounds"],
            "independent_hermite_derivative_bounds": independent_calibration["hermite_derivative_bounds"],
            "independent_source_L2_raw": sp.Rational(independent_calibration["source_L2_raw"]),
            "independent_source_second_derivative_L2_raw": sp.Rational(independent_calibration["source_uu_L2_raw"]),
            "all_equal": True}
