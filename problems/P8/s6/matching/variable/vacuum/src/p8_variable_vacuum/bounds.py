"""Continuous rational sign margins, independent of a time/lapse grid."""

from functools import cache

import sympy as sp


@cache
def derive():
    c, z = sp.symbols("c z", real=True)
    dmax = sp.Rational(101, 100)
    Q = (87*c**2-525*c+243)/50
    margin = -sp.Rational(459, 50)-Q
    shifted = sp.Poly(sp.expand(margin.subs(c, 2+2*z)), z)
    power = [shifted.nth(j) for j in range(3)]
    bernstein = [power[0], power[0]+power[1]/2, sum(power)]
    return {"c": c, "Q": Q, "endpoint_margin": margin,
            "endpoint_margin_factor": 3*(c-2)*(117-29*c)/50,
            "endpoint_margin_Bernstein": bernstein,
            "strict_rational_margins": {
                "d12_below8over7": sp.Rational(8, 7)-dmax**12,
                "hprime_lower_above3": 4*sp.Rational(99, 100)/dmax**2-3,
                "root_cube_below_9over4_cubed": sp.Rational(9, 4)**3-sp.Rational(200, 23),
                "theta_below_one": 1-sp.Rational(2, 25),
                "linear_factor_117minus29c_lower": sp.S.One,
                "c1_g_U_magnitude": sp.Rational(123, 25),
                "positive_branch_g_U_magnitude": sp.Rational(459, 800)},
            "constants": {"h_squared_upper": sp.Rational(4, 25), "hprime_strict_lower": 3,
                "y_cubed_strict_lower": 7, "y_upper": 2, "theta_strict_upper_positive_branch": sp.Rational(2, 25),
                "root_strict_upper_positive_branch": sp.Rational(9, 4),
                "g_U_upper_c1": -sp.Rational(123, 25), "g_U_upper_positive_branch": -sp.Rational(459, 800)}}


def checks():
    d = derive()
    if any(value.is_positive is not True for value in d["strict_rational_margins"].values()):
        raise ValueError("A strict continuous rational margin failed")
    if any(value < 0 for value in d["endpoint_margin_Bernstein"]):
        raise ValueError("The continuous lapse polynomial bound failed")
    c = d["c"]
    values = {
        "continuous_lapse_factorization": sp.factor(d["endpoint_margin"]-d["endpoint_margin_factor"]),
        "positive_branch_bound_polynomial": sp.expand(d["Q"]-((sp.Rational(6, 25)+sp.Rational(3, 2))*c**2
                                                           -sp.Rational(21, 2)*c+sp.Rational(243, 50))),
        "c1_bound_arithmetic": sp.Rational(6, 25)*(1+16)-9+sp.Rational(123, 25),
        "positive_branch_final_division_at_c4": -sp.Rational(459, 50)/16+sp.Rational(459, 800),
    }
    if any(value != 0 for value in values.values()):
        raise ValueError("A continuous bound identity failed")
    return values
