"""Physical normalization of the prepared source-free subspace.

The closed g ODE is a two-dimensional solution-space statement at fixed K.
It is not an off-shell action for arbitrary sources or a propagation-cone test.
"""

from functools import cache

import sympy as sp

from . import analytic, majorants
from .exact import (
    DELTA_MAX,
    DELTA_RADIUS,
    RADIUS,
    SLICE,
    nonnegative,
    number,
    parameters,
)


def center_remainder(delta, momentum_squared=1):
    delta, _ = parameters(delta, momentum_squared)
    return 48001*delta**2


def endpoint_inclusion_error(delta, momentum_squared=1):
    delta, _ = parameters(delta, momentum_squared, transfer=True)
    return 200*delta


def fixed_light_error(delta, amplitude=1, momentum_squared=1):
    delta, _ = parameters(delta, momentum_squared, transfer=True)
    return 8600*delta*nonnegative(amplitude, "amplitude")


def fixed_source_error(delta, target_norm, source_l1, momentum_squared=1):
    delta, _ = parameters(delta, momentum_squared, transfer=True)
    return delta*(8600*nonnegative(target_norm, "target norm")
                  +126000000*nonnegative(source_l1, "source L1"))


def relative_jet_remainder(delta, coordinate, parity, momentum_squared=1):
    """Full Q minus its displayed leading inner polynomial, not q alone."""
    delta, _ = parameters(delta, momentum_squared)
    coordinate = number(coordinate, "u")
    if (SLICE-abs(coordinate)).is_nonnegative is not True or parity not in ("even", "odd"):
        raise ValueError("Require |u|<=1/100 and parity 'even' or 'odd'")
    t_squared = sp.Max(delta/DELTA_RADIUS, coordinate**2/RADIUS**2)
    d_norm = majorants.calibration()["D"]
    if parity == "odd":
        return d_norm*RADIUS/5*t_squared**sp.Rational(5, 2)
    return 3*DELTA_RADIUS*d_norm*t_squared**3


def closed_equation(g_even, g_odd, coordinate):
    wronskian = g_even*sp.diff(g_odd, coordinate)-sp.diff(g_even, coordinate)*g_odd
    return {"W_g": wronskian, "F": -sp.diff(wronskian, coordinate)/wronskian,
            "G": (sp.diff(g_even, coordinate)*sp.diff(g_odd, coordinate, 2)
                  -sp.diff(g_even, coordinate, 2)*sp.diff(g_odd, coordinate))/wronskian}


@cache
def calibration():
    delta, v, radius = DELTA_MAX, DELTA_RADIUS, SLICE
    a0 = sp.Rational(7, 2640)
    # Cauchy coefficient-l1 bounds; the center has no delta^0 coefficient.
    omega_error = 3*delta**3/5
    g_even_lower, g_oddprime_lower = sp.Rational(79, 100), sp.Rational(79, 100)
    g_evenprime_upper, g_odd_upper = sp.Rational(27, 100), sp.Rational(11, 1000)
    w_lower = g_even_lower*g_oddprime_lower-g_evenprime_upper*g_odd_upper
    w_upper = sp.Rational(501, 500)*sp.Rational(101, 100)+g_evenprime_upper*g_odd_upper
    # A safe upper bound on (160-B)/delta, retaining every denominator term.
    b_loss = (104+(20+1280*a0)*delta+(640*a0+384000)*delta**2+192000*delta**3)
    even = {"l": sp.Rational(2, 5), "lp": 80, "Q": sp.Rational(9, 2500),
            "Qp": sp.Rational(74403, 100000)}
    odd = {"l": sp.Rational(1, 250), "lp": sp.Rational(6, 5),
           "Q": sp.Rational(3, 625), "Qp": sp.Rational(3201, 5000)}
    def endpoint_column(c):
        # l1 controls the Euclidean norm, and ||Psi^-1||<2.
        return 2*(c["l"]+c["lp"]+10*c["Q"]+(radius*c["Qp"]+c["Q"]/2)/(3*sp.Rational(1, 10)))
    cp, dp = 2+delta, 1+radius**2
    return {"delta_max": delta, "v": v, "omega_error": omega_error,
            "ag_squared_lower": 8/(cp*dp**12+8), "ag_squared_upper": 4*cp*dp**6/10,
            "w1_upper": cp*dp**12/10, "D_pointwise": cp*dp**4-2,
            "Dprime_pointwise": 8*cp*radius*dp**3,
            "D_delta_coefficient": dp**4, "Dprime_delta_coefficient": 8*radius*dp**3,
            "g_even_lower": g_even_lower, "g_oddprime_lower": g_oddprime_lower,
            "g_evenprime_upper": g_evenprime_upper, "g_odd_upper": g_odd_upper,
            "W_lower": w_lower, "W_upper": w_upper,
            "kinetic_lower": sp.Rational(9, 20), "kinetic_upper": sp.Rational(11, 6),
            "q_even_K_delta_coefficient": -a0, "q_even_K_remainder_coefficient": 300,
            "B_loss_per_delta": b_loss, "center_linear_coefficient": sp.Rational(14, 33),
            "center_remainder_coefficient": 48000+105*a0,
            "endpoint_even_column": endpoint_column(even), "endpoint_odd_column": endpoint_column(odd),
            "endpoint_matrix_coefficient": 200, "old_full_transfer_norm": 42,
            "fixed_light_coefficient": 8600, "fixed_source_L1_coefficient": 126000000}


@cache
def checks():
    d = analytic.derive()
    u, delta = d["u"], d["delta"]
    qe, qo1 = sp.symbols("qe qo1", real=True)
    q_at_origin = delta*qe
    ag0 = (1/d["f_sum"]).subs(u, 0)
    bg0 = (-d["w2"]/d["f_relative"]).subs(u, 0)
    g0 = ag0+bg0*q_at_origin
    fr0 = d["f_relative"].subs(u, 0)
    mu0 = 128/(delta*(delta+2))
    multiplier = 32*(delta+10)/((delta+2)*(sp.sqrt(1+delta/2)-2*delta*qe))
    out = {"actual_g_center_multiplier": sp.simplify(mu0*delta/(fr0*g0)-multiplier),
           "symplectic_center": sp.expand(1+(delta*qe)*(delta*qo1)-1-delta**2*qe*qo1),
           "center_coefficient": 160*sp.Rational(7, 2640)-sp.Rational(14, 33)}
    # The Wronskian construction solves the exact two-dimensional ODE.
    ge, go = sp.Function("ge")(u), sp.Function("go")(u)
    closed = closed_equation(ge, go, u)
    for key, field in (("even", ge), ("odd", go)):
        out[f"closed_physical_g_{key}"] = sp.simplify(
            sp.diff(field, u, 2)+closed["F"]*sp.diff(field, u)+closed["G"]*field)
    locked = (2+delta)*(5+2*delta)/(10+delta)
    out["locked_center_distinct"] = sp.factor(
        locked-1-4*delta/5-sp.Rational(6, 5)*delta**2/(10+delta))
    out["not_same_linear_coefficient"] = sp.Rational(4, 5)-sp.Rational(14, 33)-sp.Rational(62, 165)
    c = calibration()
    margins = {"symplectic_positive": sp.Rational(1, 10)-c["omega_error"],
               "ag_lower": c["ag_squared_lower"]-sp.Rational(16, 25),
               "ag_upper": 1-c["ag_squared_upper"], "w1_half": sp.Rational(1, 2)-c["w1_upper"],
               "D_pointwise": sp.Rational(1, 1000)-c["D_pointwise"],
               "Dprime_pointwise": sp.Rational(1, 5)-c["Dprime_pointwise"],
               "D_delta": 2-c["D_delta_coefficient"],
               "Dprime_delta": sp.Rational(1, 10)-c["Dprime_delta_coefficient"],
               "W_lower": c["W_lower"]-sp.Rational(3, 5), "W_upper": 2-c["W_upper"],
               "q_even_sign": sp.Rational(7, 2640)-300*DELTA_MAX,
               "B_loss": 105-c["B_loss_per_delta"],
               "center_remainder": 48001-c["center_remainder_coefficient"],
               "endpoint_even": 162-c["endpoint_even_column"],
               "endpoint_odd": 3-c["endpoint_odd_column"],
               "endpoint_matrix": 200-c["endpoint_even_column"]-c["endpoint_odd_column"]}
    if any(value.is_positive is not True for value in margins.values()):
        raise ValueError("A physical normalization, coefficient or endpoint margin failed")
    return {"residuals": out, "strict_margins": margins}
