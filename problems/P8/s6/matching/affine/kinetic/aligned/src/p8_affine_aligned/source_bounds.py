"""Uniform nonlinear source bounds from the unchanged lower-coefficient ODE.

These bounds control the specified source and its spatial curl, not a
heavy-field inverse, omitted-action remainder, loop error or EFT cutoff.
"""
from functools import cache

import sympy as sp
from p8_affine import dictionary, lower
from p8_affine_retuned import bounds as old_bounds

from . import alignment, modes


@cache
def ode_bridge():
    u, x = dictionary.u, dictionary.x
    data, equation = dictionary.lift(), lower.ode()
    h, p = data["h"], data["p"]
    coefficient = 1/(2*x)+3/(16*h*p**2)
    forcing = -3*sp.diff(h, u)*(1+x)/(16*h**3*p**2)
    return {"coefficient": coefficient, "forcing": forcing,
            "coefficient_residual": sp.factor(equation["coefficient"]-coefficient),
            "forcing_residual": sp.factor(equation["forcing"]-forcing),
            "basepoint": equation["basepoint"],
            "forcing_zero_at_basepoint": sp.factor(forcing.subs(x, -1))}


@cache
def electric_identity():
    s, h = sp.symbols("s h", positive=True)
    DK, Q, Qx, Ds, DDK = sp.symbols("Delta_K Q Q_x D_s D_Delta_K", real=True)
    delta = (s**2-1)/h
    # Coordinate source is S0=Snormal/s. The physical electric component
    # is -s*D_i(S0); its overall minus sign is irrelevant to the norm.
    coordinate = delta*DK/s+sp.Rational(3, 2)*Q
    actual = s*(sp.diff(coordinate, s)*Ds+sp.diff(coordinate, DK)*DDK
                +sp.diff(coordinate, Q)*Qx*(-2*s*Ds))
    expected = delta*DDK+((s+1/s)*DK/h-3*s**2*Qx)*Ds
    return {"actual": actual, "expected": expected, "residual": sp.factor(actual-expected)}


@cache
def proof_checks():
    # For x in [-11/10,-9/10], p²>=9/40, h>=1 and |h'|/h<=3.
    amax = sp.Rational(5, 9)+sp.Rational(5, 6)
    forcing_slope = sp.Rational(3, 1)*3/(16*sp.Rational(9, 40))
    exponent = amax/10
    exponential_bound = 1/(1-exponent)
    q_bound = forcing_slope*exponential_bound/2
    qx_bound = forcing_slope+amax*q_bound/10
    return {"ODE_absolute_coefficient_upper": amax == sp.Rational(25, 18),
            "forcing_linear_in_distance_slope": forcing_slope == sp.Rational(5, 2),
            "Volterra_exponent_less_than_one": bool(0 < exponent < 1),
            "Volterra_exponential_majorant": exponential_bound == sp.Rational(36, 31),
            "Q_second_order_bound": q_bound == sp.Rational(45, 31),
            "Q_bound_below_three_halves": bool(q_bound < sp.Rational(3, 2)),
            "Qx_first_order_bound": qx_bound == sp.Rational(335, 124),
            "Qx_bound_below_three": bool(qx_bound < 3),
            "normal_clock_upper": bool(sp.Rational(21, 20)**2 > sp.Rational(11, 10)),
            "normal_source_quadratic_constant": bool(sp.Rational(3, 2)*sp.Rational(21, 20)*sp.Rational(3, 2) < sp.Rational(5, 2)),
            "normal_plus_inverse_below_three": bool(sp.Rational(21, 20)+sp.Rational(10, 9) < 3),
            "electric_Qx_constant_below_ten": bool(3*sp.Rational(11, 10)*3 < 10),
            "source_bound_is_not_a_remainder_or_inverse": True}


def bound_values(clock_distance_over_h, K_deviation, spatial_K_deviation, spatial_s_over_h):
    eta, dk, ddk, ds = [old_bounds.exact(value, name) for value, name in
                        ((clock_distance_over_h, "eta"), (K_deviation, "Delta_K"),
                         (spatial_K_deviation, "D_Delta_K"), (spatial_s_over_h, "D_s_over_h"))]
    if any(value.is_nonnegative is not True for value in (eta, dk, ddk, ds)) or eta > sp.Rational(1, 10):
        raise ValueError("Nonnegative declared source norms and eta<=1/10 required")
    return {"normal_source_upper": eta*dk+sp.Rational(5, 2)*eta**2,
            "electric_source_upper": eta*ddk+(3*dk+10*eta)*ds,
            "physical_curl_factor": "divide normalized electric bound by tau²",
            "full_action_or_retarded_remainder_claim": False}


@cache
def checks():
    ode, electric = ode_bridge(), electric_identity()
    bg = alignment.parent.old.background()
    return {"frozen_ODE_coefficient": ode["coefficient_residual"],
            "frozen_ODE_forcing": ode["forcing_residual"],
            "ODE_zero_forcing_at_clock": ode["forcing_zero_at_basepoint"],
            "h_derivative_bridge": sp.factor(sp.diff(bg["h"], modes.u)-3*bg["H"]*bg["h"]/2),
            "exact_spatial_electric_source": electric["residual"]}
