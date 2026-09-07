"""Explicit Euclidean operator-norm bounds for all four coupled components.

Every norm is subordinate to the Euclidean vector norm, except a specified
source L1 norm. Entrywise absolute sums are only upper bounds on that norm.
The large uniform constant is a convergence bound, not a small-error pass
on the whole delta<=10^-9 domain. See the smaller exact example.
"""

from functools import cache

import sympy as sp

from . import operator
from .exact import DELTA_MAX, SLICE, nonnegative, parameters


def exponential_majorant(x):
    x = nonnegative(x, "x")
    if (sp.Rational(2, 3)-x).is_nonnegative is not True:
        raise ValueError("This elementary exponential cap requires 0<=x<=2/3")
    return 1+x+x**2/(2*(1-x/3))


def transfer_error(delta, kbar=1):
    delta, _ = parameters(delta, kbar)
    return 40000*delta**sp.Rational(1, 3)


def source_error(delta, source_l1, kbar=1):
    delta, _ = parameters(delta, kbar)
    return 3000000*delta*nonnegative(source_l1, "source_l1")


def prepared_error(delta, target_norm, source_l1, kbar=1):
    error = transfer_error(delta, kbar)
    return error*nonnegative(target_norm, "target_norm")+(2+error)*source_error(delta, source_l1, kbar)


def outer_picard_tail(radius, order):
    """Certified remainder of the interaction-picture outer Volterra series.

    The n-th kernel has L1 norm at most (64*r)^n/n!. Norms/inverses <=2
    on r<=SLICE; <=4 up to 2*SLICE. No point sampling is substituted.
    """
    radius = nonnegative(radius, "radius")
    if (2*SLICE-radius).is_nonnegative is not True:
        raise ValueError("The response/source outer chart is r<=2/100")
    if isinstance(order, bool) or not isinstance(order, int):
        raise TypeError("The Picard truncation order must be an integer")
    if order < 0:
        raise ValueError("The Picard truncation order must be nonnegative")
    factor = 2 if (SLICE-radius).is_nonnegative is True else 4
    return factor*(64*radius)**(order+1)/sp.factorial(order+1)


@cache
def coefficient_checks():
    """Polynomial pole isolation plus exact margins used by the continuous proof.

    Pointwise jet-profile interval bounds are independently replayed in
    independent.analytic_boxes. These elementary bounds also capture the
    vanishing E=O(delta+u^2), which a plain absolute interval bound would lose.
    """
    q = sp.Rational
    d = operator.derive()
    v, delta = sp.symbols("v delta", nonnegative=True)
    radius_den = (2+delta)*(1+v)**4-2
    base_den = delta+8*v
    numerator = 16*(1-v)*(8/((2+delta)*(1+v)**14)+1/(1+v)**2)
    expression = numerator/radius_den
    actual = d["mass"].subs({d["c"]: 2+delta, d["u"]: sp.sqrt(v)})
    residuals = {
        "exact_mass_numerator_denominator": sp.factor(actual-expression),
        "positive_denominator_difference": sp.expand(radius_den-base_den
            -delta*v*(4+6*v+4*v**2+v**3)-v**2*(12+8*v+2*v**2)),
        "mixed_denominator_square": sp.expand(base_den**2-32*delta*v-(delta-8*v)**2),
        "mass_outer_leading_coefficient": sp.limit(expression.subs(delta, 0)*v, v, 0)-10,
        "mass_remainder_center_path": sp.limit(expression.subs(v, 0)-80/delta, delta, 0)+32,
        "mass_remainder_punctured_path": sp.limit(expression.subs(delta, 0)-10/v, v, 0)+141,
    }
    margins = {
        "d12_below_8over7": q(8, 7)-q(101, 100)**12,
        "denominator_delta_v_coefficient": q(41, 10)-(4+6*q(1, 100)+4*q(1, 100)**2+q(1, 100)**3),
        "denominator_v2_coefficient": q(121, 10)-(12+8*q(1, 100)+2*q(1, 100)**2),
        "theta_derivative_below_7": 7-(6+q(12, 100)+q(6, 10)*q(128, 105)),
        "normalization_N_below_8": 8-(7+36*q(1, 100)),
        "omega_over_u_below_7": 7-q(48, 7),
        "omega_derivative_below_9": 9-(7+140*q(1, 100)),
        "gradient_square_below_5": 5-q(32, 7),
        "A_below_28": 28-(4*q(32, 7)+8),
        "C_below_28": 28-(8+18+84*q(1, 100)),
        "E_ratio_below_160": 160-(74+84),
        "analytic_B_correction_below_30": 30-(20+8+196*q(1, 100)),
        "total_bounded_B_below_200": 200-(152+30),
        "mass_delta_derivative_below_3_rminus4": 3-(4*q(1, 100)+q(5, 2)),
        "endpoint_fsum_squared_above_quarter": 8/(4*q(201, 100)*q(8, 7))-q(1, 4),
        "endpoint_fsum_squared_below_4": 4-(q(201, 100)*q(8, 7)+8)/8,
        "endpoint_frelative_squared_above_sixteenth": 2/(q(201, 100)*q(8, 7)+8)-q(1, 16),
        "endpoint_frelative_squared_below_1": 1-q(2, 7),
    }
    if any(sp.simplify(value) != 0 for value in residuals.values()):
        raise ValueError("An exact mass or chart identity failed")
    if any(sp.sympify(value).is_positive is not True for value in margins.values()):
        raise ValueError("A continuous coefficient bound lost its strict margin")
    return {"residuals": residuals, "strict_continuous_margins": margins}


@cache
def calibration():
    q = sp.Rational
    # All expressions in this table are independently replayed with Fraction.
    mass_remainder = 126+80*(q(41, 10)/32+q(121, 10)/64)
    outer_row = (29+(35+14*q(13, 4))*q(1, 3)
                 +q(200, 3)*q(1, 10)+q(160, 3)*q(1, 300)
                 +q(14, 3)*q(1, 30))
    outer_difference = (1+8*q(1, 10)**3
                        +(12+8*q(13, 4)+q(8, 3))*q(1, 3000)
                        +q(8, 3)*q(1, 30000)+q(8, 3)*q(1, 10000))
    inner_row = (29+(35+14*q(13, 4))*q(3, 8)+q(200, 3)*q(1, 8)
                 +q(1280, 3)*q(3, 512)+q(14, 3)*q(3, 64))
    telescoping = 128*2*4*2*2+q(1, 8)*4*2*2+1152*2*2+q(1, 8)*2*2+128*2
    whole = 4*9000+48+2
    source = 8+160000*16
    return {
        "delta_max": DELTA_MAX, "fixed_slice": SLICE,
        "fixed_momentum_min": sp.S.One, "fixed_momentum_max": sp.Integer(2),
        "t_max": q(1, 1000), "match_radius": "t=delta^(1/3)",
        "inner_eta_over_t": q(1, 8),
        "mass_remainder_upper": mass_remainder,
        "mass_remainder_margin_below_152": 152-mass_remainder,
        "bounded_B_minus_pole_upper": sp.Integer(200),
        "analytic_c_derivative_upper": sp.Integer(8),
        "outer_remainder_over_radius_upper": outer_row,
        "outer_remainder_margin_below_64": 64-outer_row,
        "outer_difference_over_delta_rminus2_upper": outer_difference,
        "outer_difference_margin_below_2": 2-outer_difference,
        "inner_remainder_over_radius_upper": inner_row,
        "inner_remainder_margin_below_72": 72-inner_row,
        "exp_two_thirds_upper": exponential_majorant(q(2, 3)),
        "exp_margin_below_2": 2-exponential_majorant(q(2, 3)),
        "inner_perturbation_exponent_max": q(576, 1000),
        "inner_exponent_margin_below_two_thirds": q(2, 3)-q(576, 1000),
        "middle_telescoping_coefficient": telescoping,
        "middle_margin_below_9000": 9000-telescoping,
        "full_telescoping_coefficient": sp.Integer(whole),
        "full_margin_below_40000": sp.Integer(40000-whole),
        "source_coefficient_before_enlargement": sp.Integer(source),
        "source_margin_below_3000000": sp.Integer(3000000-source),
        "nontrivial_delta": q(1, 10**21),
        "nontrivial_adapted_transfer_error_upper": transfer_error(q(1, 10**21)),
        "upper_box_error_not_small": transfer_error(DELTA_MAX),
        "physical_C_map_norm_upper": sp.Integer(32),
        "physical_C_inverse_norm_upper": sp.Integer(10),
        "fixed_slice_W_norm_upper": sp.Integer(40),
        "fixed_slice_W_inverse_norm_upper": sp.Integer(14),
        "adapted_endpoint_P_norm_upper": sp.Integer(3000),
        "adapted_endpoint_P_inverse_norm_upper": sp.Integer(280),
    }


def checks():
    d = calibration()
    for key, value in d.items():
        if "margin" in key and value.is_positive is not True:
            raise ValueError(f"A strict rational norm margin failed: {key}")
    if d["nontrivial_adapted_transfer_error_upper"] != sp.Rational(1, 250):
        raise ValueError("The stated small-error subrange failed")
    return d
