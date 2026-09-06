"""Exact rational initial-data screens, not a bound on a retarded Green operator."""

from fractions import Fraction

import sympy as sp

from . import background

ROOT_LOWER = Fraction(10511, 500)
ROOT_UPPER = Fraction(106, 5)


def fraction(value):
    value = sp.sympify(value)
    if value.has(sp.Float) or value.is_Rational is not True:
        raise TypeError("An exact finite rational coefficient is required")
    return Fraction(int(value.p), int(value.q))


def linear_root_enclosure(value):
    """Exact affine-Q(sqrt442) bound on a strict certified root interval."""
    value = sp.sympify(value)
    if value.has(sp.Float):
        raise TypeError("Floating inputs do not certify outward enclosures")
    root = sp.Symbol("root_442")
    polynomial = sp.Poly(sp.radsimp(value).subs(background.RADICAL, root), root)
    if polynomial.degree() > 1:
        raise ValueError("Only affine expressions in sqrt442 are certified")
    intercept, slope = fraction(polynomial.nth(0)), fraction(polynomial.nth(1))
    endpoints = (intercept+slope*ROOT_LOWER, intercept+slope*ROOT_UPPER)
    return {"constant": intercept, "root_coefficient": slope,
            "lower": min(endpoints), "upper": max(endpoints)}


def build():
    d = background.initial_jets(2)
    mass = linear_root_enclosure(d["mass_squared"])
    speed = linear_root_enclosure(d["c_light_squared"])
    ns = sp.Poly(d["N_sum"], background.A)
    rr = sp.Poly(d["Routh_frequency_squared"], background.A)
    if ns.degree() != 2 or rr.degree() != 2:
        raise ValueError("The specified asymmetric CD second-jet polynomials changed")
    ns_bounds = {str(j): linear_root_enclosure(ns.nth(j)) for j in range(3)}
    rr_bounds = {str(j): linear_root_enclosure(rr.nth(j)) for j in range(3)}
    q = Fraction
    # For 0<=A<=1/4, sqrt442<22 gives a positive lower bound on
    # 1037-(376+16sqrt442)A and an upper bound on its squared denominator.
    numerator_lower = q(1037)-(376+16*q(22))*q(1, 4)
    denominator_upper = q(16)*(60996+3047*q(22))
    omega_lower = numerator_lower**2/denominator_upper
    margins = {"strict_root_lower": q(442)-ROOT_LOWER**2,
               "strict_root_upper": ROOT_UPPER**2-q(442),
               "positive_mass_lower": mass["lower"],
               "locked_cone_positive": speed["lower"],
               "locked_cone_subluminal": 1-speed["upper"],
               "N_sum_constant_above_half": ns_bounds["0"]["lower"]-q(1, 2),
               "N_sum_linear_positive": ns_bounds["1"]["lower"],
               "N_sum_quadratic_positive": ns_bounds["2"]["lower"],
               "Routh_constant_below_three_halves": q(3, 2)-rr_bounds["0"]["upper"],
               "Routh_linear_negative": -rr_bounds["1"]["upper"],
               "Routh_quadratic_negative": -rr_bounds["2"]["upper"],
               "positive_omega_numerator_for_small_A": numerator_lower,
               "omega_squared_above_one_third_for_small_A": omega_lower-q(1, 3)}
    if mass["upper"] != q(9, 200) or any(value <= 0 for value in margins.values()):
        raise ValueError("A specified-data hierarchy enclosure failed")
    return {"root_interval": {"lower": ROOT_LOWER, "upper": ROOT_UPPER, "square": q(442)},
            "asymmetric_mass_squared_over_m_squared": mass,
            "conditional_locked_speed_squared": speed,
            "N_sum_coefficients_in_A": ns_bounds, "Routh_frequency_coefficients_in_A": rr_bounds,
            "small_A_omega_squared_lower": omega_lower,
            "small_A_domain": "0<=A<=1/4",
            "all_A_domain": "A>=0 at y0=2,rho0=1/2,CD h0=hsecond0=0",
            "strict_mass_led_ratio_lower": q(50, 9),
            "ratio_meaning": "max{Hdot_eff,omega^2}/m_alg^2 > 50/9; not an obstruction to every nonadiabatic reduction",
            "Routh_to_locked_normalization_bound": "N_sum>1/2,Omega_Routh^2<3/2, hence Omega_Routh^2<3N_sum in units m^2",
            "strict_positive_margins": margins}
