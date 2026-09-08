"""Exact domain bounds and strict input controls for the named candidate."""
from functools import cache

import sympy as sp

from . import cones, degeneracy, dynamics, geometry

ZETA_MAX = sp.Rational(1, 2000)


def positive_even(value):
    poly = sp.Poly(value, dynamics.u)
    return bool(poly.TC() > 0 and all(powers[0] % 2 == 0 and coefficient >= 0
                                    for powers, coefficient in poly.terms()))


@cache
def proof_checks():
    low, high = sp.Rational(47, 100), sp.Rational(53, 100)
    data = geometry.update()
    time, space = data["old_time"], data["old_space"]
    p = geometry.P
    principal = dynamics.principal()
    numerator, denominator = sp.fraction(dynamics.old.background()["J"])
    cone_numerator, cone_denominator = sp.fraction(cones.characteristic()["kappa"])
    ratio_bound = 3*sp.Rational(1, 100)/(8*sp.Rational(18, 19)*sp.Rational(9, 40))
    return {"rational_p_interval_contains_closed_tube": bool(low**2 < sp.Rational(9, 40)
                                                           and high**2 > sp.Rational(11, 40)),
            "trace_time_response_in_minus_six_minus_three": bool(time.subs(p, low) > -6 and time.subs(p, high) < -3),
            "trace_space_response_in_four_six": bool(space.subs(p, high) > 4 and space.subs(p, low) < 6),
            "time_response_strictly_increasing": sp.factor(sp.diff(time, p)-12*p-3/p**2) == 0,
            "space_response_strictly_decreasing": sp.factor(sp.diff(space, p)+1/p**2+sp.Rational(5, 4)/p**3) == 0,
            "positive_new_temporal_inverse_mass_lower": geometry.MU-sp.Rational(1, 3) == sp.Rational(8, 9),
            "positive_new_spatial_inverse_mass_lower": geometry.MU-sp.Rational(1, 4) == sp.Rational(35, 36),
            "both_new_mass_responses_above_18_over_19": geometry.MU-sp.Rational(1, 6) == sp.Rational(19, 18),
            "trace_map_infinity_norm": data["N"].norm(sp.oo) == 12,
            "transpose_trace_map_infinity_norm": data["N"].T.norm(sp.oo) == 2,
            "Woodbury_middle_inverse_upper": 1/(3-1/geometry.MU) == sp.Rational(11, 24),
            "uniform_quotient_inverse_below_6900": 25+25**2*12*2*sp.Rational(11, 24) == 6900,
            "trace_kinetic_relative_bound": ratio_bound == degeneracy.chart()["uniform_relative_upper"],
            "trace_kinetic_cannot_change_rank_on_tube": bool(ratio_bound < sp.Rational(1, 50)),
            "original_J_positive_numerator": positive_even(numerator),
            "original_J_positive_denominator": positive_even(denominator),
            "gradient_comparison_denominator_positive": positive_even(principal["bound_denominator"]),
            "gradient_comparison_strictly_below_1000": positive_even(principal["bound_witness"]),
            "chosen_curl_keeps_at_least_half_gradient_form": 1000*ZETA_MAX == sp.Rational(1, 2),
            "cone_kappa_positive_away_from_center_numerator": positive_even(sp.cancel(cone_numerator/dynamics.u**2)),
            "cone_kappa_positive_denominator": positive_even(cone_denominator),
            "point_transformation_regular_at_center": dynamics.source()["f"].subs(dynamics.u, 0) == 5,
            "first_order_uses_no_Theta_inverse": not sp.denom(sp.together(dynamics.first_order()["H"])).has(dynamics.old.theta),
            "primary_not_full_secondary_constraint_claim": degeneracy.chart()["nonlinear_secondary_constraint_rank_claim"] is False}


def exact(value, name):
    if isinstance(value, (bool, float, str)):
        raise TypeError(f"{name} must be an exact finite real")
    value = sp.sympify(value)
    if (not isinstance(value, sp.Expr) or value.has(sp.Float) or value.free_symbols
            or value.is_real is not True or value.is_finite is not True):
        raise ValueError(f"{name} must be an exact finite real constant")
    return value


def require_domain(p, coupling, time, momentum, *, punctured=True):
    p, zeta, point, q = [exact(value, name) for value, name in
                         ((p, "p"), (coupling, "zeta"), (time, "time"), (momentum, "q"))]
    if p.is_positive is not True or not bool(sp.Rational(9, 40) <= p**2 <= sp.Rational(11, 40)):
        raise ValueError("Require the original closed positive-p tube")
    if zeta.is_zero is True:
        return {"branch": "exact_auxiliary_control", "p": p, "zeta": zeta}
    if not bool(0 < zeta <= ZETA_MAX):
        raise ValueError("Require 0<zeta<=1/2000 for the certified principal branch")
    if q.is_positive is not True:
        raise ValueError("Require q>0")
    if punctured and point.is_zero is not False:
        raise ValueError("The v-velocity chart requires u!=0")
    return {"branch": "positive_principal_rolling_candidate", "p": p, "zeta": zeta,
            "u": point, "q": q, "nonlinear_or_UV_health_claim": False}


def units(mass_squared, time_scale, physical_curl):
    mass, tau, curl = [exact(value, name) for value, name in
                       ((mass_squared, "M²"), (time_scale, "tau"), (physical_curl, "zeta_physical"))]
    if mass.is_positive is not True or tau.is_positive is not True or curl.is_positive is not True:
        raise ValueError("Require positive M², tau and physical curl")
    return {"mu": geometry.MU, "normalized_zeta": curl/(mass*tau**2),
            "isolated_Proca_mass_squared_physical": mass/curl,
            "isolated_Proca_mass_squared_normalized": mass*tau**2/curl,
            "isolated_mass_is_not_a_coupled_gap_or_cutoff": True}
