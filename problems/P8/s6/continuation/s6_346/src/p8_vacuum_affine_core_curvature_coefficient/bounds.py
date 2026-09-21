"""Exact finite-domain inequalities, not a numerical large-mass fit."""

from functools import cache
from math import factorial

import sympy as s

from . import core, source


@cache
def data():
    n = core.N
    L = s.Symbol("L", real=True)
    n0 = s.Integer(10) ** 6
    denominator = 315 * n**4 * (n - 2) ** 4 * (n - 1) ** 9
    numerator = s.cancel(-core.original_coefficient() * denominator).expand()
    poly = s.Poly(numerator.subs(s.log(n), L), n, L)
    tail = sum(
        abs(co) * 462 ** ex[1] / n0 ** (13 - ex[0])
        for ex, co in poly.terms()
        if ex != (13, 0)
    )
    d0 = (1 - 2 / n0) ** 4 * (1 - 1 / n0) ** 9
    exp_lower = sum(s.Rational(7, 3) ** j / factorial(j) for j in range(9))
    checks = {
        "entire_numerator_reconstruction": s.expand(
            numerator - poly.as_expr().subs(L, s.log(n))
        ),
        "exact_leading_numerator_coefficient": poly.coeff_monomial(n**13) - 69,
        "logarithmic_original_upper_budget": s.Integer(198) * s.Rational(7, 3) - 462,
        "original_contact_not_asymptotically_retuned": s.factor(
            source.CONTACT / source.CUBIC**2 - core.contact_ratio(source.HEAVY_MASS2)
        ),
    }
    return {
        "checks": checks,
        "gates": {
            "every_tail_degree_at_most12": all(
                ex[0] <= 12 for ex, co in poly.terms() if ex != (13, 0)
            ),
            "only_first_power_of_log": poly.degree(L) == 1,
            "all_tail_powers_nonnegative": all(ex[1] >= 0 for ex, co in poly.terms()),
            "strict_finite_positive_numerator": 69 - tail > 0,
            "finite_magnitude_lower_above_one_fifth": (69 - tail) / 315
            > s.Rational(1, 5),
            "finite_magnitude_upper_below_one_quarter": (69 + tail) / (315 * d0)
            < s.Rational(1, 4),
            "positive_Taylor_polynomial_proves_exp_7_over3_above10": exp_lower > 10,
            "original_mass_above_million": source.HEAVY_MASS2 > n0,
            "original_mass_below_10_power198": source.HEAVY_MASS2 < 10**198,
            "original_Born_normalized_bound_below_10_minus207": source.CUBIC**2
            / (2304 * source.HEAVY_MASS2)
            < s.Rational(1, 10) ** 207,
            "no_float_estimate_or_asymptotic_only_bound": True,
            "no_bound_on_independent_extra_parent_chi": True,
        },
        "whole_exact_numerator": numerator,
        "whole_exact_denominator": denominator,
        "whole_exact_tail_budget": tail,
        "whole_exact_denominator_lower": d0,
        "whole_positive_exponential_partial_sum": exp_lower,
        "whole_finite_domain_proof": "For n>=10^6 and0<=log(n)<=462, write numerator=69n^13+R. Every other term has degree<=12 and log power0 or1. The displayed rational tail sums every absolute coefficient times462^logpower/(10^6)^(13-degree), so |R|/n^13<=tail. The denominator divided by315n^17 is(1-2/n)^4(1-1/n)^9, between d0 and1. Exact rational gates give1/5<(69-tail)/315 and(69+tail)/(315d0)<1/4. Therefore -1/(4n^4)<c_core(n)<-1/(5n^4).",
        "whole_original_parameter_and_Born_bound": "Original10^6<n<10^198. The positive degree8 Taylor polynomial for exp(7/3) exceeds10 exactly, hence log10<7/3 and logn<462. Thus chi_core<0 and |chi_core|<g^4/(64pi^2 n^4)<g^4/(576n^4). With the existing original A0>4g^2/n^3, |chi_core|/A0<g^2/(2304n)<10^-207. This is only a local matching-coefficient bound; it is not an above-threshold Taylor approximation, a physical radiative remainder bound, or a bound on unmatched parent chi.",
    }
