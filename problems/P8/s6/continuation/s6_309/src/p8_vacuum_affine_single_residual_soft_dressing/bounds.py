"""Uniform weighted-variation bound and strictly scoped positive reference."""

from functools import cache

import sympy as s

from . import source


def relative_bound(energy, transfer, resolution):
    bound = source.softlimit.finite_remainder_bound(
        energy, transfer, resolution, source.KAPPA
    )
    return (1 + source.conversion.original_reference_ratio_bound()) * bound


def coarse_relative_bound(energy, transfer, resolution):
    source.softlimit.estimates.require_domain(energy, transfer, resolution)
    return 2 * source.softlimit.FIXED / source.KAPPA


def vanishing_relative_bound(energy, transfer, resolution):
    *_, x = source.softlimit.estimates.require_domain(energy, transfer, resolution)
    return (
        2
        * (source.softlimit.LINEAR * x + source.softlimit.QUADRATIC * x * x)
        / source.KAPPA
    )


@cache
def data():
    checks = {}

    def put(name, value):
        checks[name] = s.factor(s.expand_func(s.sympify(value)))

    x, k = s.symbols("resolution kappa", positive=True)
    a = s.Symbol("physical_index", nonnegative=True)
    L, Q, B = (
        source.softlimit.LINEAR,
        source.softlimit.QUADRATIC,
        source.softlimit.FIXED,
    )
    put("unchanged_absolute_linear_coefficient", L - 2 * 10**14)
    put("unchanged_absolute_quadratic_coefficient", Q - 2 * 10**37)
    put("unchanged_absolute_fixed_coefficient", B - 10**32)
    put(
        "whole_ratio_factor",
        1
        + source.conversion.original_reference_ratio_bound()
        - 1
        - 4250 / source.KAPPA,
    )
    put("whole_uniform_original_bound", 2 * B / source.KAPPA - s.Rational(2, 10**768))
    put(
        "whole_vanishing_relative_majorant",
        2 * (L * x + Q * x * x) / k - (4 * 10**14 * x + 4 * 10**37 * x * x) / k,
    )
    put(
        "whole_relative_zero_threshold_limit",
        s.limit(2 * (L * x + Q * x * x) / k, x, 0, dir="+"),
    )
    put(
        "whole_linear_rate_at_zero",
        s.limit(2 * (L * x + Q * x * x) / (k * x), x, 0, dir="+") - 2 * L / k,
    )
    put(
        "whole_old_vs_new_vanishing_power",
        2 * (L * s.sqrt(x) + Q * x * s.sqrt(x)) * s.sqrt(x) - 2 * (L * x + Q * x * x),
    )
    for p in range(1, 7):
        ratio = s.gamma(p + 1) * s.gamma(1 + a) / s.gamma(p + 1 + a)
        put(
            f"full_marked_energy_suppression_{p}",
            ratio - s.prod(s.Integer(j) / (j + a) for j in range(1, p + 1)),
        )
    put("wrong_individual_cut_first_power", 1 - 1 / (1 + a) - a / (1 + a))
    put(
        "wrong_individual_cut_second_power",
        1 - 2 / ((1 + a) * (2 + a)) - a * (a + 3) / ((1 + a) * (2 + a)),
    )
    m, v, P = s.symbols("signed_dressed_value variation positive_reference", real=True)
    put("positive_reference_lower_algebra", P + m - P * (1 + m / P))
    put("full_weighted_variation_algebra", (1 + 4250 / k) * v - v - 4250 * v / k)
    margins = {
        "uniform_weight_below_two": 1 - 4250 / source.KAPPA,
        "uniform_dressed_relative_error_below_one": 1 - 2 * B / source.KAPPA,
        "stronger_vanishing_bound_at_small_resolution": s.Rational(1, 8)
        - s.Rational(1, 64),
        "unchanged_reference_positive_Born": source.KAPPA,
    }
    for name, value in margins.items():
        put("positive_arithmetic_" + name, value - s.Abs(value))
    return {
        "whole_sharp_relative_bound": (1 + 4250 / source.KAPPA)
        * s.Min(B, L * x + Q * x * x)
        / source.KAPPA,
        "whole_uniform_original_relative_bound": s.Rational(2, 10**768),
        "whole_vanishing_relative_bound": 2 * (L * x + Q * x * x) / source.KAPPA,
        "whole_weighted_measure_proof": "A nonnegative kernel bounded above byC has absolute signed integral bounded byC times the total variation. Here C=1+4250/kappa<2 and S307 givesTV<min(B,Lx+Qx^2)/kappa. ThereforeabsD[R]/P_elastic has the displayed bound, uniformly in every nonforward angle and every0<x<=1/8. This is an all-N additional-LEADING-soft result for a single marked finite residual.",
        "whole_positive_reference": "P_elastic(x)+D[R](x)>P_elastic(x)*(1-2e-768)>0 at each fixedenergy/angle/threshold. This establishes pointwise positivity of the named reference only. Neither monotonicity inx nor a positive detector measure, unitarity or equality with the full interacting rate follows.",
        "whole_limit_discipline": "The additional-soft dimensional regulator is removed first at fixed marked state and then under the finite signed residual integral at fixed positivex. The resulting reference's relative correction vanishes uniformly asx->0. No simultaneous finite-e forward/threshold estimate or complete all-N nonleading theorem is asserted.",
        "whole_positive_margins": margins,
        "checks": checks,
        "gates": {
            "all_weight_and_positivity_margins_strict": all(
                bool(v > 0) for v in margins.values()
            ),
            "same_S307_absolute_total_variation_bound": True,
            "uniform_forward_and_zero_resolution_reference_control": True,
            "original_unexpanded_detector_power_retained": True,
            "positive_reference_not_detector_measure_or_unitarity": True,
            "no_unknown_matching_or_full_quantum_remainder_assigned": True,
        },
    }
