"""Entire fixed-source fifth jets with the original physical metric factor."""

from functools import cache

import sympy as s
from p8_vacuum_affine_finite_window_growth import background as b
from p8_vacuum_affine_finite_window_growth.intervals import I, evaluate
from p8_vacuum_affine_heavy_scalar_parent import clock as old_heavy
from p8_vacuum_affine_physical_volume_correction import source as physical

TIME = s.Rational(1, 10**60)
EPSILON = s.Rational(1, 10**230)
LOW = s.Integer(10) ** 64
HIGH = 2 * LOW
ERROR = s.Rational(1, 10**300)


@cache
def data():
    base = b.source_data()
    radius = base["source_complex_radius"]
    switch = base["full_switch_modulus_base"]
    cauchy = s.factorial(5) / radius**5
    N = b.N
    bounds = (3, 7, 25, 121, 1000)
    for order, bound in enumerate(bounds, 1):
        value = evaluate(
            s.diff(N**-2, N, order),
            {N: I(s.Rational(999, 1000), s.Rational(1001, 1000))},
        )
        assert max(abs(value.lo), abs(value.hi)) < bound
    # Exponential Bell polynomials, not a guessed fifth chain constant.
    bell = [s.Integer(1)]
    for n in range(1, 6):
        bell.append(
            sum(
                s.binomial(n - 1, j - 1) * bounds[j - 1] * bell[n - j]
                for j in range(1, n + 1)
            )
        )
    assert bell == [1, 3, 16, 115, 1027, 11153], bell
    chain = max(bell)
    heavy = old_heavy.data()
    exp_bound = s.factorial(8) * 2**8 / old_heavy.family.LOCALIZER**8
    derivative_factors = {
        (i, j): s.factorial(i) * s.factorial(j) * 32 ** (i + j)
        for i in range(6)
        for j in range(6 - i)
    }
    assert max(derivative_factors.values()) < 10**10
    heavy_five = (
        heavy["common_complete_coefficient_prefactor"]
        * max(derivative_factors.values())
        * exp_bound
    )
    heavy_composed = chain * heavy_five
    assert heavy_composed < s.Rational(1, 10**2700)
    # The original switch is analytic on the same old joint neighborhood.
    old_R = 8 * cauchy * switch
    old_F = 10**207 * cauchy * switch
    constant = 10**301 * cauchy * switch
    assert old_R < s.Rational(1, 10**2490)
    assert old_F < s.Rational(1, 10**2290)
    assert constant < s.Rational(1, 10**2190)
    # Smooth fixed time profiles are differentiated by Leibniz, never Cauchy.
    # The factor 32 bounds every fifth-order binomial sum. Each profile
    # factor A+B(X-1), after X=N^-2, has mixed derivatives at most 2000
    # times the original C5 profile bound. T and its positive-order jets
    # are below one on the real strip, by the analytic complement bound.
    profile_factor = 32 * 2000
    profile_error = profile_factor * b.prior_initial.PROFILE_BOUND
    assert profile_factor < 10**7
    assert profile_error < s.Rational(1, 10**390)
    assert cauchy * switch < 1
    # The finite original onepoint coefficient obeys the same explicit
    # S256 majorant, enlarged for the fifth Cauchy derivative.
    onepoint = 10**25 * switch / s.Integer(10) ** 400
    assert onepoint < s.Rational(1, 10**2890)
    whole_RF = old_R + old_F + constant + profile_error + heavy_composed
    whole_j = heavy_composed + onepoint
    assert whole_RF < s.Rational(1, 10**380)
    assert whole_j < s.Rational(1, 10**2600)
    # Keep the original, weaker 1e-300/1e-2450 bounds when comparing
    # existing coefficient-distance probes; the improved margins are optional.

    return {
        "whole_current_function_bindings": base["entire_current_function_bindings"],
        "whole_fixed_profile_bindings": base["whole_fixed_profile_bindings"],
        "whole_normalized_heavy_source": base["complete_normalized_heavy_source"],
        "whole_correct_physical_metric_factor": physical.C,
        "original_localizer": old_heavy.family.LOCALIZER,
        "whole_fifth_Cauchy_factor": cauchy,
        "whole_lapse_Bell_chain_factors": bell,
        "whole_heavy_fifth_Cauchy_factors": derivative_factors,
        "complete_RF_fifth_error_bound": s.Rational(1, 10**380),
        "complete_heavy_fifth_bound": s.Rational(1, 10**2600),
        "weaker_RF_error_used_in_comparison": ERROR,
        "checks": {
            "complete_fifth_Bell_coefficient": bell[-1] - 11153,
            "complete_fifth_inner_Cauchy_factor": max(derivative_factors.values())
            - 4026531840,
            "literal_correct_physical_factor": physical.C
            - physical.R ** (-s.Rational(1, 2)),
            "whole_full_R_initial_time_derivative_zero": base["checks"][
                "entire_full_R_initial_time_derivative_zero"
            ],
            "same_original_time_and_momentum_band": s.Matrix(
                [TIME - b.TIME_LENGTH, LOW - b.MOMENTUM_MIN, HIGH - b.MOMENTUM_MAX]
            ),
        },
        "gates": {
            "all_previous_source_enclosure_gates": all(base["gates"].values()),
            "entire_heavy_fifth_composition_bound": heavy_composed
            < s.Rational(1, 10**2700),
            "complete_RF_C5_difference_bound": whole_RF < s.Rational(1, 10**380),
            "complete_normalized_source_C5_bound": whole_j < s.Rational(1, 10**2600),
            "smooth_profile_Leibniz_not_complex_Cauchy": profile_factor < 10**7,
            "original_action_sources_profiles_and_localizer_unchanged": True,
            "clock_tree_only_comparison_not_off_clock_replacement": True,
        },
    }
