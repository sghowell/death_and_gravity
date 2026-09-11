"""Explicit small coefficient error budget and unchanged algebraic constraint pivot."""

from functools import cache

import sympy as s
from p8_vacuum_affine_coupled_response import inverse as old_inverse
from p8_vacuum_analytic_affine_parent import dynamics

from . import profile

EPS = profile.EPS


@cache
def data():
    radius = s.Rational(1, 32)
    ratio = s.Rational(2, 5)
    exact = {k: 2 * s.factorial(k) * 32**k * ratio**profile.N for k in range(6)}
    switch_upper = s.Rational(1, 10**390)
    mixed = {
        (i, k): (2 * EPS if k <= 1 else 10 * EPS * switch_upper)
        for i in range(6)
        for k in range(6)
    }
    global_radius = s.Rational(1, 64 * profile.N)
    global_T = {k: 2 * s.factorial(k) * (64 * profile.N) ** k for k in range(6)}
    global_mixed = {
        (i, k): ((5 * EPS) if k == 0 else (4 * global_T[k] + k * global_T[k - 1]) * EPS)
        for i in range(6)
        for k in range(6)
    }
    d = dynamics.clock_lapse()
    lower = s.Rational(1215, 800) * s.Rational(4, 5) ** 18
    delta = s.symbols("delta", positive=True)
    A, B = s.symbols("A B", real=True)
    lapse = s.symbols("N", positive=True)
    # The full chart preserves X=N^-2 and has volume R^-3/4.
    R = 1 + 2 * delta * (lapse**-2 - 1)
    density = lapse * R ** -s.Rational(3, 4) * (A + B * (lapse**-2 - 1))
    shift = (21 * delta**2 - 3 * delta) * A / 2 + (1 - 6 * delta) * B
    C1 = old_inverse.first_row()["physical_lapse_reconstruction_C1_upper"] + 19 * EPS
    checks = {
        "complex_disc_denominator_floor": s.Rational(3, 4)
        - radius
        - s.Rational(23, 32),
        "complex_disc_numerator_ceiling": s.Rational(1, 4) + radius - s.Rational(9, 32),
        "complex_disc_ratio_bound": (s.Rational(9, 32) / s.Rational(23, 32))
        - s.Rational(9, 23),
        "Cauchy_radius_power": radius**-5 - 32**5,
        "positive_lapse_denominator_at_slab_endpoint": d["lapse_denominator"].subs(
            d["u"], s.Rational(1, 2)
        )
        - 800 * s.Rational(5, 4) ** 18,
        "positive_lapse_numerator_constant": d[
            "positive_lapse_numerator_in_u_squared"
        ].subs(next(iter(d["positive_lapse_numerator_in_u_squared"].free_symbols)), 0)
        - 1215,
        "literal_new_lapse_pivot_shift": s.simplify(
            s.diff(density, lapse, 2).subs(lapse, 1) / 2 - shift
        ),
        "bounce_lapse_pivot_shift": s.expand(
            shift.subs(
                {
                    delta: s.Rational(1, 2),
                    A: -profile.P,
                    B: -(profile.rho + profile.P) / 2,
                }
            )
            - (profile.rho - s.Rational(7, 8) * profile.P)
        ),
        "lapse_A_coefficient_at_upper_delta": ((21 * delta**2 - 3 * delta) / 2).subs(
            delta, s.Rational(1, 2)
        )
        - s.Rational(15, 8),
        "lapse_A_coefficient_increasing": s.diff((21 * delta**2 - 3 * delta) / 2, delta)
        - (21 * delta - s.Rational(3, 2)),
        "lapse_B_coefficient_at_upper_delta": (6 * delta - 1).subs(
            delta, s.Rational(1, 2)
        )
        - 2,
        "new_local_recovery_quantum_sum": s.Integer(7)
        + s.Rational(3, 2)
        + 12
        + 9
        - s.Rational(59, 2),
        "new_local_recovery_minus_previous_quantum": s.Rational(59, 2)
        - s.Rational(21, 2)
        - 19,
        "global_middle_power_phase_upper": 16 * profile.N * global_radius
        - s.Rational(1, 4),
        "global_Cauchy_radius": 1 / global_radius - 65536,
    }
    return {
        "quantitative_domain": "|u|<=1/2 and3/4<=X<6/5; all36 mixed derivatives with0<=u_order<=5 and0<=X_order<=5",
        "reference_stress_derivative_upper": EPS,
        "complex_X_disc_radius": radius,
        "complex_ratio_upper": ratio,
        "Cauchy_T_minus_one_derivative_envelopes": exact,
        "uniform_T_minus_one_derivative_upper_through_five": switch_upper,
        "mixed_correction_derivative_upper": mixed,
        "higher_X_derivative_upper": 10 * EPS * switch_upper,
        "common_mixed_derivative_strict_upper": 2 * EPS,
        "full_X_strip_domain": "|u|<=1/2 and-1/4096<X<6/5; all36 mixed derivatives through five in each variable",
        "full_X_strip_complex_disc_radius": global_radius,
        "full_X_strip_T_derivative_envelopes": global_T,
        "full_X_strip_mixed_correction_derivative_upper": global_mixed,
        "full_X_strip_common_mixed_derivative_strict_upper": s.Rational(1, 10**742),
        "classical_original_J_lower_on_slab": lower,
        "new_classical_J_shift": shift,
        "new_classical_J_shift_absolute_upper": 4 * EPS,
        "new_classical_J_lower_on_slab": lower - 4 * EPS,
        "new_retained_first_row_C1_upper": C1,
        "new_first_row_quantum_added_upper": s.Rational(59, 2) * EPS,
        "checks": checks,
        "gates": {
            "disc_ratio_strictly_below_two_fifths": s.Rational(9, 23) < ratio,
            "complex_q_power_below_half": ratio**profile.N < s.Rational(1, 2),
            "all_Cauchy_envelopes_below_one_e_minus_390": all(
                v < switch_upper for v in exact.values()
            ),
            "normalized_vacuum_constant_below_one_e_minus_789": profile.PV
            < s.Rational(1, 10**789),
            "common_coefficient_budget_valid": s.Rational(5, 4) * EPS + 2 * profile.PV
            < 2 * EPS,
            "first_X_budget_valid": EPS + 2 * EPS * switch_upper < 2 * EPS,
            "all_higher_X_jets_fit_common_budget": 10 * EPS * switch_upper < 2 * EPS,
            "lapse_shift_coefficient_sum_below_four": s.Rational(15, 8) + 2 < 4,
            "actual_positive_lapse_margin_preserved": lower - 4 * EPS
            > s.Rational(1, 100),
            "new_local_lapse_recovery_still_below_15000": C1 < 15000,
            "global_radius_below_one_eighth": global_radius < s.Rational(1, 8),
            "global_outer_disc_ratio_below_half": (s.Rational(1, 4) + global_radius)
            / (s.Rational(3, 4) - global_radius)
            < s.Rational(1, 2),
            "global_middle_power_phase_below_pi_over_two": s.Rational(1, 4) < s.pi / 2,
            "global_outer_power_below_half": s.Rational(1, 2) ** profile.N
            < s.Rational(1, 2),
            "all_full_strip_mixed_bounds_below_one_e_minus_742": all(
                v < s.Rational(1, 10**742) for v in global_mixed.values()
            ),
        },
    }
