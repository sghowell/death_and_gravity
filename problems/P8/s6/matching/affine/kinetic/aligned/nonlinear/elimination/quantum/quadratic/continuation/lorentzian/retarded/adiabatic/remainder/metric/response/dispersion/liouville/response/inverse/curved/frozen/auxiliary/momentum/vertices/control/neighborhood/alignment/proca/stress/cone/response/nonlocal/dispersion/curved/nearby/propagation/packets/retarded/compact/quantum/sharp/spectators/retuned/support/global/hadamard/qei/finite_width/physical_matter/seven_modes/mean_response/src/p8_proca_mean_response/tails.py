"""Global coefficient, weighted-response, parity and geometric bounds."""

from functools import cache

import sympy as sp
from p8_proca_seven_modes import model

from . import mean


def half_line_integral_upper(n):
    if isinstance(n, bool) or not isinstance(n, int) or n < 1:
        raise TypeError("Require a positive integer integral index")
    return sp.Rational(11, 7) * sp.binomial(2 * n - 2, n - 1) / 4 ** (n - 1)


@cache
def data():
    d = mean.data()
    u = mean.u
    t = 1 + u * u
    poly = sp.Poly(sp.factor(800 * t**18 * (d["J"] - 1 / (5 * t))), u)
    cn = (
        sp.Rational(3, 110)
        + sp.Rational(241, 200) * half_line_integral_upper(3)
        + 3 * half_line_integral_upper(9)
        + sp.Rational(3, 40) * half_line_integral_upper(20)
    )
    fn = (
        sp.Rational(25, 6)
        + 100 * half_line_integral_upper(1)
        + sp.Rational(25, 4) * half_line_integral_upper(12)
    )
    r0 = sp.Integer(1_500_004)
    C = sp.Integer(2_400_000_000)
    L = sp.Integer(550_000_000)
    actual_L = sp.Rational(9, 40) * C + sp.Rational(25, 4) * r0
    scalar = (L / 10 + 3 * C / 10) * half_line_integral_upper(6)
    T = sp.diag(sp.eye(3), -sp.eye(3))
    MP = model.data()["Proca_cartesian_generator"]
    parity = {
        mean.xi: mean.xi,
        mean.dp: -mean.dp,
        mean.rho: mean.rho,
        mean.pressure: mean.pressure,
    }
    checks = {
        "elementary_positive_integral_certifies_pi_below_22_over_7": sp.simplify(
            sp.integrate(u**4 * (1 - u) ** 4 / (1 + u * u), (u, 0, 1))
            - (sp.Rational(22, 7) - sp.pi)
        ),
        "elementary_global_diagonal_integral": sp.integrate(u / t**12, (u, 0, sp.oo))
        - sp.Rational(1, 22),
        "elementary_global_source_integral": sp.integrate(u / t**4, (u, 0, sp.oo))
        - sp.Rational(1, 6),
        "global_J_bound_remainder_is_actual_polynomial": sp.factor(
            d["J"] - 1 / (5 * t) - poly.as_expr() / (800 * t**18)
        ),
        "global_alpha_has_explicit_decaying_form": sp.factor(
            d["alpha"] - 3 * u * (4 * t**3 - 1) / t**4
        ),
        "global_Proca_reversal_is_anti_Hamiltonian": model.clean(T * MP * T + MP),
        "global_actual_Proca_generator_even_in_clock": model.clean(
            model.data()["actual_full_fourteen_phase_generator"][-6:, -6:].subs(u, -u)
            - model.data()["actual_full_fourteen_phase_generator"][-6:, -6:]
        ),
        "relative_mean_scale_flow_has_even_solution_parity": sp.factor(
            d["linearized_hat_scale_trace_flow"][0]
            .subs(u, -u)
            .subs(parity, simultaneous=True)
            + d["linearized_hat_scale_trace_flow"][0]
        ),
        "relative_mean_trace_flow_has_odd_solution_parity": sp.factor(
            d["linearized_hat_scale_trace_flow"][1]
            .subs(u, -u)
            .subs(parity, simultaneous=True)
            - d["linearized_hat_scale_trace_flow"][1]
        ),
        "relative_lapse_has_even_solution_parity": sp.factor(
            d["linearized_lapse"].subs(u, -u).subs(parity, simultaneous=True)
            - d["linearized_lapse"]
        ),
        "weighted_trace_damping_keeps_six_over_time": sp.factor(
            -3 * d["H"] + sp.diff(t**3, u) / t**3 + 6 * u / t
        ),
    }
    for n in range(2, 21):
        checks["elementary_half_line_integral_recurrence_" + str(n)] = sp.factor(
            half_line_integral_upper(n)
            - sp.Rational(2 * n - 3, 2 * n - 2) * half_line_integral_upper(n - 1)
        )
    return {
        "global_J_lower": 1 / (5 * t),
        "positive_J_remainder_polynomial_coefficients": poly.all_coeffs(),
        "global_positive_relative_vector_density_bound": "rho_P(u)<=1500004*eta/(1+u^2)^4",
        "weighted_phase": "Y=(xi,100*(1+u^2)^3*delta_p)",
        "weighted_generator_integrable_entry_sum_upper": cn,
        "weighted_forcing_L1_upper_per_initial_density": fn,
        "weighted_global_Gronwall_factor_upper": sp.Integer(9),
        "global_weighted_phase_upper_per_eta": C,
        "global_lapse_upper_per_eta": L,
        "derived_global_lapse_bound_before_rounding": actual_L,
        "global_linear_physical_log_scale_upper_per_eta": sp.Integer(2_700_000_000),
        "global_actual_frame_log_scale_upper_per_eta": sp.Integer(3_000_000_000),
        "global_matter_field_upper_per_eta": sp.Integer(400_000_000),
        "global_scalar_density_upper_per_eta": sp.Integer(81_000_000),
        "global_relative_scalar_density_upper_per_eta": sp.Integer(16_200_000_000),
        "global_lapse_decay_bound_per_eta": sp.Rational(3, 10) * C * sp.Abs(u) / t**3
        + sp.Rational(3, 40) * C / t**11
        + sp.Rational(25, 4) * r0 / t**3,
        "positive_tail_scale_limit_error_bound_per_eta_for_u_ge_one": sp.Rational(
            3, 220
        )
        * C
        / t**11
        + sp.Rational(241, 1000) * C / u**5
        + sp.Rational(25, 6) * r0 / t**3,
        "complete_first_order_frame_representative": "N_eta=1+n_eta; hat_a_eta=a*exp(xi_eta); a_eta=hat_a_eta*((h-1+N_eta^-2)/h)^(-1/4)",
        "representative_is_not_an_exact_semiclassical_solution": True,
        "checks": checks,
        "bounds": {
            "strict_positive_polynomial_remainder_all_global_times": all(
                c >= 0 for c in poly.all_coeffs()
            )
            and poly.nth(0) > 0,
            "weighted_generator_integral_below_two": cn < 2,
            "weighted_source_integral_below_170": fn < 170,
            "global_weighted_phase_bound": 9 * 170 * r0 < C,
            "global_lapse_bound": actual_L < L,
            "global_linear_physical_scale_bound": C + L / 2 < 2_700_000_000,
            "global_actual_frame_physical_scale_bound": C + L < 3_000_000_000,
            "global_matter_field_bound": scalar < 400_000_000,
            "global_scalar_relative_density_bound": 6 * 2_700_000_000 == 16_200_000_000,
            "calibrated_lapse_is_globally_below_one_e_minus_ten": L
            * sp.Rational(1, 10**20)
            < sp.Rational(1, 10**10),
            "calibrated_actual_frame_scale_is_globally_below_one_e_minus_ten": sp.Rational(
                3_000_000_000, 10**20
            )
            < sp.Rational(1, 10**10),
        },
    }
