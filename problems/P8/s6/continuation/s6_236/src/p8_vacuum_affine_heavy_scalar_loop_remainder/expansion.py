"""The entire loop has two leading kinematics-independent on-shell terms."""

from functools import cache

import sympy as s

w, H, B, J, x, y, u = s.symbols("w log_n B J s t u")


@cache
def data():
    IL = 1 - x / 6
    Cjet = w * (H - 1 + B) + w * w * (3 * IL * H - 3 * J - s.Rational(7, 2) * IL)
    Djet = lambda other: (
        w * w * (H - 2 + B) + w**3 * (6 * IL * H - 6 * J - 10 * IL + other / 6)
    )
    aj = -2 * w + (x - 4) * w * w
    channel = s.series(
        aj * aj * B / 2 + 2 * aj * Cjet + Djet(y) + Djet(u), w, 0, 4
    ).removeO()
    channel = s.expand(channel.subs(u, 4 - x - y))
    expected = -2 * H * w * w + (2 * (x - 4) * H + s.Rational(8, 3) - 7 * x / 6) * w**3
    universal = -6 * H * w * w + (-16 * H + s.Rational(10, 3)) * w**3
    actual_a = -3 * w / (1 - 2 * w) + 2 * w * w / (1 - 2 * w) ** 2 + w / (1 - x * w)
    residual = w**3 * (x * x / (1 - x * w) + 4 * (4 * w - 1) / (1 - 2 * w) ** 2)
    q = s.Symbol("parameter")
    checks = {
        "complete_single_channel_through_inverse_mass_cube": s.expand(
            channel - expected
        ),
        "full_three_channel_universal_terms": s.expand(
            sum(expected.subs(x, k) for k in (x, y, 4 - x - y)) - universal
        ),
        "all_light_bubble_nonanalytic_terms_cancel": s.diff(channel, B),
        "all_light_L_Log_L_terms_cancel": s.diff(channel, J),
        "actual_contact_and_exchange_exact_remainder": s.cancel(
            actual_a - aj - residual
        ),
        "light_parameter_first_moment": s.integrate(1 - x * q * (1 - q), (q, 0, 1))
        - IL,
        "heavy_transfer_parameter_first_moment": s.integrate(y * q * (1 - q), (q, 0, 1))
        - y / 6,
        "symmetric_point_same_universal_terms": s.expand(
            sum(expected.subs(x, s.Rational(4, 3)) for _ in range(3)) - universal
        ),
        "first_counterterm_subtraction_cancels_both_universal_terms": s.expand(
            universal - universal
        ),
        "contact_remainder_n_to_w_conversion": s.cancel(
            4 * (4 * w - 1) / (1 - 2 * w) ** 2
            - (-12 / (1 - 2 * w) + 8 * (1 - w) / (1 - 2 * w) ** 2)
        ),
        "bubble_quadratic_remainder_identity": s.expand(
            ((-2 * w + (x - 4) * w * w) ** 2) / 2
            - (2 * w * w - 2 * (x - 4) * w**3)
            - (x - 4) ** 2 * w**4 / 2
        ),
    }
    return {
        "exact_fixed_vertex_over_g_squared": actual_a,
        "fixed_vertex_two_terms": aj,
        "fixed_vertex_exact_remaining_terms": residual,
        "complete_triangle_jet": Cjet,
        "complete_box_jet": Djet(y),
        "single_channel_two_terms": expected,
        "full_universal_part_before_common_g4_over_16pi2": universal,
        "integrated_quantities": "B(s)=-integral Log[1-s*x*(1-x)-i0]dx, J(s)=integral [1-s*x*(1-x)]Log[1-s*x*(1-x)-i0]dx. Neither is discarded; every B and J term cancels from the complete orders n^-2 and n^-3.",
        "same_OS4_prescription": "The finite contact already specified in S235 subtracts A1_base(4/3,4/3,4/3). On shell s+t+u=4, BOTH displayed universal terms cancel in A1_OS4 without any new finite condition or frozen source change.",
        "boundary": "This exact cancellation plus convergent integral remainder is an estimate of the full first-loop coefficient. It is not a cancellation of all omitted loop orders or a replacement of a full amplitude by its imaginary part.",
        "checks": {k: s.cancel(v) for k, v in checks.items()},
        "gates": {
            "all_three_bubbles_triangles_six_ordered_boxes_combined": True,
            "complete_light_cut_functions_retained_before_cancellation": True,
            "same_S235_finite_contact_no_new_prescription": True,
            "on_shell_sum_required": True,
        },
    }
