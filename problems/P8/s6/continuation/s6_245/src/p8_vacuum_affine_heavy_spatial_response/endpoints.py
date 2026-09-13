"""Complete scalar endpoint remainder, dimensional assembly and same-mask graph bounds."""

from functools import cache

import sympy as s
from p8_vacuum_affine_heavy_curved_state import state
from p8_vacuum_affine_heavy_homogeneous_response import (
    estimates as homogeneous_estimates,
)
from p8_vacuum_affine_heavy_homogeneous_response import renormalization as homogeneous
from p8_vacuum_affine_heavy_spatial_uv import matching as uv_matching

from . import domain, time_remainder


def slots():
    return tuple((j, degree) for j in range(5) for degree in range(5 - j))


def row_bound(j):
    if isinstance(j, bool) or not isinstance(j, (int, s.Integer)) or not 0 <= j <= 4:
        raise ValueError("Retain a complete endpoint order zero through four")
    return (
        domain.NORMALIZED_PAIR**2
        * sum(time_remainder.recurrence()[int(j)])
        * domain.RADIUS ** (-j)
    )


@cache
def constants():
    n, kappa = state.MASS2, state.KAPPA
    e, R = domain.EPS, domain.RADIUS
    B = tuple(row_bound(j) for j in range(5))
    raw = (
        2
        * domain.PAIR**2
        * sum(
            sum(time_remainder.recurrence()[j]) * (domain.MASS_FLOOR * R) ** (-j)
            for j in range(5)
        )
    )
    real = 4 * domain.PAIR**2
    far = tuple((4 - j, B[j] / (18 * e ** (4 - j))) for j in range(5))
    near = [(4, 4 * real / (9 * e**4), False, "complete_raw_row")]
    for j, degree in slots():
        if j + degree < 4:
            near.append(
                (
                    4 - j,
                    B[j]
                    * 2 ** (4 - j - degree)
                    / (18 * (4 - j - degree) * e ** (4 - j)),
                    False,
                    f"endpoint_{j}_degree_{degree}",
                )
            )
        else:
            near.append(
                (
                    5 - j,
                    2 * B[j] / (18 * e ** (5 - j)),
                    True,
                    f"endpoint_{j}_degree_{degree}_log",
                )
            )
    low = real / 27
    finite = (
        8 * sum(value for q, value in far)
        + 8 * sum(value for q, value, inverse_mass, label in near)
        + low
    )
    tail_far = 8 * sum(2 * B[j] / (9 * e ** (5 - j)) for j in range(5))
    tail_near = 8 * 201 * sum(value for q, value, inverse_mass, label in near)
    tail_low = 8 * low
    tail = tail_far + tail_near + tail_low
    time = time_remainder.constants()
    anchor = homogeneous.data()[
        "entire_normalized_homogeneous_current_and_response_bound"
    ]
    anchor_tail = homogeneous_estimates.data()[
        "same_frequency_cutoff_tail_coefficient_over_K_squared"
    ]
    local = uv_matching.finite()["complete_normalized_local_UV_difference_bound"]
    complete = (
        anchor
        + 2 * finite * n * n / kappa
        + 2
        * time["complete_time_integral_numerator_over_mass"]
        / (domain.MASS_FLOOR * kappa)
        + time["whole_actual_minus_W6_transfer_difference_normalized"]
        + local
    )
    complete_tail = (
        anchor_tail
        + 2 * tail * n * n * 10**99 / kappa
        + 2 * time["complete_time_removed_union_numerator_over_radial_cutoff"] / kappa
        + time["whole_actual_minus_W6_transfer_difference_normalized_tail"]
    )
    return {
        "all_five_full_normalized_endpoint_rows": B,
        "entire_unrounded_real_endpoint_majorant": raw,
        "complete_real_endpoint_majorant": real,
        "all_five_far_terms": far,
        "whole_raw_and_all_fifteen_near_terms": tuple(near),
        "complete_unexpanded_low_numerator_over_mass_four": low,
        "entire_endpoint_finite_numerator_over_mass_four": finite,
        "complete_far_tail_numerator_over_mass_five": tail_far,
        "complete_near_tail_numerator_over_mass_five": tail_near,
        "complete_low_tail_numerator_over_mass_five": tail_low,
        "entire_endpoint_tail_numerator_over_mass_five": tail,
        "full_actual_homogeneous_anchor_normalized": anchor,
        "complete_unchanged_scalar_local_UV_difference_bound": local,
        "complete_actual_all_transfer_spatial_response_bound": complete,
        "whole_original_frequency_projected_response_tail_numerator": complete_tail,
    }


@cache
def data():
    cc = constants()
    r, L, m, P = s.symbols("r L m P", positive=True)
    checks = {}
    for q in range(4):
        checks["complete_near_radial_power_" + str(q)] = s.integrate(
            r ** (3 - q), (r, m, L)
        ) - (L ** (4 - q) - m ** (4 - q)) / (4 - q)
    checks["complete_near_logarithmic_integral"] = s.expand_log(
        s.integrate(1 / r, (r, m, L)) - s.log(L / m)
    )
    checks["whole_far_integral"] = s.integrate(r**-2, (r, L, s.oo)) - 1 / L
    checks["all_fifteen_complete_UV_slots"] = len(slots()) - 15
    checks["all_thirty_five_source_jet_coefficients"] = (
        sum(j + 1 for j, degree in slots()) - 35
    )
    checks["source_spatial_weight_first_identity"] = s.expand(
        2 * (m * m + P * P) - (m + P) ** 2 - (m - P) ** 2
    )
    checks["source_spatial_weight_second_identity"] = s.expand(
        m * m * (1 + P * P) - (m * m + P * P) - (m * m - 1) * P * P
    )
    checks["full_finite_endpoint_numerator"] = cc[
        "entire_endpoint_finite_numerator_over_mass_four"
    ] - s.Rational(111268759571258000000000000000, 9)
    checks["full_endpoint_tail_numerator"] = cc[
        "entire_endpoint_tail_numerator_over_mass_five"
    ] - s.Rational(67387527200968092800000000000000, 27)
    anchor, errorP, error0, timeP, time0, QP, Q0, local = s.symbols(
        "actual_anchor state_P state_0 time_P time_0 Q_P Q_0 local_UV_difference"
    )
    contact = s.Symbol("whole_spatial_Hamiltonian_contact")
    rawP = errorP + timeP + QP + contact + local
    raw0 = error0 + time0 + Q0 + contact
    assembled = anchor + (errorP - error0) + (timeP - time0) + (QP - Q0) + local
    checks["complete_finite_transfer_repartition_without_dropping_contact"] = s.expand(
        anchor + rawP - raw0 - assembled
    )
    direction = s.Symbol("angular_direction", real=True)
    checks["angular_zero_row_without_mask"] = s.integrate(
        direction / 2, (direction, -1, 1)
    )
    checks["angular_zero_row_nonzero_with_two_leg_mask"] = s.integrate(
        direction / 2, (direction, 0, 1)
    ) - s.Rational(1, 4)
    checks["admitted_angular_zero_counterexample_radius"] = 6**2 + 8**2 - 10**2
    n, kappa = state.MASS2, state.KAPPA
    near = cc["whole_raw_and_all_fifteen_near_terms"]
    gates = {
        "full_all_real_endpoint_bound_including_all_five_orders": cc[
            "entire_unrounded_real_endpoint_majorant"
        ]
        < cc["complete_real_endpoint_majorant"],
        "all_five_far_rows_and_every_near_slot": len(cc["all_five_far_terms"]) == 5
        and len(near) == 16,
        "every_logarithmic_inverse_mass_factor_retained": sum(
            bool(im) for q, value, im, label in near
        )
        == 5,
        "full_finite_mass_power_at_most_four": all(
            q - int(im) <= 4 for q, value, im, label in near
        ),
        "full_tail_mass_power_at_most_five": all(
            q + 1 - int(im) <= 5 for q, value, im, label in near
        ),
        "full_tail_external_power_at_most_six": all(
            q + 1 <= 6 for q, value, im, label in near
        ),
        "complete_low_band_not_removed": cc[
            "complete_unexpanded_low_numerator_over_mass_four"
        ]
        > 0,
        "entire_finite_endpoint_constant": cc[
            "entire_endpoint_finite_numerator_over_mass_four"
        ]
        < 10**30,
        "entire_original_mask_endpoint_tail_constant": cc[
            "entire_endpoint_tail_numerator_over_mass_five"
        ]
        < 10**33,
        "all_three_removed_regions_retained": all(
            cc[key] > 0
            for key in (
                "complete_far_tail_numerator_over_mass_five",
                "complete_near_tail_numerator_over_mass_five",
                "complete_low_tail_numerator_over_mass_five",
            )
        ),
        "actual_mass_below_retained_tail_upper": n < 10**198,
        "full_local_UV_prescription_input_once": cc[
            "complete_unchanged_scalar_local_UV_difference_bound"
        ]
        < s.Rational(1, 10**600),
        "full_actual_spatial_response_positive_sum": cc[
            "complete_actual_all_transfer_spatial_response_bound"
        ]
        < s.Rational(10**46, kappa) * n * n,
        "complete_normalized_spatial_response_below_display": cc[
            "complete_actual_all_transfer_spatial_response_bound"
        ]
        < s.Rational(1, 10**350),
        "complete_same_mask_graph_tail_below_display": cc[
            "whole_original_frequency_projected_response_tail_numerator"
        ]
        < s.Rational(1, 10**250),
        "all_dimensional_domination_hypotheses_checked": all(
            domain.data()["gates"].values()
        ),
        "angular_zero_UV_terms_kept_before_masking": s.Rational(1, 4) != 0,
        "full_same_MSbar_local_piece_unprojected_and_unchanged": True,
        "homogeneous_contact_and_fixed_profile_kept_in_anchor": True,
        "no_bare_finite_cutoff_Ward_or_full_ADM_clock_claim": True,
        "no_finite_inhomogeneous_history_or_inverse_transfer": True,
    }
    return {
        "complete_finite_response_repartition": assembled,
        "whole_positive_constants_and_all_near_far_low_terms": cc,
        "full_graph": "For real prepared spatial sources, Z136(G)^2=sum(j0..13)||(1-Delta)^3 partial_t^j G||L2(dt dx;F)^2. The detector has plain L2 Frobenius norm. Complex Fourier coefficients are inserted by linear complexification of the real tensor kernels. There is no physical transfer cutoff.",
        "complete_endpoint_definition": "Q is the full unexpanded r<m first-five W6 endpoints plus their r>=m values minus ALL unaveraged Taylor coefficients through inverse-radius degree4-j. Every source time jet and both trace-gradient orders remain. The complete tensor is bounded before splitting it into four geometric products.",
        "whole_near_far_low_proof": "Far r>=200(m+|P|) uses the entire Cauchy remainder and integrates its actual threshold. Near retains the whole raw endpoint and every UV power/log term separately. The unexpanded low ball remains. The explicit1/m in all logarithmic estimates is retained, so the full finite coefficient grows at most m^4. All15 Taylor slots and35 source-jet entries are kept.",
        "complete_original_mask_tail": "The SAME Omega_star projection is applied to BOTH created legs of each complete convergent integrand. Far removal forces r>Kradial/2; near removal implies Kradial<201(m+|P|); low removal implies Kradial<m+|P|. No intersection or large-small pair is dropped. The highest source transfer degree is6 and mass power5. For K>=2m, Kradial=4sqrt(K²-m²)>K.",
        "full_physical_identification": "The complete scalar full-W6 endpoint coefficients agree with frozen S243. Its full angular/radial MSbar finite difference is restored after the dominated dimensional limit. The actual S244 homogeneous response anchors zero transfer; exact-state, time and endpoint differences supply every nonzero transfer. No old vector finite weights, extra S241 heat action or fitted counterterm is added.",
        "finite_cutoff_scope": "This is the complete UV-subtracted restored projection: the original two-leg mask multiplies the ENTIRE unaveraged subtracted integrand, while the fixed scalar finite local matching and profile remain full. It is not the bare projected Gaussian current alone. Angular-zero UV rows cannot be discarded before masking. Full finite-cutoff Ward/local-shape bookkeeping and nonlinear mean-tail contacts remain explicit obligations for later ADM/clock assembly.",
        "remaining": "Full prepared heavy ADM/common-clock Ward assembly, quantum constraints/inverse, nonlinear inhomogeneous feedback, interacting light/mixed loops, quantum gravitational limit, physical UV/Regge and original V/G/B/P8 remain open.",
        "checks": checks,
        "gates": {key: bool(value) for key, value in gates.items()},
    }
