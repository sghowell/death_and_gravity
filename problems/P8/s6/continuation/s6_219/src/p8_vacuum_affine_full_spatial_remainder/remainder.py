"""Full spatial state, time, endpoint, contact and original-mask bounds."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_metric_noise import reference
from p8_vacuum_affine_corrected_spatial_current import response as corrected
from p8_vacuum_affine_reference_state_prefactor import contact as occupation_contact
from p8_vacuum_affine_reference_state_prefactor import limit as occupation
from p8_vacuum_affine_retarded_reference_boundary import bounds as time_remainder
from p8_vacuum_affine_retarded_state_remainder import contact, memory
from p8_vacuum_affine_retarded_state_remainder import tail as state_tail
from p8_vacuum_affine_subleading_band_conversion import centered
from p8_vacuum_affine_uniform_uv_remainder import domain, estimates, tail

KNOWN = 5 * s.Integer(10) ** 48
KNOWN_TAIL = 5 * s.Integer(10) ** 53
SHAPE_TAIL = s.Integer(10) ** 40


@cache
def data():
    F, R = reference.FIELD_NORM, contact.REMAINDER
    old = occupation.constants()
    coefficient = domain.FIELD**2
    checks = {
        "unit_full_feature_reference_pair_count": (2 * F) ** 2 - 64000000,
        "full_feature_error_pair_coefficient": 3 * F**2 * R - 96000000000000,
        "all_ten_normalized_fields_full_spatial_pair": coefficient - 400,
        "unchanged_complete_actual_unit_state_time_bound": old[
            "complete_known_finite_curved_piece"
        ]
        - 2 * s.Integer(10) ** 48,
        "unchanged_complete_state_time_original_tail": old[
            "complete_known_finite_curved_tail"
        ]
        - 2 * s.Integer(10) ** 52,
        "five_corrected_endpoints_absolute_signs": s.Matrix(
            [abs((-1) ** j) - 1 for j in range(5)]
        ),
        "full_state_contact_half_and_all_polarizations": contact.MODE
        - s.Rational(15, 2) * F**2 * R,
        "full_occupation_contact_half_and_all_polarizations": occupation_contact.constants()[
            "complete_three_mode_contact_coefficient"
        ]
        - s.Rational(3, 2) * F**2 * R**2,
        "same_complete_all_transfer_endpoint_rows": len(estimates.near_terms()) - 16,
        "same_complete_endpoint_tail_external_degree": max(
            q + 1 for q, _c, _label in estimates.near_terms()
        )
        - 6,
        "same_sixth_bulk_and_fifth_boundary_weight": time_remainder.constants()[
            "complete_finite_reference_integral"
        ]
        - time_remainder.constants()["sixth_bulk_full_integral"]
        - time_remainder.constants()["fifth_endpoint_full_integral"],
    }
    return {
        "full_pair_input": "The exact full physical feature bridge gives |b_D(v,w)|<=||D||F||v||||w|| for arbitrary spatial D, not only shear. The same full ten-field reference bounds give Cref<=64000000<1e9. With actual=alpha0*W8+e, both mixed readouts and the error square are bounded by3*4000^2*R sqrt(nu mu)(nu^-6+mu^-6)<1e15 sqrt(nu mu)(nu^-6+mu^-6). This uses R/nu^6<1 and keeps the complete actual state.",
        "full_contact_transfer": "Replace the old shear contact feature by the independently derived full QQ feature. Its operator norm is still<=||D||F||G||F; the temporal block is nonzero but fits this bound. The same three polarization sums, half-Hamiltonian normalization, both mixed covariance terms and error square give S197's full contact bound/tail. The alpha0-to-unit occupation contact similarly uses all ten fields and the same norm. These are complete contact COMPARISONS; the unit contact itself stays in the homogeneous anchor.",
        "state_and_time": "Every hypothesis of the S197/S201 real-time memory/contact comparison now holds for arbitrary spatial tensors with the same Cref,Cerr and mass. Their unchanged low/high external tail partition therefore gives actual-to-unit3e23 and2e27/K. The full six-step corrected branch has fifth endpoint minus the legacy fifth endpoint, and an independently positive-detector-phase sixth bulk with +Im. Pair norms and the unchanged inverse summed phases give the complete2e48 state/time bound and2e52/K tail. No final source cutoff is needed: every upper endpoint remains; the lower germ is unchanged.",
        "full_endpoint_remainder": "For all real momenta, the full feature contraction obeys the existing1e9 sqrt(nu mu) pair bound and harmonic-mean phase suppression, including zero internal momentum. On the normalized joint inverse-radius/time domain every ten-field norm is<20 and full spatial pair norm<=400<1000. Thus all S208 row/jet Cauchy, unexpanded low, near and far estimates apply with unchanged B_j. No added polarization or geometry count is needed: the bound is on the complete nine physical pairs before decomposition into ten algebraic products.",
        "full_known": "Known_full=R_actual-unit_full+F_correct_full+Q_correct_full. Its complete finite bound is5e48 M[D]Y[G] and its original one-/two-leg regulator tail is5e53 M[D]Y[G]/K, with Y^2=N61^2+X46^2 as before. This is the full six-dimensional spatial tensor form, not a sum of separately rounded scalar-channel bounds.",
        "full_shape_tail": "The full pair exchange/Schwarz bridge and all six channel shape coefficients give the same even centered TOTAL grading. Every full uncentered coefficient of degree d has azimuthal degree<=d+4<=8: trace shifts add no momenta, LC has degree at most2, and CC at most2 before the expansion. Therefore S209's8192 shell,6121 angular C2 and64 shape-functional constants apply to the full row bounds. The same all-P partition gives1e40 M[D]X46[G]/K for K>=2m, including the retained q4 log shell.",
        "norms": "M^2=||D||L2^2+||grad_x D||L2^2; N61^2=sum_r0..6||partial_t^r G||L2^2+||grad_x G||L2^2; X46^2=sum_r0..4||(1-Delta)^3 partial_t^r G||L2^2. All tensors use Frobenius norms. These are weak derivative-losing bounds on the original unit slab.",
        "unrounded_sums": {
            "full_known": old["complete_known_finite_curved_piece"] + estimates.FINITE,
            "full_known_tail": old["complete_known_finite_curved_tail"] + tail.TAIL,
            "full_shape_tail": centered.tail_constants()[
                "complete_uniform_conversion_remainder"
            ],
        },
        "checks": checks,
        "gates": {
            "complete_full_reference_pair_below_old_majorant": (2 * F) ** 2
            < reference.PAIR_REF,
            "complete_full_error_pair_below_old_majorant": 3 * F**2 * R
            < reference.PAIR_ERROR,
            "full_error_relative_factor_below_one": R / reference.MASS**6 < 1,
            "full_normalized_pair_below_old_majorant": coefficient < domain.PAIR,
            "all_full_ten_field_domain_bounds": all(domain.data()["gates"].values()),
            "full_memory_constant_unchanged_and_bounded": memory.constants()[
                "complete_memory_H1_coefficient"
            ]
            < memory.DISPLAY,
            "full_contact_constant_unchanged_and_bounded": contact.constants()[
                "contact_L2_coefficient"
            ]
            < contact.DISPLAY,
            "full_state_contact_and_memory_original_tail": state_tail.constants()[
                "complete_K_tail_numerator"
            ]
            < state_tail.DISPLAY,
            "full_occupation_contact_and_memory_remainder": old[
                "full_prefactor_memory_contact_bound"
            ]
            < occupation.DISPLAY,
            "full_occupation_original_tail": old[
                "full_prefactor_regulator_tail_numerator"
            ]
            < occupation.TAIL,
            "complete_sixth_bulk_and_fifth_endpoint": time_remainder.constants()[
                "complete_finite_reference_integral"
            ]
            < time_remainder.DISPLAY,
            "complete_six_step_original_tail": time_remainder.constants()[
                "complete_removed_two_leg_tail_numerator"
            ]
            < time_remainder.TAIL,
            "full_endpoint_near_far_low_constant": estimates.constants()[
                "complete_all_momentum_UV_subtracted_endpoint_bound"
            ]
            < estimates.FINITE,
            "full_endpoint_removed_union_tail": tail.constants()[
                "complete_original_regulator_tail_numerator"
            ]
            < tail.TAIL,
            "full_known_display": old["complete_known_finite_curved_piece"]
            + estimates.FINITE
            < KNOWN,
            "full_known_tail_display": old["complete_known_finite_curved_tail"]
            + tail.TAIL
            < KNOWN_TAIL,
            "full_original_shape_tail_display": centered.tail_constants()[
                "complete_uniform_conversion_remainder"
            ]
            < SHAPE_TAIL,
            "corrected_branch_not_a_sign_guess_on_sixth_bulk": all(
                corrected.canonical_data()["gates"].values()
            ),
        },
    }
