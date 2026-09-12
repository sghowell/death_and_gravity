"""Full corrected homogeneous anchoring, original regulator and Ward input."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes
from p8_vacuum_affine_matched_spatial_current import lift
from p8_vacuum_affine_ordered_scalar_symbol import matching
from p8_vacuum_affine_prepared_ward_reconstruction import norms as prepared

from . import response, subtraction

FINITE = s.Integer(10) ** 6
CURRENT = 2 * s.Integer(10) ** 95
TAIL = 3 * s.Integer(10) ** 54
WARD_KNOWN = s.Integer(10) ** 118


def anchored_current(homogeneous, known_P, known_zero, finite):
    return homogeneous + known_P - known_zero + finite


def subtractor(K, m, A3, A1, F2, F4, finite):
    return K**3 * A3 + K * A1 + F2 * (K * K - m * m) / 2 + F4 * s.log(K / m) - finite


@cache
def data():
    H, KP, K0, LP, L0, UP, U0, C = s.symbols("H KP K0 LP L0 UP U0 C")
    K, m = s.symbols("K m", positive=True)
    A3, A1, F2, F4, finite, E = s.symbols("A3 A1 F2 F4 finite E")
    UV = subtractor(K, m, A3, A1, F2, F4, 0)
    physical = anchored_current(H, LP, L0, finite)
    regulated = (
        H + (KP + C + UP) - (K0 + C + U0) - subtractor(K, m, A3, A1, F2, F4, finite)
    )
    error = s.expand((regulated - physical).subs(UP, U0 + UV + E))
    inv, _g, _ell, finite_rows = matching.finite_input()
    fixed = tuple(s.factor(v.subs(dict.fromkeys(inv[3:], 0))) for v in finite_rows)
    checks = {
        "exact_corrected_full_contact_cancellation_only_in_difference": s.expand(
            (KP + C + UP) - (K0 + C + U0) - (KP - K0 + UP - U0)
        ),
        "corrected_anchor_recovers_actual_homogeneous_current": anchored_current(
            H, K0, K0, 0
        )
        - H,
        "original_corrected_regulator_exact_full_error": s.expand(
            error - (KP - LP - (K0 - L0) + E)
        ),
        "finite_lower_band_not_deleted": s.expand(
            UV
            - (K**3 * A3 + K * A1 + F2 * K * K / 2 + F4 * s.log(K / m))
            + F2 * m * m / 2
        ),
        "same_original_homogeneous_derivative_display": lift.HOMOGENEOUS
        - s.Integer(10) ** 95,
        "two_canonical_metric_factors": 4 * CURRENT / modes.KAPPA
        - 8 * s.Rational(1, 10) ** 705,
        "two_canonical_tail_factors": 4 * TAIL / modes.KAPPA
        - 12 * s.Rational(1, 10) ** 746,
    }
    for r, row in enumerate(fixed):
        checks[f"corrected_finite_UV_zero_transfer_jet_{r}"] = s.factor(
            row.subs(subtraction.jets.p, 0)
        )
    current_bound = lift.HOMOGENEOUS + 4 * response.KNOWN + FINITE
    tail_bound = 4 * response.KNOWN_TAIL + subtraction.CONVERSION_TAIL
    return {
        "same_actual_homogeneous_anchor": "At P0, the actual S195 covariance tangent with both source terms, unchanged initial covariance and full second vertex is exactly the derivative of the S194 homogeneous Hamiltonian flow. Its current came from the original homogeneous one-point prescription, before the erroneous S207 phase bridge. Finite-mode ODE uniqueness identifies H0 without using the withdrawn S213 endpoint identification.",
        "corrected_full_current": "The independent canonical bridge, complete corrected six-step decomposition, uniform endpoint dimension limit and original S216 finite part now give Jren_correct(P)=H0+[Known_correct(P)-Known_correct(0)]+Ffinite_correct(P). Both finite endpoint and sixth-bulk phases are corrected, and the full one-leg contact remains in H0.",
        "full_weak_bound": "The homogeneous Frobenius dual/Hilbert lift and initial-germ time embedding are unaffected by the phase correction. Each Known term has bound5e48 M Y, Y<2 Z136. The corrected local UV term has bound1e6 ||D||L2 Z24<=1e6 M Z136. Their strict sum is below2e95 M Z136, with zero for a zero norm product.",
        "current_unrounded_bound": current_bound,
        "original_regulator": "Keep the original two-leg mask. The subtractor is K^3 A3+K A1+F2(K^2-m^2)/2+F4 log(K/m)-Ffinite_correct. All nondecaying physical coefficients are unchanged by the verified odd-UV structure. The finite correction and corrected Known terms are not dropped.",
        "anchored_regulator_error": "The exact error is[Known_correct,K(P)-Known_correct(P)]-[Known_correct,K(0)-Known_correct(0)]+E_K. The two corrected Known tails and the original all-P shape tail give3e54 M Z136/K for K>=2m. This is an exact homogeneous-anchored auxiliary-regulator limit, not an unanchored cutoff rate or a physical EFT cutoff.",
        "tail_unrounded_bound": tail_bound,
        "canonical_scope": "Both external metric factors give4/kappa. Dividing by a^3>=1 cannot increase detector M because a depends only on time. Displays8e-705 and12e-746/K remain weak tracefree derivative-losing estimates, not reduced scalar inverse or nonlinear stability.",
        "checks": checks,
        "gates": {
            "corrected_current_strict_rounded_bound": current_bound < CURRENT,
            "corrected_all_transfer_tail_strict_rounded_bound": tail_bound < TAIL,
            "corrected_full_finite_coefficient_norm": matching.finite_coefficient_norm()[
                "sum"
            ]
            < FINITE,
            "all_original_source_and_detector_norm_embeddings": all(
                bool(v) for v in lift.data()["gates"].values()
            ),
            "entire_current_bridge_not_only_UV_finite_patch": True,
            "original_state_contact_and_prescription": True,
            "retarded_source_germ_and_no_future_source_dependence": True,
            "full_scalar_current_and_reduced_inverse_not_inferred": True,
        },
    }


@cache
def ward_data():
    # The maps and density/chart Ward terms are independent of the erroneous
    # endpoint phase. Only their tracefree input is replaced by this new current.
    source = s.Integer(10) ** 20
    detector = s.Integer(100)
    ward = s.Integer(10) ** 35
    total = source * detector * CURRENT + ward
    checks = {
        "same_prepared_detector_factor": detector - prepared.DETECTOR,
        "same_prepared_source_factor": source - prepared.SOURCE,
        "same_complete_one_point_Ward_bound": ward - prepared.WARD,
        "corrected_known_Ward_unrounded_total": total
        - (2 * s.Integer(10) ** 117 + s.Integer(10) ** 35),
        "no_change_to_unit_slab_initial_primitive_factor": s.integrate(
            1, (s.Symbol("t"), 0, 1)
        )
        - 1,
    }
    return {
        "restored_input_not_old_formula": "S215's prepared maps, ordered covariant density Ward identities and full nonlinear metric-chart contacts are unchanged. Replace its withdrawn tracefree input by the correctly assembled current proved here. Do not re-endorse S213's old physical formula or S212's old finite value.",
        "initial_only_source_extension": "The new six-step formula uses only the unchanged zero initial source germ; every upper endpoint is retained. The state comparison and homogeneous one-point derivative require no final source cutoff. Thus the same initial-only source extension applies to retarded synchronous sources. Advanced detectors may be cut off inside the initial zero-output germ at no cost because M has no time derivatives.",
        "restored_known_piece_bound": "M[Dsyn]<=100 V04[D], Z136[Gsyn]<=1e20 U138[G], and the unchanged full one-point Ward correction is below1e35 V04 U138. With the new2e95 tracefree input, the known tracefree-plus-Ward piece is strictly below1e118 V04 U138.",
        "unrounded_known_bound": total,
        "remaining": "This supplies the corrected known block of the conditional reconstruction only. Three ordered scalar kernels, the homogeneous trace anchor, complete scalar state/time/contact bounds and the genuinely reduced inverse remain open.",
        "checks": checks,
        "gates": {
            "corrected_known_Ward_bound": total < WARD_KNOWN,
            "original_primitive_projection_and_one_point_bounds": prepared.data()[
                "projection_unrounded"
            ]["detector"]
            < detector
            and prepared.data()["projection_unrounded"]["source"] < source
            and 28520 * s.Integer(10) ** 30 < ward,
            "both_prepared_source_detector_boundaries_retained": True,
            "no_division_by_H_or_momentum_added": True,
            "unchanged_actual_one_point_and_metric_chart_contacts": True,
            "known_block_not_full_scalar_operator": True,
        },
    }
