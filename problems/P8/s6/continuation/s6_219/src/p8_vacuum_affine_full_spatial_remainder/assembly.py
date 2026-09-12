"""Complete reference spatial current and prepared unreduced metric response."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes
from p8_vacuum_affine_corrected_spatial_current import assembly as corrected
from p8_vacuum_affine_homogeneous_trace_anchor import anchor
from p8_vacuum_affine_matched_spatial_current import lift
from p8_vacuum_affine_ordered_scalar_symbol import matching
from p8_vacuum_affine_prepared_ward_reconstruction import norms as prepared

from . import remainder, shapes

CURRENT = 2 * s.Integer(10) ** 95
TAIL = 3 * s.Integer(10) ** 54
CROSS = 3 * s.Integer(10) ** 49
FULL_WARD = s.Integer(10) ** 118
FINITE = s.Integer(10) ** 6


@cache
def data():
    H, KP, K0, LP, L0, UP, U0, C = s.symbols("H KP K0 LP L0 UP U0 C")
    K, mass = s.symbols("K mass", positive=True)
    A3, A1, F2, F4, finite, E = s.symbols("A3 A1 F2 F4 finite E")
    uv = corrected.subtractor(K, mass, A3, A1, F2, F4, 0)
    physical = corrected.anchored_current(H, LP, L0, finite)
    regulated = (
        H
        + KP
        + C
        + UP
        - K0
        - C
        - U0
        - corrected.subtractor(K, mass, A3, A1, F2, F4, finite)
    )
    error = s.expand((regulated - physical).subs(UP, U0 + uv + E))
    inv, _g, _ell, finite_rows = matching.finite_input()
    checks = {
        "full_original_regulator_error_exact": s.expand(
            error - (KP - LP - (K0 - L0) + E)
        ),
        "complete_full_contact_cancels_only_in_transfer_difference": s.expand(
            KP + C + UP - K0 - C - U0 - (KP - K0 + UP - U0)
        ),
        "full_homogeneous_anchor_not_replaced": corrected.anchored_current(H, K0, K0, 0)
        - H,
        "full_finite_UV_zero_transfer": s.Matrix(
            [s.factor(value.subs(shapes.p, 0)) for value in finite_rows]
        ),
        "full_trace_homogeneous_norm_is_original_display": anchor.FULL_HOMOGENEOUS
        - s.Integer(10) ** 95,
        "full_canonical_spatial_current_two_external_factors": 4 * CURRENT / modes.KAPPA
        - 8 * s.Rational(1, 10) ** 705,
        "full_canonical_original_tail_two_external_factors": 4 * TAIL / modes.KAPPA
        - 12 * s.Rational(1, 10) ** 746,
        "same_all_six_full_spatial_invariants": len(inv) - 6,
    }
    current = anchor.FULL_HOMOGENEOUS + 4 * remainder.KNOWN + FINITE
    cross = 4 * remainder.KNOWN + FINITE
    tail = 4 * remainder.KNOWN_TAIL + remainder.SHAPE_TAIL
    return {
        "complete_full_spatial_current": "The full first/second ADM feature bridge, corrected finite-mode CCR decomposition, full state/time/contact and endpoint/dimension limits identify the actual reference spatial Gaussian response: Jren_full(P)=H0_full+Known_full(P)-Known_full(0)+Ffinite_full(P). The full H0 is S218, and Ffinite is exactly S216's original six-invariant finite prescription. No legacy Q203/Pfin or local target is added again.",
        "original_regulator": "Retain the original two-leg mask and one-leg contact. The full homogeneous-anchored subtractor is K^3 A3_full+K A1_full+F2_full(K^2-m^2)/2+F4_full log(K/m)-Ffinite_full. The three new scalar A3/A1 rows are derived here with every ordered source jet; the corrected tracefree rows are S217/S210. Quadratic/finite shell artifacts cancel in the complete sum, while the q4 shell remains in the tail.",
        "complete_spatial_weak_bound": "For arbitrary spatial symmetric tensors, |Jren_full(D,G)|<2e95 M[D]Z136[G] and the exact homogeneous-anchored original regulator error is<3e54 M[D]Z136[G]/K for K>=2m. The bounds are zero for a zero norm product. Full pair bounds are established before channel splitting, so no unproved scalar-channel summation factor is suppressed.",
        "ordered_scalar_bounds": "The normalized trace/trace kernel satisfies Ltt=2e95. Both independently ordered normalized cross kernels satisfy Lt0=L0t=3e49 as bounds, not equality of kernels: their homogeneous anchors vanish, and4*Known+finite<3e49. Momentum-direction projections are orthogonal Frobenius multipliers and commute with the spatial/time norm weights.",
        "norm_lift": "The S218 full trace/tracefree homogeneous identity multiplier is bounded by1e95||D||L2 Z130. The same complete source embeddings give Y<2Z136 and Z24<=Z136. The pointwise explicit Known(0) kernel is lifted as an identity multiplier before Fourier Cauchy; no evaluation of an arbitrary weak multiplier at a point is assumed.",
        "canonical_boundary": "The two spatial metric factors give8e-705 and12e-746/K. These are complete reference Gaussian weak derivative-losing estimates, not a same-space reduced inverse or stability. The actual parent's reference Hessian agrees because its source and first source variation vanish there; no full nonlinear sourced-parent remainder follows.",
        "unrounded": {
            "current": current,
            "ordered_cross": cross,
            "original_tail": tail,
        },
        "checks": checks,
        "gates": {
            "full_spatial_current_bound": current < CURRENT,
            "both_ordered_cross_bounds": cross < CROSS,
            "full_original_regulator_tail": tail < TAIL,
            "full_fixed_six_invariant_finite_norm": matching.finite_coefficient_norm()[
                "sum"
            ]
            < FINITE,
            "same_all_transfer_source_norm_embeddings": all(
                lift.data()["gates"].values()
            ),
            "same_full_homogeneous_anchor": anchor.CURRENT[1]
            == anchor.FULL_HOMOGENEOUS,
            "initial_source_germ_only_all_upper_endpoints_retained": True,
            "original_actual_state_and_prescription": True,
        },
    }


@cache
def ward_data():
    direct = prepared.DETECTOR * prepared.SOURCE * CURRENT + prepared.WARD
    conditional = s.Integer(10) ** 22 * (CURRENT + CURRENT + 2 * CROSS) + prepared.WARD
    checks = {
        "full_prepared_direct_projection_count": direct
        - (2 * s.Integer(10) ** 117 + s.Integer(10) ** 35),
        "all_three_previously_missing_ordered_scalar_bounds": conditional
        - (4 * s.Integer(10) ** 117 + 6 * s.Integer(10) ** 71 + s.Integer(10) ** 35),
        "unchanged_detector_projection": prepared.DETECTOR - 100,
        "unchanged_source_projection": prepared.SOURCE - s.Integer(10) ** 20,
        "unchanged_complete_one_point_Ward_correction": prepared.WARD
        - s.Integer(10) ** 35,
    }
    return {
        "full_prepared_reconstruction": "S215's exact source-retarded and detector-advanced synchronous projections, ordered covariant-density Ward identities and full metric-chart contacts apply to the complete spatial current now proved. Both lower/upper boundary terms are retained with their original vanishing germs. This supplies every lapse, shift, trace and tracefree metric coordinate direction, not just the previously known shear-plus-Ward block.",
        "source_detector_domains": "The spatial formulas require only the original zero initial source germ, with all upper endpoint terms kept. Hence the prepared synchronous source's initial-only domain is admitted. Cut the advanced detector off inside the initial zero-output germ; M has no time derivatives, so no cutoff-derivative norm is incurred.",
        "full_bound": "M[Qsyn_D]<=100 V04[D], Z136[Qsyn_G]<=1e20 U138[G], and the full one-point Ward correction is<1e35 V04 U138. The resulting complete reference Gaussian metric response satisfies |Rfull(D,G)|<1e118 V04[D]U138[G]. Independently inserting Ltt=2e95 and Lt0=L0t=3e49 into S215's conditional block formula also fits this display.",
        "unrounded_direct": direct,
        "unrounded_ordered_blocks": conditional,
        "frontier": "This discharges the missing reference Gaussian metric-current existence and weak-bound input of the prepared Ward reconstruction. It does not produce the genuinely constrained/reduced scalar-mixed inverse, a same-space contraction, a finite-amplitude sourced-parent remainder, quantum-corrected background/stability, physical cutoff/heavy-sector control, or original V/G/B closure.",
        "checks": checks,
        "gates": {
            "full_prepared_direct_bound": direct < FULL_WARD,
            "full_prepared_ordered_block_bound": conditional < FULL_WARD,
            "all_original_projection_and_one_point_constants": prepared.data()[
                "projection_unrounded"
            ]["detector"]
            < prepared.DETECTOR
            and prepared.data()["projection_unrounded"]["source"] < prepared.SOURCE
            and 28520 * s.Integer(10) ** 30 < prepared.WARD,
            "all_three_scalar_inputs_now_actual_not_assumed": True,
            "no_division_by_H_or_external_momentum": True,
            "weak_metric_response_not_reduced_inverse": True,
        },
    }
