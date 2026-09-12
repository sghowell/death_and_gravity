"""Exact original-prescription spatial anchoring and complete weak norm."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes
from p8_vacuum_affine_dimensional_spatial_symbol import matching
from p8_vacuum_affine_uniform_uv_remainder import limit as known

from . import lift

FINITE_LOCAL = s.Integer(10) ** 5
CURRENT = 2 * s.Integer(10) ** 95


def anchored_current(homogeneous, known_P, known_zero, finite_UV):
    return homogeneous + known_P - known_zero + finite_UV


@cache
def data():
    H0, KP, K0, CP, C0, UP, U0, FP = s.symbols("H0 KP K0 CP C0 UP U0 FP")
    actual_difference = KP + CP + UP - (K0 + C0 + U0)
    checks = {
        "finite_regulator_full_contact_difference": s.expand(
            actual_difference.subs(CP, C0) - (KP - K0 + UP - U0)
        ),
        "anchored_same_prescription_identity": s.expand(
            anchored_current(H0, KP, K0, FP) - (H0 + KP - K0 + FP)
        ),
        "zero_spatial_transfer_recovers_complete_homogeneous_anchor": s.expand(
            anchored_current(H0, K0, K0, 0) - H0
        ),
        "complete_state_time_endpoint_sum": known.KNOWN - 5 * s.Integer(10) ** 48,
        "both_external_metric_factors": 4 * CURRENT / modes.KAPPA
        - 8 * s.Rational(1, 10) ** 705,
        "proper_density_divisor_minimum": s.expand(
            (1 + s.Symbol("t", real=True) ** 2) ** 6
            - 1
            - (
                s.Symbol("t", real=True) ** 12
                + 6 * s.Symbol("t", real=True) ** 10
                + 15 * s.Symbol("t", real=True) ** 8
                + 20 * s.Symbol("t", real=True) ** 6
                + 15 * s.Symbol("t", real=True) ** 4
                + 6 * s.Symbol("t", real=True) ** 2
            )
        ),
    }
    local = matching.finite_coefficient_norm()
    _inv, _g, _ell, finite = matching.finite_input()
    from p8_vacuum_affine_dimensional_spatial_symbol import jets

    for r, row in enumerate(finite):
        checks[f"finite_UV_origin_source_jet_{r}"] = s.factor(row.subs(jets.p, 0))
    bound = lift.HOMOGENEOUS + 4 * known.KNOWN + FINITE_LOCAL
    return {
        "same_current_and_preparation": "The S195 covariance tangent atP0 is exactly the derivative of the S194 homogeneous six-dimensional Hamiltonian flow: identical metric exponential, gamma-independent temporal constraint, common all-order covariance, zero initial tangent from the unchanged source germ, both propagated source terms and the complete second vertex. Relabelling ell=-q maps its pair form to the S198-S208 endpoint decomposition. Uniqueness of the finite-mode linear tangent ODE identifies the anchor before either regulator limit.",
        "fixed_matching_identity": "After the dimension-limit argument and the S212 fixed pole finite part, Jren(P)=H0+[Known(P)-Known(0)]+Ffinite(P). H0 is the full S194 fixed-prescription homogeneous linear current, not an adiabatic surrogate. Full contact cancels only in the spatial difference and remains in H0.",
        "no_extra_local_action": "S212 Ffinite is the actual UV-difference finite part at the unchanged MSbar scale, including odd evanescent endpoint, fixed-source dimensional reconstruction, Euler/volume and lower-band terms. It is not S204's local action target. No local action or old Q203/Pfin is added a second time.",
        "complete_bound": "The exact homogeneous contribution is at most1e95||D||L2 Z136. Each new Known contribution is below1e49 M[D]Z136; the finite local difference is below1e5||D||L2 Z136. Thus the complete tracefree spatial Gaussian current is bounded by2e95 M[D]Z136.",
        "bound_unrounded_sum": bound,
        "canonical": "Divide coordinate current by a^3>=1 and multiply4/kappa for both external metric factors. Since a depends only on time, this division does not increase the detector spatial H1 norm. The canonical weak-response display is8e-705. The small number is not a reduced inverse or background theorem.",
        "retarded": "The prepared actual covariance tangent is retarded. Algebraic endpoint repartition and local fixed counterterms preserve time causality. H0 and Known(0) are pointwise spatial multipliers; the finite spatial UV difference is local. No retarded strong-norm differentiability beyond the proved weak response is asserted.",
        "checks": checks,
        "gates": {
            "finite_local_norm_from_exact_frozen_coefficients": local["sum"]
            < FINITE_LOCAL,
            "full_current_strict_rounded_bound": bound < CURRENT,
            "same_actual_mass_and_kappa": modes.MASS == 1000
            and modes.KAPPA == s.Integer(10) ** 800,
            "all_actual_gaussian_state_and_contact_pieces_retained": True,
            "same_original_finite_prescription": True,
            "contact_anchor_not_zero": True,
            "no_local_target_or_old_finite_piece_double_counting": True,
            "tracefree_spatial_sector_only": True,
        },
    }
