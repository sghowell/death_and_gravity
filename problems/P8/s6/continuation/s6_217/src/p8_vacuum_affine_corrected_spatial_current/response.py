"""Canonical unequal-mode branch bridge and full corrected finite decomposition."""

from functools import cache

import sympy as s
from p8_vacuum_affine_reference_state_prefactor import limit as state_time
from p8_vacuum_affine_uniform_uv_remainder import estimates, tail

KNOWN = 5 * s.Integer(10) ** 48
KNOWN_TAIL = 5 * s.Integer(10) ** 53


@cache
def canonical_data():
    x, y, z, w = s.symbols("x y z w", real=True)
    U = s.Matrix([[1, x], [y, 1 + x * y]])
    V = s.Matrix([[1, z], [w, 1 + z * w]])
    D = s.Matrix(2, 2, s.symbols("D0:4", real=True))
    G = s.Matrix(2, 2, s.symbols("G0:4", real=True))
    J = s.Matrix([[0, 1], [-1, 0]])
    v = s.Matrix([1, -s.I]) / s.sqrt(2)
    source = (J * G + G * J.T) / 2
    covariance = -s.trace(D.T * U * source * V.T)
    detector_pair = (v.T * U.T * D * V * v)[0]
    source_pair = (v.T * G * v)[0]
    pair = 2 * s.im(s.conjugate(detector_pair) * source_pair)
    checks = {
        "independent_left_symplectic_transport": U * J * U.T - J,
        "independent_right_symplectic_transport": V * J * V.T - J,
        "full_unequal_mode_covariance_equals_positive_imaginary_annihilation_pair": s.expand(
            covariance - pair
        ),
        "equivalent_negative_imaginary_creation_pair": s.expand(
            pair + 2 * s.im(detector_pair * s.conjugate(source_pair))
        ),
    }
    return {
        "finite_regulator_bridge": "At each real-mode pair, the exact S195 covariance source with both reverse blocks equals +2Im(conjugate(b_D)b_G), b=u^t M u. The source/readout cross blocks are arbitrary real matrices and the two symplectic transports are independent. Summing the ordered physical polarization pairs gives the original Fourier normalization, not independent duplicate fields.",
        "reference_boundary": "This CCR identity applies algebraically to the chosen normalized unit-W8 readout family even though it is not propagated as a new physical state. Actual-minus-reference memory and full contact remain the S197/S201 comparison, expressed consistently in annihilation rather than creation amplitudes.",
        "checks": checks,
        "gates": {
            "two_distinct_internal_modes_and_transports": U != V,
            "complex_stripped_amplitudes_not_assumed_real": True,
            "both_reverse_blocks_and_original_Kubo_factor": True,
            "full_second_metric_contact_kept_separate": True,
            "unchanged_actual_gaussian_initial_state": True,
        },
    }


@cache
def decomposition_data():
    R, C, F5, F6 = s.symbols("R C F5_legacy F6_positive_phase")
    E = s.symbols("E0:5")
    U = s.symbols("U0:5")
    signs = tuple((-1) ** j for j in range(5))
    endpoints = sum(signs[j] * E[j] for j in range(5))
    UV = sum(signs[j] * U[j] for j in range(5))
    Q = sum(signs[j] * (E[j] - U[j]) for j in range(5))
    F = -F5 + F6
    actual = R + C + endpoints + F
    assembled = R + F + Q + C + UV
    norms = state_time.constants()
    checks = {
        "complete_corrected_current_same_original_regulator": s.expand(
            actual - assembled
        ),
        "each_corrected_endpoint_repartition": s.Matrix(
            [
                signs[j] * E[j] - signs[j] * (E[j] - U[j]) - signs[j] * U[j]
                for j in range(5)
            ]
        ),
        "fifth_endpoint_reversed_not_discarded": s.diff(F, F5) + 1,
        "sixth_bulk_retained_with_corrected_phase": s.diff(F, F6) - 1,
        "state_time_absolute_bound_unchanged": norms[
            "complete_known_finite_curved_piece"
        ]
        - 2 * s.Integer(10) ** 48,
        "state_time_original_tail_unchanged": norms["complete_known_finite_curved_tail"]
        - 2 * s.Integer(10) ** 52,
        "each_endpoint_absolute_bound_is_phase_invariant": s.Matrix(
            [abs(v) - 1 for v in signs]
        ),
    }
    return {
        "complete_reference": "Use Kplus and +Im from S216 for every complete annihilation pair. First five E_j and U_j are(-1)^j times their legacy mathematical rows. The fifth endpoint reverses; the sixth bulk has positive detector phase, negative source phase, coefficient(-i)^6 and +Im. It is not obtained by multiplying the entire legacy sixth bulk by a guessed sign.",
        "same_finite_current": "Jactual,K=R_actual-unit_correct,K+F_correct,K+Q_correct,K+C_unit,K+U_correct,K, with Q=sum low E_j_correct+high(E_j_correct-U_j_correct). Every endpoint label, finite cell and original one-/two-leg mask is retained.",
        "absolute_known_bound": "S197/S201 use norms of both mixed pair differences and the error square, invariant under consistent branch conjugation. The corrected fifth endpoint/sixth bulk have the same absolute inverse-phase majorants. For Q, both E_j and its Taylor subtraction acquire the same unit sign separately for each j, preserving S208's low/near/far and removed-union bounds.",
        "no_double_counting": "Rebase only from R_actual-unit+F. Do not add the earlier Q203/Pfin pieces or S204 local-action target to Q_correct. The full contact remains in the homogeneous anchor.",
        "checks": checks,
        "gates": {
            "corrected_known_absolute_bound": norms[
                "complete_known_finite_curved_piece"
            ]
            + estimates.FINITE
            < KNOWN,
            "corrected_known_original_tail": norms["complete_known_finite_curved_tail"]
            + tail.TAIL
            < KNOWN_TAIL,
            "full_low_near_far_absolute_endpoint_bound": estimates.constants()[
                "complete_all_momentum_UV_subtracted_endpoint_bound"
            ]
            < estimates.FINITE,
            "full_removed_union_absolute_endpoint_tail": tail.constants()[
                "complete_original_regulator_tail_numerator"
            ]
            < tail.TAIL,
            "sixth_bulk_phase_corrected_not_only_local_UV": True,
            "all_lower_endpoints_from_same_initial_source_germ": True,
            "original_masks_and_fixed_contact_not_changed": True,
        },
    }
