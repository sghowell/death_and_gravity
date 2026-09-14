"""Correct physical onepoint contacts and explicit limits of the old replay."""

from functools import cache

import sympy as s
from p8_vacuum_affine_gauge_mean_transport import crossing as old_crossing
from p8_vacuum_affine_gauge_mean_transport import geometry as old_geometry
from p8_vacuum_affine_gauge_mean_transport import modes
from p8_vacuum_affine_nonlinear_auxiliary_measure import canonical

from . import source

ALPHA = s.Rational(9, 4)


@cache
def whole_volume_jets():
    v = s.Symbol("log_hat_volume", real=True)
    a = s.Symbol("homogeneous_hat_scale", positive=True)
    N, R = source.N, source.R
    F = a**3 * s.exp(3 * v) * canonical.U
    gradient = s.Matrix([s.diff(F, v), s.diff(F, N)])
    hessian = s.hessian(F, (v, N))
    slope = -s.Rational(3, 4) * s.diff(R, N) / R
    second = (
        s.Rational(21, 16) * s.diff(R, N) ** 2 / R**2
        - s.Rational(3, 4) * s.diff(R, N, 2) / R
    )
    return {
        "whole_physical_spatial_volume_density": F,
        "whole_physical_spatial_volume_gradient": gradient,
        "whole_physical_spatial_volume_Hessian": hessian,
        "whole_actual_logarithmic_lapse_slope": slope,
        "whole_actual_normalized_second_lapse_derivative": second,
        "checks": {
            "whole_actual_first_source_derivatives": (
                gradient - F * s.Matrix([3, slope])
            ).applyfunc(s.factor),
            "whole_actual_second_source_derivatives": (
                hessian - F * s.Matrix([[9, 3 * slope], [3 * slope, second]])
            ).applyfunc(s.factor),
            "whole_actual_parent_density_not_shifted_R": s.simplify(
                F / (a**3 * s.exp(3 * v)) - source.C ** s.Rational(3, 2)
            ),
        },
        "gates": {
            "full_R_first_and_second_lapse_derivatives_retained": hessian.has(
                s.diff(R, N)
            )
            and hessian.has(s.diff(R, N, 2)),
            "nonlinear_physical_source_ordering_not_selected": True,
        },
    }


@cache
def corrected_contact():
    Cvv, Cvn = s.symbols("whole_fixed_Cvv whole_fixed_Cvn", real=True)
    p = s.Symbol("scalar_lapse_density_power", real=True)
    slope = source.clock_jets()["whole_actual_volume_lapse_slope_at_bounce"]
    first = s.Matrix([-ALPHA * Cvv, -ALPHA * Cvn])
    dcov = s.Matrix(
        [[s.Rational(3, 2) * Cvv, s.Rational(3, 4) * Cvn], [s.Rational(3, 4) * Cvn, 0]]
    )
    grad = s.Matrix([3, slope])
    Hess = s.Matrix([[9, 3 * slope], [3 * slope, s.Rational(3, 4)]])
    onepoint = grad.dot(first)
    twopoint = s.trace(Hess * dcov) / 2
    wrong_slope = old_geometry.physical_volume_contact()[
        "whole_bounce_physical_volume_lapse_slope"
    ]
    wrong_onepoint = s.Matrix([3, wrong_slope]).dot(first)
    arbitrary_slope = s.Symbol("arbitrary_positive_C_log_slope", real=True)
    generic_cancel = (
        s.Matrix([3, arbitrary_slope]).dot(first)
        + s.trace(s.Matrix([[9, 3 * arbitrary_slope], [3 * arbitrary_slope, p]]) * dcov)
        / 2
    )
    packet = modes.packet()
    dv1 = packet["whole_first_log_volume_variation"]
    dv2 = packet["whole_second_log_volume_variation"]
    dn2 = packet["whole_second_lapse_variation"]
    v, n = (
        packet["whole_three_direction_scalar_fixture"],
        packet["whole_nonindependent_lapse_fixture"],
    )
    finite_one = 3 * modes.mean(dv2) + slope * modes.mean(dn2)
    finite_two = 9 * modes.mean(modes.mul(v, dv1)) + 3 * slope * modes.mean(
        modes.mul(n, dv1)
    )
    return {
        "whole_fixed_reference_mean_transport": first,
        "whole_fixed_reference_covariance_transport": dcov,
        "whole_correct_actual_volume_gradient_at_bounce": grad,
        "whole_correct_actual_volume_Hessian_at_bounce": Hess,
        "whole_correct_actual_onepoint_contact": onepoint,
        "whole_correct_actual_twopoint_contact": twopoint,
        "frozen_wrong_physical_lapse_slope": wrong_slope,
        "frozen_wrong_minus_actual_onepoint_contact": wrong_onepoint - onepoint,
        "whole_corrected_three_direction_finite_onepoint_contact": finite_one,
        "whole_corrected_three_direction_finite_twopoint_contact": finite_two,
        "whole_replay_limitation": "The generic cancellation is zero for every supplied logarithmic conformal slope, including the wrong one. It cannot select the physical parent. Independent binding to the original metric, U, Cchi and M is essential. Passing S261 algebra and replay tests does not endorse its false physical-map identification.",
        "checks": {
            "whole_corrected_physical_volume_contact_cancellation": s.expand(
                onepoint + twopoint
            ),
            "whole_corrected_onepoint_lapse_coefficient": s.expand(
                onepoint + s.Rational(27, 4) * Cvv + s.Rational(27, 8) * Cvn
            ),
            "whole_corrected_twopoint_lapse_coefficient": s.expand(
                twopoint - s.Rational(27, 4) * Cvv - s.Rational(27, 8) * Cvn
            ),
            "whole_frozen_wrong_minus_actual_slope": wrong_slope
            - slope
            + s.Rational(15, 2),
            "whole_frozen_wrong_minus_actual_onepoint": s.expand(
                wrong_onepoint - onepoint - s.Rational(135, 8) * Cvn
            ),
            "whole_generic_cancellation_does_not_select_actual_parent": s.expand(
                generic_cancel
            ),
            "whole_corrected_three_direction_contact": s.factor(
                finite_one + finite_two
            ),
        },
        "gates": {
            "frozen_actual_minus_six_identification_is_refuted": wrong_slope != slope,
            "wrong_physical_onepoint_is_not_actual_onepoint": wrong_onepoint
            != onepoint,
            "generic_identity_can_pass_with_wrong_parent_binding": generic_cancel == 0,
            "corrected_full_three_direction_contacts_are_nonzero": finite_one != 0
            and finite_two != 0,
            "coordinate_mean_and_canonical_covariance_not_reprepared": True,
        },
    }


@cache
def retained_symbols_and_erratum():
    old = old_crossing.bounce_symbols()
    cross = old["whole_R3_Cvn_Lambda_squared_coefficient"]
    lapse_mean = old["whole_gauge_mean_N_Lambda_squared_coefficient"]
    correct = s.Rational(3, 2) * lapse_mean
    wrong = -6 * lapse_mean
    return {
        "retained_whole_bounce_reference_symbols": {
            key: value for key, value in old.items() if key not in ("checks", "gates")
        },
        "whole_corrected_physical_volume_onepoint_quadratic_coefficient_from_lapse": correct,
        "whole_corrected_physical_volume_twopoint_quadratic_coefficient_from_lapse": -correct,
        "frozen_incorrect_onepoint_quadratic_coefficient_from_lapse": wrong,
        "explicit_erratum": {
            "frozen_report_sha256": "4ca9381026403ff0b759c68c56c4683a5abf2610ccb2a3d949dad6d2e63e557e",
            "refuted_identification": "S6.261 uses R-1/2 as the physical conformal factor and identifies c1(0)=-6 with the actual parent. Both identifications are false.",
            "correct_binding": "C=R**(-1/2)=Cchi**2=M**(-2), spatial volume U=R**(-3/4), c1(0)=3/2.",
            "retained_results": "The generic density-weight orbit, both-tensor periodic map, full Fourier displacement, Weyl coordinate mean shifts, canonical lapse row and commutator, and bounce coordinate covariance symbols do not use the erroneous physical C binding and are retained within their original limitations.",
            "not_retained": "The S261 numerical physical-volume lapse contact, its actual-parent labels and any derived identification of its -6-weighted contact with a physical source.",
            "immutability": "No S261 or earlier source/report byte is edited. The old matching payload is historical reproducibility evidence, explicitly overridden on this physical-binding issue by the current correction; its unchanged bytes are not renewed physical certification.",
        },
        "checks": {
            "retained_actual_cross_symbol_not_modified": cross
            - old["whole_R3_Cvn_Lambda_squared_coefficient"],
            "whole_corrected_lapse_onepoint_symbol": s.factor(
                correct + s.Rational(27, 8) * cross
            ),
            "whole_corrected_lapse_symbol_cancellation": correct - correct,
            "whole_wrong_minus_actual_lapse_symbol": s.factor(
                wrong - correct - s.Rational(135, 8) * cross
            ),
        },
        "gates": {
            "old_canonical_state_and_covariance_symbols_retained": True,
            "corrected_source_coefficient_not_new_quantum_mean": True,
            "old_frozen_report_reproducibility_is_not_physical_endorsement": True,
            "original_quantum_and_cutoff_obligations_unchanged": True,
        },
    }
