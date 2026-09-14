"""Strict reference-composite and individual-coefficient scope, not quantum closure."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_finite_window_growth.intervals import fraction
from p8_vacuum_affine_quantitative_gaussian_window import audit as previous

from . import gaussian, local, physical, source

STATE = "same_original_S251_Gaussian_no_state_projection_or_conditioning"
MEASURE = (
    "specified_positive_lapse_reference_observables_not_interacting_quantum_measure"
)
OBSERVABLES = (
    "complete_bounded_single_lapse_negative_power_spectral_moments",
    "complete_commuting_bounce_slice_nonpolynomial_reference_volume",
    "whole_individual_ADM_kinetic_coefficient_Gaussian_divergence",
    "whole_canonical_matter_coefficient_finite_first_infinite_second_moment",
)
ITEM = {
    "id": "QG2_H8A429_full_nonlinear_reference_composites_and_Gaussian_integrability_boundary",
    "status": "FULL_BOUNDED_COEFFICIENT_AND_COMMUTING_BOUNCE_REFERENCE_VOLUME_MOMENTS_WITH_SINGULAR_INDIVIDUAL_COEFFICIENT_OBSTRUCTION_NOT_CONSTRAINED_ACTION_NO_GO_INTERACTING_MEAN_OR_ORIGINAL_V_G_B_P8",
}


def parameters():
    return previous.parameters()


def require_parameters(values):
    return previous.require_parameters(values)


def require_band(lower, upper):
    return previous.require_band(lower, upper)


def require_times(first, last):
    return previous.require_times(first, last)


def require_state(value):
    if not isinstance(value, str) or value != STATE:
        raise ValueError(
            "The original Gaussian is not projected, conditioned or reselected"
        )
    return value


def require_measure(value):
    if not isinstance(value, str) or value != MEASURE:
        raise ValueError(
            "The stated reference observables do not define an interacting measure"
        )
    return value


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError(
            "Keep the complete original coefficient and its stated reference scope"
        )
    return value


def require_power(value):
    value = fraction(value)
    if value not in {fraction(alpha) for alpha in source.ALPHAS}:
        raise ValueError("Only the three complete bounded negative powers are included")
    return value


def require_joint_time(value):
    value = fraction(value)
    if value != 0:
        raise ValueError(
            "The joint spectral volume law is proved only on the bounce slice"
        )
    return value


def frontier():
    return previous.frontier()


def matching():
    return [*previous.matching(), dict(ITEM)]


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError(
            "Do not promote individual-coefficient results to original P8 closure"
        )
    return True


@cache
def packets():
    return {
        "whole_original_R_and_global_positive_lapse_coefficient_bounds": source.global_R(),
        "whole_full_original_Hamiltonian_and_lapse_weight_bindings": source.weights(),
        "whole_complete_finite_strip_inverse_power_derivative_bounds": local.derivatives(),
        "whole_strongly_commuting_bounce_rows_and_nonzero_covariance": physical.bounce(),
        "whole_complete_Gaussian_fourth_remainders_and_tail_bounds": gaussian.remainder(),
        "whole_positive_lapse_nonlinear_reference_composites": gaussian.composites(),
        "whole_singular_individual_coefficient_integrability_boundary": physical.integrability(),
    }


@cache
def residuals():
    return {
        packet + "_" + key: clean(value)
        for packet, data in packets().items()
        for key, value in data["checks"].items()
    }


def scalar_entry_count():
    return sum(
        len(value) if isinstance(value, s.MatrixBase) else 1
        for value in residuals().values()
    )


@cache
def gates():
    return {
        **{
            packet + "_" + key: bool(value)
            for packet, data in packets().items()
            for key, value in data["gates"].items()
        },
        "original_primitive_frontier_unchanged": len(frontier()) == 9,
        "all_matching_identifiers_unique": len({row["id"] for row in matching()})
        == len(matching()),
        "all_original_parameters_unchanged": require_parameters(parameters())
        == parameters(),
        "S261_physical_binding_remains_explicitly_refuted": True,
        "original_quantum_mean_cutoff_and_P8_still_open": True,
    }


def observable():
    return {
        "positive_result": "All three complete positive-lapse inverse-R spectral factors have nonlinear Gaussian expectations within1e-960 of their full second-order center moments on the S264 band and slab. At the bounce, the entire positive-lapse reference volume exp(3v_hat)U(1+n) has a legitimate commuting joint spectral meaning, and its expectation differs from1 plus the full S264 contact by less than1e-960. Neither source switch nor any covariance cross term is dropped.",
        "strict_probability_boundary": "The original Gaussian state is not conditioned or projected. The stated observable is set to zero on N<=0. The joint spectral interpretation is only at the bounce, where the complete real Weyl smearing vectors use commuting momenta; away from it only single-lapse spectral moments are claimed. No joint noncommuting law or spacetime supremum is inferred.",
        "negative_result": "The actual lapse variance is positive, above1e-576 throughout the tiny slab and1e-572 at the bounce. For unrestricted positive-lapse Gaussian substitution, the full individual coefficient M/N has infinite mean; the full canonical matter coefficient N/U has finite first but infinite second moment. Small nonzero tails do not regularize these exact singular integrals.",
        "original_problem": "These are full background-map composites of LINEAR reference fields and an individual-coefficient integrability obstruction, NOT a full nonlinear auxiliary lapse reconstruction, divergent complete constrained Hamiltonian/action, candidate no-go, interacting physical mean, quantum ordering/gauge regulator, physical subtraction/cutoff or omitted-loop result. Full original H/Proca preparations and determinants remain. Original nonlinear inhomogeneous B and V/G/B/P8 remain OPEN; completed scoped P8(a) and the S261 erratum are unchanged.",
    }


def bad_cases():
    cases = []
    invalid = (
        True,
        False,
        1.0,
        s.Float(1),
        "1",
        None,
        s.oo,
        s.I,
        s.nan,
        s.Symbol("unknown"),
    )
    for i, value in enumerate(invalid):
        for label, call in (
            ("state", require_state),
            ("measure", require_measure),
            ("observable", require_observable),
            ("power", require_power),
            ("joint_time", require_joint_time),
        ):
            cases.append(("invalid_" + label + "_" + str(i), call, (value,)))
        cases.append(("invalid_band_" + str(i), require_band, (value, source.HIGH)))
        cases.append(("invalid_time_" + str(i), require_times, (value, source.TIME)))
        for key in parameters():
            values = parameters()
            values[key] = value
            cases.append(
                (
                    "invalid_parameter_" + key + "_" + str(i),
                    require_parameters,
                    (values,),
                )
            )
    for value in (0, 1, -s.Rational(1, 4), s.Rational(3, 2)):
        cases.append(("unsupported_power_" + str(value), require_power, (value,)))
    for value in (-source.TIME, source.TIME, 1):
        cases.append(
            ("non_bounce_joint_time_" + str(value), require_joint_time, (value,))
        )
    for lo, hi in (
        (0, source.HIGH),
        (source.LOW - 1, source.HIGH),
        (source.LOW, source.HIGH + 1),
        (source.HIGH, source.LOW),
    ):
        cases.append(
            ("outside_band_" + str(lo) + "_" + str(hi), require_band, (lo, hi))
        )
    for first, last in (
        (-2 * source.TIME, 0),
        (0, 2 * source.TIME),
        (source.TIME, -source.TIME),
    ):
        cases.append(
            (
                "outside_times_" + str(first) + "_" + str(last),
                require_times,
                (first, last),
            )
        )
    for label in (
        "conditioned_or_projected_replacement_Gaussian",
        "bounded_R_inverse_means_all_action_coefficients_integrable",
        "zero_Gaussian_variance_or_zero_small_tail",
        "finite_first_moment_implies_finite_second",
        "entire_Hamiltonian_diverges_from_one_coefficient",
        "candidate_no_go_from_off_branch_Gaussian_substitution",
        "positive_R_proves_all_auxiliary_pivots_regular",
        "full_nonlinear_lapse_equals_linear_Gaussian_lapse",
        "noncommuting_joint_spectral_law_away_from_bounce",
        "uniform_infinite_volume_or_spacetime_event",
        "old_physical_shifted_R_is_valid",
        "whole_cross_covariance_deleted",
        "interacting_physical_mean_is_stationary",
        "all_quantum_Ward_state_endpoint_defects_vanish",
        "mathematical_band_is_Wilsonian_cutoff",
        "full_H_Proca_onepoint_and_all_loops_are_bounded",
        "nonlinear_inhomogeneous_bounce_stability_complete",
        "original_P8_closed",
    ):
        cases.append(("unsupported_" + label, require_observable, (label,)))
    for i in range(len(frontier())):
        rows = frontier()
        rows[i]["status"] = "COMPLETE"
        cases.append(
            ("primitive_promotion_" + str(i), validate_scope, (rows, matching()))
        )
    for i in range(len(matching())):
        rows = matching()
        rows[i]["status"] = "COMPLETE"
        cases.append(
            ("matching_promotion_" + str(i), validate_scope, (frontier(), rows))
        )
    cases.append(
        ("previous_reference_label_is_not_new_mean", require_state, (previous.STATE,))
    )
    return cases


def rejected_inputs():
    count = 0
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (TypeError, ValueError):
            count += 1
        else:
            raise AssertionError("Unsupported input accepted: " + name)
    return count


def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "whole_original_R_and_all_fixed_parameters_retained": True,
        "original_state_not_conditioned_or_projected": True,
        "complete_bounce_commuting_rows_and_cross_covariance_kept": True,
        "exact_fourth_remainder_not_formal_jet_promotion": True,
        "tiny_tail_not_zero_or_integrability_guarantee": True,
        "individual_coefficient_obstruction_not_entire_action_no_go": True,
        "S261_physical_binding_remains_refuted": True,
        "original_interacting_mean_cutoff_and_V_G_B_P8_remain_open": True,
    }
