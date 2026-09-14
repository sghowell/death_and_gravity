"""Source-pinned correction of an unpublished frozen physical-map error."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_finite_window_growth.intervals import fraction
from p8_vacuum_affine_gauge_mean_transport import audit as previous

from . import actual, contact, source

STATE = "unchanged_full_reference_with_explicit_S261_physical_binding_refutation"
MEASURE = (
    "source_pinned_classical_ADM_and_Weyl_contacts_not_quantum_ordering_or_regulator"
)
OBSERVABLES = (
    "literal_current_metric_scalar_Maxwell_and_source_coefficient_bridge",
    "whole_actual_volume_clock_jets_and_source_contacts",
    "corrected_fixed_reference_physical_volume_cancellation",
    "explicit_S261_erratum_with_retained_generic_and_covariance_results",
)
ITEM = {
    "id": "QG2_H8A426_source_pinned_physical_conformal_volume_correction_and_explicit_S261_erratum",
    "status": "CORRECTED_FULL_PHYSICAL_ADM_VOLUME_BINDING_AND_FIXED_REFERENCE_CONTACTS_WITH_EXPLICIT_S261_REFUTATION_NOT_INTERACTING_MEAN_REGULATOR_ORDERING_CUTOFF_ORIGINAL_V_G_B_OR_P8",
}


def parameters():
    return previous.parameters()


def require_parameters(values):
    return previous.require_parameters(values)


def require_state(label):
    if not isinstance(label, str) or label != STATE:
        raise ValueError(
            "Keep the original reference and the explicit physical-binding correction"
        )
    return label


def require_measure(label):
    if not isinstance(label, str) or label != MEASURE:
        raise ValueError(
            "The corrected classical source map is not a completed quantum measure"
        )
    return label


def require_observable(label):
    if not isinstance(label, str) or label not in OBSERVABLES:
        raise ValueError(
            "Require the literal parent bridge or explicit corrected source identity"
        )
    return label


def require_bounce_slope(value):
    value = fraction(value)
    if value != fraction(s.Rational(3, 2)):
        raise ValueError(
            "The original physical spatial-volume lapse slope is exactly 3/2"
        )
    return value


def frontier():
    return previous.frontier()


def matching():
    return [*previous.matching(), dict(ITEM)]


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError(
            "Historical payloads stay frozen; the explicit erratum controls interpretation"
        )
    return True


@cache
def packets():
    return {
        "literal_physical_metric_binding_from_original_coefficients": source.coefficient_bridge(),
        "whole_actual_six_family_clock_jets": source.clock_jets(),
        "whole_physical_metric_inverse_density_and_Jacobian": actual.metric(),
        "whole_physical_scalar_and_Maxwell_ADM_bridge": actual.matter(),
        "whole_actual_volume_source_gradient_and_Hessian": contact.whole_volume_jets(),
        "whole_corrected_physical_mean_and_covariance_contacts": contact.corrected_contact(),
        "whole_explicit_erratum_and_retained_reference_symbols": contact.retained_symbols_and_erratum(),
    }


@cache
def residuals():
    return {
        packet + "_" + name: clean(value)
        for packet, data in packets().items()
        for name, value in data["checks"].items()
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
            packet + "_" + name: bool(value)
            for packet, data in packets().items()
            for name, value in data["gates"].items()
        },
        "all_nine_original_primitive_rows_unchanged": len(frontier()) == 9,
        "all_matching_identifiers_unique": len({row["id"] for row in matching()})
        == len(matching()),
        "all_original_parameters_unchanged": require_parameters(parameters())
        == parameters(),
        "S261_historical_payload_not_renewed_physical_certification": True,
        "original_P8_and_quantum_mean_still_open": True,
    }


def observable():
    return {
        "positive_result": "The actual physical metric factor is pinned independently by C=Cchi^2=M^-2=R**(-1/2). The entire ADM metric and inverse, positive densities, scalar/Maxwell/mass/source coefficients, six clock-jet families and original cotangent Jacobian agree. Correct physical volume onepoint and two-point contacts cancel with c1(0)=3/2.",
        "explicit_refutation": "The immutable S261 report incorrectly used R-1/2 and identified its -6 volume slope as the physical parent value. That identification is REFUTED, although the generic cancellation and replay remain algebraically reproducible. The original action was never changed. Its correct coefficient is +3/2, and its physical onepoint contact differs from the old claimed one by -135 Cvn/8.",
        "retained_scope": "S261 generic gauge orbits, unprojected Fourier displacement, unchanged coupled covariance/coordinate means, canonical commutator and actual bounce covariance symbols remain within their stated noninteracting/Weyl/asymptotic limits. No frozen source or report is edited; the explicit current erratum overrides the false physical-binding interpretation of the historical matching payload.",
        "original_problem": "Full quantum ordering, BRST regulator, Ward/state/endpoint defects, physical subtraction, interacting fixed mean, actual cutoff and matching, nonlinear bounce, quantum gravity limit, physical UV and finite-gravity IR/Regge remain OPEN. Original V/G/B/P8 and scoped P8(a) qualifications are unchanged.",
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
            ("slope", require_bounce_slope),
        ):
            cases.append(("invalid_" + label + "_" + str(i), call, (value,)))
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
    for value in (-6, 0, 1, 2, -s.Rational(3, 2), s.Rational(5, 2)):
        cases.append(
            ("false_actual_volume_slope_" + str(value), require_bounce_slope, (value,))
        )
    for label in (
        "R_minus_half_is_physical_metric_factor",
        "all_S261_physical_claims_certified",
        "generic_Ward_identity_selects_actual_parent",
        "all_replay_passes_imply_physical_binding",
        "source_or_report_bytes_repaired_after_freeze",
        "actual_interacting_Nielsen_vector_evaluated",
        "Weyl_ordering_is_original_nonlinear_quantum_ordering",
        "physical_source_mean_factorized",
        "Gaussian_support_inside_auxiliary_branch",
        "finite_Fourier_BRST_regulator_complete",
        "physical_volume_changes_under_spatial_coordinate_map",
        "fixed_state_reprepared",
        "cutoff_equals_heavy_mass",
        "ultraviolet_asymptotic_is_finite_loop_error",
        "original_P8_closed",
        "all_P8a_qualifications_removed",
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
            (
                "historical_matching_promotion_" + str(i),
                validate_scope,
                (frontier(), rows),
            )
        )
    cases.append(
        ("previous_state_without_explicit_erratum", require_state, (previous.STATE,))
    )
    return cases


@cache
def rejected_inputs():
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (TypeError, ValueError):
            continue
        raise ValueError("Unsupported corrected physical-map claim accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "actual_C_bound_to_independent_original_U_Cchi_M_coefficients": True,
        "frozen_S261_false_physical_identification_explicitly_refuted": True,
        "generic_cancellation_not_used_as_physical_binding_test": True,
        "whole_shift_metric_Maxwell_source_and_Jacobian_bridge_retained": True,
        "correct_onepoint_and_twopoint_contacts_not_deleted": True,
        "existing_state_and_valid_reference_symbols_not_reselected": True,
        "no_full_quantum_regulator_mean_cutoff_or_P8_promotion": True,
        "historical_payloads_unchanged_but_physical_interpretation_corrected": frontier()
        == previous.frontier()
        and matching()[:-1] == previous.matching(),
    }
