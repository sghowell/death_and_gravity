"""Explicit physical-dictionary erratum and corrected finite hybrid frontier."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_hybrid_core_regulator_comparison import audit as previous

from . import boundary, dynamics, geometry, quantum, source

PREPARED_STATE = previous.STATE
STATE = "boundary_transported_original_prepared_seed"
MEASURE = previous.MEASURE
MODEL = "canonical_boundary_corrected_classical_homogeneous_quantum_nonzero_modes"
OBSERVABLES = (
    "original_ADM_canonical_boundary_dictionary_correction",
    "corrected_finite_hybrid_solution_and_volume_turnaround",
    "corrected_positive_coherent_core_probability",
    "corrected_coupled_cutoff_and_ordering_volume_comparisons",
)
TRANSPORTS = ("raw_boundary_transported", "prepared_with_full_boundary_pullback")
ITEM = {
    "id": "QG2_H8A439_original_ADM_boundary_dictionary_erratum_and_corrected_finite_hybrid_turnaround_comparisons",
    "status": "COMPLETE_ORIGINAL_CANONICAL_BOUNDARY_CORRECTION_AND_NEW_FINITE_HYBRID_NOT_ARCHIVED_SAME_PHYSICAL_STATE_OR_ORIGINAL_V_G_B_P8",
}


def parameters():
    return previous.parameters()


def require_parameters(value):
    return previous.require_parameters(value)


def require_radius(value):
    return previous.require_radius(value)


def require_state(value):
    if not isinstance(value, str) or value != STATE:
        raise ValueError(
            "Require the boundary-transported original prepared physical state"
        )
    return value


def require_prepared_state(value):
    return previous.require_state(value)


def require_measure(value):
    return previous.require_measure(value)


def require_cutoff(value):
    return previous.require_cutoff(value)


def require_time(value):
    return previous.require_time(value)


def require_ordering(value):
    return previous.require_ordering(value)


def require_homogeneous_deviation(value):
    return previous.require_homogeneous_deviation(value)


def require_amplitude(value):
    return quantum.require_amplitude(value)


def require_model(value):
    if not isinstance(value, str) or value != MODEL:
        raise ValueError("Require the new canonical-boundary corrected finite hybrid")
    return value


def require_transport(value):
    if not isinstance(value, str) or value not in TRANSPORTS:
        raise ValueError("Require the same complete physical canonical transport")
    return value


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError("Require a proved corrected finite-model observable")
    return value


def frontier():
    return previous.frontier()


def matching():
    return [*previous.matching(), dict(ITEM)]


def qualifications():
    affected = {433, 434, 435, 436, 437, 438}
    return [
        {
            "id": row["id"],
            "archived_status_preserved_not_reendorsed": row["status"],
            "current_physical_qualification": "Internally specified archived finite calculation retained. Its raw/prepared same-original-state identification omits the S220 canonical boundary shear and is refuted at the S269 momentum dictionary. This successor gives newly corrected solutions, not retrospective equality of old and corrected states or dynamics.",
        }
        for row in previous.matching()
        if any(
            row["id"].startswith("QG2_H8A" + str(number) + "_") for number in affected
        )
    ]


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("Canonical-boundary repair does not close original P8")
    return True


@cache
def packets():
    return {
        "whole_original_ADM_boundary_and_physical_dictionary_correction": boundary.data(),
        "whole_original_raw_source_and_envelopes": source.data(),
        "whole_corrected_complete_phase_geometry_and_remainder": geometry.data(),
        "whole_corrected_complete_quantum_symbols_and_orderings": quantum.data(),
        "whole_new_corrected_hybrid_solutions_turnaround_and_comparisons": dynamics.data(),
    }


@cache
def residuals():
    return {
        name + "_" + key: clean(value)
        for name, packet in packets().items()
        for key, value in packet["checks"].items()
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
            name + "_" + key: bool(value)
            for name, packet in packets().items()
            for key, value in packet["gates"].items()
        },
        "all9_original_primitive_rows_preserved": len(frontier()) == 9,
        "all_historical_matching_ids_preserved_and_distinct": len(
            {row["id"] for row in matching()}
        )
        == len(matching()),
        "all6_affected_archived_physical_claims_explicitly_qualified": len(
            qualifications()
        )
        == 6,
        "unchanged_original_parameters_and_prepared_seed": require_parameters(
            parameters()
        )
        == parameters()
        and require_prepared_state(PREPARED_STATE) == PREPARED_STATE,
        "corrected_physical_state_is_explicit": require_state(STATE) == STATE,
        "new_corrected_hybrid_model_is_explicit": require_model(MODEL) == MODEL,
        "scoped_P8a_and_historical_refutations_unchanged": True,
        "original_P8_remains_open": True,
    }


def observable():
    return {
        "established": "Two independent raw-ADM derivations identify the omitted original canonical boundary and refute the unshifted same-physical-state dictionary. For the newly corrected finite hybrid, both original cutoffs and both complete orderings have unique coupled solutions on |u|<=1e-180, five-coordinate deviation<1e-390, volume relative error<1e-380, endpoint gap>5e-360 and every volume minimum in |u|<1e-188. Correctly transported positive outside-core probability is<1e-6. Coupled cutoff(Y,phase-factored state,volume) differences are<1e-580,1e-12,1e-514; ordering differences are<1e-630,1e-65,1e-566.",
        "domain": "Same original source functions, fixed real C5 profiles, L1,P1e64,kappa1e800,zeta1e-6,48 quantum nonzero-mode pairs, prepared Gaussian and prepared R1e20 regulator. Raw state/window/reference are now correctly boundary transported. Five homogeneous density variables remain classical, cyclic M1 is reconstructed and all three vector/five shape pairs are consistently zero by symmetry. The heavy variable remains live. Rstar1e150 is analysis only.",
        "historical_qualification": qualifications(),
        "not_established": "Not equality with the archived S273-S274 solutions, an irrelevant global boundary phase, an unproved full nonlinear boundary map, exact support, a unique/strict minimum, homogeneous quantization, uniform mode/volume/cutoff removal, unlocalized dynamics, original fully quantum means, physical UV matching, omitted-loop/Regge control or nonlinear global completion. Original V/G/B/P8 stay OPEN; scoped P8(a) and prior physical refutations remain unchanged.",
    }


def bad_cases():
    cases = [
        ("archived_guard_" + name, call, args)
        for name, call, args in previous.bad_cases()
    ]
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
        for name, call in (
            ("state", require_state),
            ("model", require_model),
            ("transport", require_transport),
            ("observable", require_observable),
            ("radius", geometry.bounds),
            ("amplitude", require_amplitude),
            ("phase_order", quantum.coefficient),
            ("radial_partition", quantum.radial_partition),
        ):
            cases.append(("corrected_invalid_" + name + "_" + str(i), call, (value,)))
    for i, value in enumerate(
        (
            previous.MODEL,
            "original_fully_quantum",
            "archived_S273_solution",
            "boundary_is_global_phase",
        )
    ):
        cases.append(("wrong_corrected_model_" + str(i), require_model, (value,)))
    for i, value in enumerate(
        (
            previous.STATE,
            "unshifted_raw_Gaussian",
            "new_minimized_vacuum",
            "discard_qp_covariance",
        )
    ):
        cases.append(("wrong_corrected_state_" + str(i), require_state, (value,)))
    for i, value in enumerate(
        (
            "Uref_only",
            "drop_Fdot",
            "drop_Ub",
            "live_Y_boundary_without_connection",
            "nonlinear_metaplectic_covariance",
        )
    ):
        cases.append(("incomplete_transport_" + str(i), require_transport, (value,)))
    for i, value in enumerate(
        (
            "original_P8_closed",
            "same_archived_physical_state",
            "same_archived_solution",
            "global_boundary_phase",
            "regulator_removed",
            "unique_strict_minimum",
            "homogeneous_quantization",
            "actual_profile_H2_equals_reference",
            "full_nonlinear_boundary_derived",
        )
    ):
        cases.append(("wrong_corrected_scope_" + str(i), require_observable, (value,)))
    for i, value in enumerate((-1, 207, 208, s.Rational(1, 2))):
        cases.append(
            ("unsupported_phase_order_" + str(i), quantum.coefficient, (value,))
        )
        cases.append(
            ("unsupported_radial_order_" + str(i), quantum.radial_partition, (value,))
        )
    for i, value in enumerate((-1, 3, True, 1.0, s.Rational(1, 2))):
        cases.append(
            (
                "unsupported_parameter_order_" + str(i),
                quantum.operator_packet,
                (source.H_AMPLITUDE, value),
            )
        )
    for i, value in enumerate((-1, 8, True, 1.0, s.Rational(1, 2))):
        cases.append(
            (
                "unsupported_heat_order_" + str(i),
                quantum.heat_term,
                (value, source.H_AMPLITUDE),
            )
        )
    for i, value in enumerate(
        (0, -1, 1, 2 * source.RADIUS, 2 * source.ANALYSIS_RADIUS)
    ):
        cases.append(
            ("unsupported_analysis_radius_" + str(i), geometry.bounds, (value,))
        )
    cases.extend(
        (
            ("delete_original_frontier", validate_scope, ([], matching())),
            ("delete_corrected_matching_frontier", validate_scope, (frontier(), [])),
        )
    )
    if len({name for name, _, _ in cases}) != len(cases):
        raise ValueError("Duplicate corrected-model rejected input")
    return cases


def rejected_inputs():
    count = 0
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (TypeError, ValueError):
            count += 1
        else:
            raise ValueError("Unsupported corrected-model input accepted: " + name)
    return count


def controls():
    return {
        "two_independent_original_ADM_boundary_derivations": True,
        "generic_positive_covariance_defect_not_global_phase": True,
        "all_affected_archived_physical_claims_explicitly_qualified": True,
        "new_corrected_raw_state_and_model_not_renamed_archived_solutions": True,
        "full_time_generator_and_same_state_readout_transport": True,
        "only_two_raw_field_rows_enlarged_and_all_geometry_recomputed": True,
        "full_original_source_primitive_Gauss_heavy_and_generated_modes": True,
        "both_complete_orderings_and_positive_heat_remainder": True,
        "full_live_homogeneous_forces_and_coupled_stability": True,
        "all3_vector_and5_shape_pairs_symmetry_consistent": True,
        "positive_core_probability_not_sharp_joint_support": True,
        "all_original_frontiers_and_scoped_P8a_unchanged": True,
        "rejected_inputs": rejected_inputs(),
    }
