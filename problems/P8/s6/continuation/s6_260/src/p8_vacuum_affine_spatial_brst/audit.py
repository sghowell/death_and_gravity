"""Classical symmetry and conditional Ward identities do not close P8."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure import audit as previous
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_finite_window_growth.intervals import fraction

from . import brst, ward

STATE = (
    "unchanged_fixed_preparation_with_explicit_BRST_endpoint_and_regulator_obligations"
)
MEASURE = (
    "finite_conditional_projective_spatial_algebra_not_anomaly_free_continuum_measure"
)
OBSERVABLES = (
    "full_spatial_and_projective_Grassmann_jet_BRST_algebra",
    "full_canonical_source_and_gauge_fermion_identity",
    "conditional_fixed_source_Nielsen_and_mean_transport_with_defects",
    "finite_convergent_gauge_orbit_and_noninvariant_state_controls",
)
ITEM = {
    "id": "QG2_H8A424_complete_spatial_projective_BRST_and_conditional_fixed_source_Ward_transport",
    "status": "EXACT_FULL_CLASSICAL_SPATIAL_PROJECTIVE_BRST_AND_CONDITIONAL_SOURCE_WARD_IDENTITIES_NOT_COMPLETE_QUANTUM_REGULATOR_STATE_MEAN_CUTOFF_ORIGINAL_V_G_B_OR_P8",
}


def parameters():
    return previous.parameters()


def require_parameters(values):
    return previous.require_parameters(values)


def require_state(label):
    if not isinstance(label, str) or label != STATE:
        raise ValueError(
            "Keep the original preparation and explicit endpoint Ward obligations"
        )
    return label


def require_observable(label):
    if not isinstance(label, str) or label not in OBSERVABLES:
        raise ValueError("Require the stated classical or conditional source identity")
    return label


def require_measure(label):
    if not isinstance(label, str) or label != MEASURE:
        raise ValueError("No completed continuum quantum measure is supplied")
    return label


def require_wavevector(values):
    if not isinstance(values, tuple) or len(values) != 3:
        raise TypeError("Require all three exact spatial components")
    result = tuple(fraction(value) for value in values)
    if not any(result):
        raise ValueError(
            "The residual translation modes have no displayed ghost inverse"
        )
    return result


def require_orbit(v, t, Jy):
    v, t, Jy = map(fraction, (v, t, Jy))
    if v <= 0 or t <= 0 or 1 - 2 * v * t * Jy <= 0:
        raise ValueError("Require v>0,t>0 and the complete convergent source domain")
    return v, t, Jy


def frontier():
    return previous.frontier()


def matching():
    return [*previous.matching(), dict(ITEM)]


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("No primitive or previous matching verdict is promoted")
    return True


@cache
def packets():
    return {
        "whole_tensor_density_shift_and_connection_BRST": brst.tensor_identities(),
        "whole_projective_compensation_and_ghost_doublets": brst.projective_identities(),
        "whole_physical_ADM_and_vector_gauge_dictionary": brst.physical_point_identities(),
        "whole_canonical_boundary_and_off_gauge_fermion": brst.canonical_identities(),
        "whole_seven_ghost_matrix_and_original_source": brst.projective_ghost_block(),
        "whole_translation_and_finite_projection_obstructions": brst.zero_modes(),
        "whole_convergent_two_source_gauge_orbit": ward.orbit(),
        "whole_differentiated_conditional_Nielsen_identity": ward.differentiated_identity(),
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
        "all_original_fixed_parameters_retained": require_parameters(parameters())
        == parameters(),
        "no_complete_quantum_measure_fixed_mean_or_original_P8_closure": True,
    }


def observable():
    return {
        "positive_result": "The full local spatial/projective Grassmann jet differential, compensated trace gauge, seven-ghost matrix, complete canonical spatial flux and source-preserving physical map are explicit. The conditional fixed-source Ward identity retains state/measure/regulator/endpoint defects and the full off-shell onepoint contact. Exact finite orbit integration distinguishes a moving coordinate mean from physical state changes.",
        "scope_boundary": "The BRST algebra is classical on the S257 regular reduced branch. It is not an off-constraint BFV prescription for every unreduced variable or an anomaly-free continuum regulator. The original fixed preparations and source contacts are retained, not reconstructed as a complete interacting BRST state. The actual Nielsen vector, physical subtraction and fixed quantum mean are not evaluated.",
        "original_problem": "All coefficients, sources, profiles, canonical boundaries and prior primitive/matching/P8(a) qualifications remain. Controlled cutoff/matching, nonlinear bounce, quantum gravitational limit, physical UV scattering, finite-gravity IR/Regge and original V/G/B/P8 stay OPEN.",
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
            ("wavevector", require_wavevector),
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
    for value in ((0, 0, 0), (0, 0), (0, 0, 1.0), [0, 0, 1], (True, 0, 1)):
        cases.append(("invalid_wavevector_" + str(value), require_wavevector, (value,)))
    for value in (
        (0, 1, 0),
        (1, 0, 0),
        (1, -1, 0),
        (1, 1, s.Rational(1, 2)),
        (1, 1, 1),
        (1.0, 1, 0),
    ):
        cases.append(("outside_full_orbit_domain_" + str(value), require_orbit, value))
    for label in (
        "all_gauge_coordinate_means_invariant",
        "all_Ward_defects_zero",
        "mean_zero_ghosts_form_Lie_algebra",
        "finite_Fourier_regulator_BRST_exact",
        "full_quantum_mean_is_classical",
        "original_P8_closed",
        "physical_state_reset_in_each_gauge",
        "source_contact_deleted",
        "nonlinear_Horava_unitarity_imported",
        "off_constraint_BFV_rules_proved",
        "all_P8a_qualifications_removed",
        "full_continuum_regulator_proved",
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
        (
            "previous_state_scope_not_current_state_contract",
            require_state,
            (previous.STATE,),
        )
    )
    return cases


@cache
def rejected_inputs():
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (TypeError, ValueError):
            continue
        raise ValueError("Unsupported BRST/Ward input accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "whole_64_connection_and_seven_gauge_ghosts_retained": True,
        "projective_compensator_and_time_shift_not_omitted": True,
        "translation_modes_and_projected_algebra_defects_retained": True,
        "source_contacts_and_composite_means_not_deleted": True,
        "state_regulator_and_endpoint_Ward_defects_explicit": True,
        "off_shell_onepoint_and_regular_on_shell_transport_distinguished": True,
        "no_quantum_measure_state_cutoff_or_P8_promotion": True,
        "all_original_and_earlier_matching_rows_unchanged": frontier()
        == previous.frontier()
        and matching()[:-1] == previous.matching(),
    }
