"""Strict original-parent domain and frontier contract for the nonlinear branch."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_finite_window_growth.intervals import fraction
from p8_vacuum_affine_nonlinear_reference_volume import audit as previous

from . import branch, response, source

STATE = "unchanged_original_preparations_no_Gaussian_invariant_assignment"
MEASURE = "classical_regular_auxiliary_branch_not_quantum_measure"
OBSERVABLES = (
    "entire_unique_nonlinear_lapse_and_normal_vector",
    "all_full_implicit_and_physical_volume_chain_contacts",
)
ITEM = {
    "id": "QG2_H8A430_full_quantitative_nonlinear_bounce_slice_auxiliary_branch",
    "status": "COMPLETE_LOCAL_FULL_NONLINEAR_LAPSE_NORMAL_VECTOR_DOMAIN_AND_IMPLICIT_RESPONSE_NOT_CANONICAL_GAUSSIAN_SUPPORT_INTERACTING_MEAN_QUANTUM_REGULATOR_OR_ORIGINAL_V_G_B_P8",
}


def parameters():
    return previous.parameters()


def require_parameters(values):
    return previous.require_parameters(values)


def require_time(value):
    value = fraction(value)
    if value != 0:
        raise ValueError(
            "The full nonlinear reconstruction is proved only on the bounce slice"
        )
    return value


def require_state(value):
    if not isinstance(value, str) or value != STATE:
        raise ValueError(
            "No original state is projected, reselected or assigned to density invariants"
        )
    return value


def require_measure(value):
    if not isinstance(value, str) or value != MEASURE:
        raise ValueError("The classical branch does not specify a quantum measure")
    return value


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError(
            "Require the whole nonlinear constraint or full physical chain"
        )
    return value


def require_coordinates(values):
    expected = {str(z) for z in source.COORDS}
    if not isinstance(values, dict) or set(values) != expected:
        raise ValueError(
            "Require every one of the twelve density and spatial jet coordinates"
        )
    result = {key: fraction(value) for key, value in values.items()}
    if any(abs(value) > fraction(source.DELTA) for value in result.values()):
        raise ValueError("The full nonlinear branch requires the stated invariant box")
    return result


def frontier():
    return previous.frontier()


def matching():
    return [*previous.matching(), dict(ITEM)]


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError(
            "Do not promote the local nonlinear branch to quantum or original P8 closure"
        )
    return True


@cache
def packets():
    return {
        "whole_original_bounce_constraint_with_every_spatial_channel": source.data(),
        "whole_quantitative_nonlinear_root_contraction_and_responses": branch.enclosure(),
        "whole_complete_two_auxiliary_regular_branch_and_coefficient_domain": branch.auxiliaries(),
        "whole_all_mixed_implicit_root_derivatives": response.jets(),
        "whole_original_physical_volume_and_canonical_chart_contacts": response.physical(),
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
        "original_quantum_domain_and_V_G_B_P8_remain_open": True,
    }


def observable():
    return {
        "established": "On u=0 and the complete twelve-invariant box of radius1e-250, the exact original lapse constraint has one root in|N-1|<1e-245, unique on the outer1e-6 strip. At the unchanged reference center its source-only displacement is below1e-375; no shift is assumed nonzero. The explicit fixed-point iteration contracts by1/30. The complete temporal root, both auxiliary pivots, all mixed root derivatives through order3 and every original physical-volume second-chain contact are retained.",
        "domain": "The variables are normalized densities and spatial invariants, with H=eta/1e100 and the actual heavy mass preserved. They are not twelve independent canonical modes. Physical data form a subset of the signed analytic box. Full scalar/spatial-gauge constraints, nonlinear canonical invariant maps and the quantum operator domain remain distinct.",
        "not_established": "There is no projection or new initial state, Gaussian invariant law, interacting quantum mean, positive Weyl ordering, regulated continuum/spatial Ward construction, Wilsonian cutoff or nonlinear inhomogeneous stability theorem. S265's Gaussian individual-coefficient divergence and S261's physical-binding refutation remain unchanged. Original V/G/B/P8 remain open; completed scoped P8(a) stays closed.",
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
    zero = {str(z): 0 for z in source.COORDS}
    for i, value in enumerate(invalid):
        for label, call in (
            ("time", require_time),
            ("state", require_state),
            ("measure", require_measure),
            ("observable", require_observable),
            ("coordinates", require_coordinates),
        ):
            cases.append(("invalid_" + label + "_" + str(i), call, (value,)))
        for key in zero:
            changed = {**zero, key: value}
            cases.append(
                (
                    "invalid_coordinate_" + key + "_" + str(i),
                    require_coordinates,
                    (changed,),
                )
            )
        for key in parameters():
            changed = parameters()
            changed[key] = value
            cases.append(
                (
                    "invalid_parameter_" + key + "_" + str(i),
                    require_parameters,
                    (changed,),
                )
            )
    for key in zero:
        for sign in (-1, 1):
            cases.append(
                (
                    "outside_" + key + "_" + str(sign),
                    require_coordinates,
                    ({**zero, key: sign * 2 * source.DELTA},),
                )
            )
        missing = {name: value for name, value in zero.items() if name != key}
        cases.append(("missing_" + key, require_coordinates, (missing,)))
    cases.append(("extra_coordinate", require_coordinates, ({**zero, "new_mode": 0},)))
    for value in (-source.fifth.TIME, source.fifth.TIME, 1):
        cases.append(("non_bounce_" + str(value), require_time, (value,)))
    for value in ("new_state", "conditioned_Gaussian", "Gaussian_density_invariants"):
        cases.append(("unsupported_state_" + value, require_state, (value,)))
    for value in (
        "flat_quantum_measure",
        "full_BRST_regulator",
        "positive_Weyl_symbol",
    ):
        cases.append(("unsupported_measure_" + value, require_measure, (value,)))
    for value in (
        "full_interacting_mean",
        "entire_action_divergence",
        "original_P8_closed",
    ):
        cases.append(("unsupported_observable_" + value, require_observable, (value,)))
    cases.append(("drop_original_frontier", validate_scope, ([], matching())))
    cases.append(("drop_original_matching", validate_scope, (frontier(), [])))
    assert len({name for name, _, _ in cases}) == len(cases)
    return cases


def rejected_inputs():
    count = 0
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (TypeError, ValueError):
            count += 1
        else:
            raise ValueError("Invalid input was accepted: " + name)
    return count


def controls():
    return {
        "whole_original_source_all_spatial_channels_and_mass_retained": True,
        "centered_exact_bound_not_independent_rounded_subtraction": True,
        "full_nonlinear_root_not_identified_with_linear_lapse": True,
        "all_canonical_invariant_chain_contacts_retained": True,
        "original_state_not_conditioned_or_assigned_to_invariants": True,
        "S265_individual_coefficient_obstruction_and_S261_refutation_unchanged": True,
        "remaining_spatial_gauge_and_quantum_domain_not_solved": True,
        "original_V_G_B_P8_still_open": True,
        "rejected_inputs": rejected_inputs(),
    }
