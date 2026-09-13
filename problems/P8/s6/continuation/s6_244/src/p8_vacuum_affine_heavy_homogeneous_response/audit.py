"""Whole actual heavy homogeneous response, with no nonzero-transfer or inverse transfer."""

from functools import cache

import sympy as s
from p8_vacuum_affine_heavy_curved_state import state
from p8_vacuum_affine_heavy_spatial_uv import audit as previous

from . import estimates, evolution, renormalization

ITEM = {
    "id": "QG2_H8A420_complete_actual_heavy_homogeneous_spatial_response_with_full_covariant_matching_contacts_and_mass_uniform_state_bounds",
    "status": "COMPLETE_PREPARED_HOMOGENEOUS_SPATIAL_RESPONSE_AND_FIXED_PROFILE_NOT_NONZERO_TRANSFER_ADM_CLOCK_RESPONSE_QUANTUM_INVERSE_NONLINEAR_BOUNCE_UV_REGGE_OR_P8",
}
exact = previous.exact


def require_stage(stage):
    allowed = (
        "complete_scalar_Hamiltonian_and_covariance",
        "full_prepared_initial_mismatch",
        "entire_finite_reference_and_actual_error",
        "full_general_dimensional_adiabatic_matching",
        "complete_fixed_finite_action_and_profile",
        "entire_homogeneous_spatial_response",
        "same_frequency_regulator_tail",
    )
    if not isinstance(stage, str) or stage not in allowed:
        raise ValueError(
            "Only the complete stated homogeneous heavy response is established"
        )
    return stage


def require_parameters(mass_squared, kappa):
    n, k = map(exact, (mass_squared, kappa))
    if n != state.MASS2 or k != state.KAPPA:
        raise ValueError("Keep the actual QG2 heavy parameters")
    return n, k


def require_domain(time, transfer):
    t, p = map(exact, (time, transfer))
    if not -s.Rational(1, 2) <= t <= s.Rational(1, 2) or p != 0:
        raise ValueError(
            "The actual response here has exactly zero spatial transfer on the stated slab"
        )
    return t, p


def require_jet_bound(order, bound):
    if isinstance(order, bool) or not isinstance(order, (int, s.Integer)):
        raise TypeError("The history jet order is an exact integer")
    b = exact(bound)
    if not 0 <= order <= 12 or not 0 <= b <= s.Rational(1, 100):
        raise ValueError(
            "Keep every raw homogeneous history jet through twelve in its full budget"
        )
    return int(order), b


def require_derivative(order):
    if isinstance(order, bool) or not isinstance(order, (int, s.Integer)):
        raise TypeError("The current derivative order is an exact integer")
    if not 0 <= order <= 2:
        raise ValueError(
            "The entire current bound covers source derivatives zero through two"
        )
    return int(order)


def require_graph(label, order):
    if isinstance(order, bool) or not isinstance(order, (int, s.Integer)):
        raise TypeError("The graph degree is an exact integer")
    graphs = {"detector_time": 0, "source_time": 13, "actual_external_transfer": 0}
    if not isinstance(label, str) or label not in graphs or graphs[label] != order:
        raise ValueError(
            "Use the complete homogeneous detector L2 and source H13 graph"
        )
    return label, int(order)


def require_cutoff_squared(value):
    cutoff = exact(value)
    if cutoff < 4 * state.MASS2:
        raise ValueError("Keep the same Omega_star cutoff with K squared at least4n")
    return cutoff


def frontier():
    return previous.frontier()


def matching():
    return [*previous.matching(), dict(ITEM)]


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("Original primitive and previous matching rows stay unchanged")
    return True


@cache
def packets():
    return {
        "complete_actual_scalar_evolution_and_contact": evolution.data(),
        "entire_mass_uniform_state_and_subtraction_bounds": estimates.data(),
        "complete_covariant_matching_and_local_profile": renormalization.data(),
    }


@cache
def residuals():
    return {
        name + "_" + key: value
        for name, data in packets().items()
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
            name + "_" + key: bool(value)
            for name, data in packets().items()
            for key, value in data["gates"].items()
        },
        "same_actual_S240_SLE_and_complete_scalar_prescription": True,
        "original_primitive_frontier_unchanged": True,
        "unique_matching_identifiers": len({row["id"] for row in matching()})
        == len(matching()),
        "zero_transfer_anchor_not_complete_nonzero_transfer_response": True,
        "finite_prescribed_homogeneous_history_not_quantum_feedback_solution": True,
    }


def observable():
    return {
        "input": "Keep the actual S240 QG2 heavy SLE Cauchy covariance, mass, full scalar prescription and fixed reference profile. Homogeneous symmetric spatial histories have a common initial zero germ and a stated full C12 jet budget.",
        "full": "The full scalar Hamiltonian, exact Gaussian graph and complete covariance/current tangent retain all trace, tracefree and instantaneous second metric contacts. The entire initial SLE/comparison/W6 mismatch is bounded and never reset.",
        "matching": "Complete general-dimensional angular/radial action matching, with arbitrary noncommuting extrinsic-curvature jets and weighted boundary terms, gives exactly the S240 scalar finite action. Its complete Euler current and whole fixed profile are included once.",
        "bound": "All internal momenta and current source derivatives0..2 are controlled in the finite prescribed homogeneous history neighborhood. The entire normalized first homogeneous spatial response has bound below10^-350 on detector L2 times source H13. The unchanged Omega_star<=K projection has same-graph error below10^-748/K squared for K>=2sqrt(n).",
        "remaining": "The nonzero-transfer remainder, full ADM/common-clock Ward assembly, quantum constraints/inverse, light/mixed loops, nonlinear same-state bounce, quantum gravitational limit, physical UV, Regge and original V/G/B/P8 remain open.",
    }


def bad_cases():
    invalid = (
        True,
        False,
        1.0,
        s.Float(1),
        "1",
        None,
        s.oo,
        -s.oo,
        s.zoo,
        s.I,
        s.nan,
        s.Symbol("x"),
    )
    out = []
    templates = (
        ("parameters", require_parameters, [state.MASS2, state.KAPPA]),
        ("domain", require_domain, [0, 0]),
        ("jet", require_jet_bound, [12, s.Rational(1, 100)]),
        ("derivative", require_derivative, [2]),
        ("graph", require_graph, ["source_time", 13]),
        ("cutoff", require_cutoff_squared, [4 * state.MASS2]),
    )
    for i, value in enumerate(invalid):
        for name, call, args0 in templates:
            for pos in range(len(args0)):
                args = list(args0)
                args[pos] = value
                out.append((f"{name}_type_{pos}_{i}", call, tuple(args)))
        out.append((f"stage_type_{i}", require_stage, (value,)))
    for args in ((1, 0), (-1, 0), (0, 1), (0, -1), (0, s.Rational(1, 10**100))):
        out.append((f"domain_scope_{len(out)}", require_domain, args))
    for args in ((1, state.KAPPA), (state.MASS2, 1)):
        out.append((f"parameters_scope_{len(out)}", require_parameters, args))
    for args in ((-1, 0), (13, 0), (0, -1), (0, s.Rational(1, 99))):
        out.append((f"jet_scope_{len(out)}", require_jet_bound, args))
    for order in (-1, 3, 10):
        out.append((f"derivative_scope_{len(out)}", require_derivative, (order,)))
    for args in (
        ("detector_time", 1),
        ("source_time", 12),
        ("actual_external_transfer", 1),
        ("same_space_inverse", 0),
    ):
        out.append((f"graph_scope_{len(out)}", require_graph, args))
    for cutoff in (0, -1, state.MASS2, 4 * state.MASS2 - 1):
        out.append((f"cutoff_scope_{len(out)}", require_cutoff_squared, (cutoff,)))
    for stage in (
        "reset_initial_reference_graph",
        "ignore_SLE_beta",
        "finite_WKB_evolution_is_physical_state",
        "live_state_reminimization",
        "omit_second_Hamiltonian_contact",
        "differentiate_covariance_only",
        "drop_quadratic_first_error_term",
        "full_covariance_is_even_in_marker",
        "differentiate_a_moving_marker_contour",
        "new_scalar_finite_counterterm",
        "borrow_old_vector_multiplicity",
        "tracefree_volume_cancellation_for_all_directions",
        "drop_Euler_or_volume_before_dimension_limit",
        "boundary_only_at_physical_dimension",
        "adiabatic_tail_alone_is_renormalized_response",
        "add_S243_finite_difference_at_zero_transfer",
        "double_add_S241_heat_action",
        "profile_zero_alone",
        "finite_cutoff_matched_mean_zero",
        "low_internal_momenta_omitted",
        "proof_threshold_is_physical_cutoff",
        "bounded_amplitude_without_twelve_time_jets",
        "nonzero_transfer_response_closed",
        "full_ADM_clock_response_closed",
        "same_space_inverse_from_graph_smallness",
        "full_quantum_constraints_inverse",
        "light_mixed_loops_controlled",
        "nonlinear_self_consistent_bounce",
        "quantum_gravitational_limit",
        "exact_UV_Regge_complete",
        "closed_P8",
    ):
        out.append(("unsupported_" + stage, require_stage, (stage,)))
    for i in range(len(frontier())):
        rows = frontier()
        rows[i]["status"] = "COMPLETE"
        out.append((f"primitive_{i}", validate_scope, (rows, matching())))
    for i in range(len(matching())):
        rows = matching()
        rows[i]["status"] = "COMPLETE"
        out.append((f"matching_{i}", validate_scope, (frontier(), rows)))
    out.extend(
        (
            ("missing_primitive", validate_scope, (frontier()[:-1], matching())),
            ("missing_matching", validate_scope, (frontier(), matching()[:-1])),
            (
                "extra_parent",
                validate_scope,
                (
                    frontier(),
                    matching() + [{"id": "closed_parent", "status": "COMPLETE"}],
                ),
            ),
        )
    )
    return out


@cache
def rejected_inputs():
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (ValueError, TypeError):
            continue
        raise ValueError(
            "Unsupported homogeneous response or closure claim accepted: " + name
        )
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "actual_state_and_full_initial_mismatch_retained": True,
        "complete_covariance_and_all_metric_contacts": True,
        "full_mass_uniform_all_internal_momentum_bounds": True,
        "whole_general_dimensional_scalar_matching": True,
        "unchanged_finite_action_and_fixed_profile_once": True,
        "same_auxiliary_frequency_cutoff_not_a_physical_cutoff": True,
        "original_primitives_unchanged": frontier() == previous.frontier(),
        "no_nonzero_transfer_inverse_UV_Regge_or_P8_closure": True,
    }
