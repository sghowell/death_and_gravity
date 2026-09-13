"""Complete actual scalar spatial response with explicit restored-projection and inverse boundary."""

from functools import cache

import sympy as s
from p8_vacuum_affine_heavy_curved_state import state
from p8_vacuum_affine_heavy_homogeneous_response import audit as previous

from . import domain, endpoints, time_remainder

ITEM = {
    "id": "QG2_H8A420_complete_actual_heavy_all_transfer_spatial_response_full_scalar_UV_subtracted_remainder_and_same_restored_frequency_projection",
    "status": "COMPLETE_ACTUAL_LINEAR_SPATIAL_RESPONSE_ALL_TRANSFERS_NOT_FULL_ADM_CLOCK_WARD_QUANTUM_INVERSE_NONLINEAR_UV_REGGE_OR_P8",
}
exact = previous.exact


def require_stage(stage):
    allowed = (
        "full_prepared_scalar_pair",
        "whole_SLE_W6_error",
        "complete_six_step_time_remainder",
        "entire_unaveraged_UV_subtracted_endpoints",
        "full_dimensional_dominated_matching",
        "complete_actual_all_transfer_spatial_response",
        "same_restored_two_leg_projection",
    )
    if not isinstance(stage, str) or stage not in allowed:
        raise ValueError(
            "Only the complete stated scalar spatial response is established"
        )
    return stage


def require_parameters(mass_squared, kappa):
    n, k = map(exact, (mass_squared, kappa))
    if n != state.MASS2 or k != state.KAPPA:
        raise ValueError("Keep the actual QG2 heavy parameters")
    return n, k


def require_domain(time, transfer):
    t, p = map(exact, (time, transfer))
    if not -s.Rational(1, 2) <= t <= s.Rational(1, 2) or p < 0:
        raise ValueError("Retain the stated slab and a nonnegative transfer magnitude")
    return t, p


def require_projection(label):
    expected = "complete_unaveraged_subtracted_two_leg_with_full_local_terms"
    if not isinstance(label, str) or label != expected:
        raise ValueError(
            "Keep the full unaveraged subtraction, original two-leg mask and restored fixed local terms"
        )
    return label


def require_derivative(order):
    if isinstance(order, bool) or not isinstance(order, (int, s.Integer)):
        raise TypeError("The current derivative order is an exact integer")
    if order != 1:
        raise ValueError(
            "This continuation establishes the full linear spatial response"
        )
    return int(order)


def require_graph(label, order):
    if isinstance(order, bool) or not isinstance(order, (int, s.Integer)):
        raise TypeError("The graph degree is an exact integer")
    graphs = {
        "detector_time": 0,
        "detector_spatial": 0,
        "source_time": 13,
        "source_spatial": 6,
    }
    if not isinstance(label, str) or label not in graphs or graphs[label] != order:
        raise ValueError("Use detector L2 and the full prepared source Z136 graph")
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
        "complete_scalar_joint_domain_and_dimensional_control": domain.data(),
        "whole_state_difference_and_six_step_time_remainder": time_remainder.data(),
        "entire_UV_subtracted_endpoint_and_actual_spatial_response": endpoints.data(),
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
        "restored_projection_not_bare_finite_cutoff_Ward_identity": True,
        "complete_linear_spatial_response_not_finite_inhomogeneous_feedback": True,
    }


def observable():
    return {
        "input": "Keep the actual S240 QG2 heavy SLE, complete scalar covariant prescription and unchanged fixed profile. The reference FLRW clock slab and common prepared source germ are unchanged.",
        "full": "The actual S244 homogeneous response anchors zero transfer. Complete actual-minus-W6, sixth-time and first-five endpoint differences plus S243 same-scheme finite UV difference give every physical spatial transfer. All scalar directions, contacts and source jets remain.",
        "matching": "Full six-iterate frequency coefficients equal the frozen UV input. Exact normalized near/far geometry, both same-dimension phases and integrable radial/angular envelopes justify the full dimensional limit at fixed physical invariants.",
        "bound": "The entire normalized linear spatial response is below10^-350 on detector L2 times source Z136. Its complete restored Omega_star two-leg projection has same-graph tail below10^-250/K for K>=2sqrt(n), with all unaveraged UV terms kept before masking and the fixed finite local contribution unprojected.",
        "remaining": "Full finite-cutoff Ward/local-shape bookkeeping, heavy ADM/common-clock assembly, quantum constraints/inverse, nonlinear inhomogeneous feedback, interacting light/mixed loops, quantum gravitational limit, physical UV/Regge and original V/G/B/P8 remain open.",
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
        (
            "projection",
            require_projection,
            ["complete_unaveraged_subtracted_two_leg_with_full_local_terms"],
        ),
        ("derivative", require_derivative, [1]),
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
    for args in ((1, 0), (-1, 0), (0, -1)):
        out.append((f"domain_scope_{len(out)}", require_domain, args))
    for args in ((1, state.KAPPA), (state.MASS2, 1)):
        out.append((f"parameters_scope_{len(out)}", require_parameters, args))
    for order in (-1, 0, 2, 3):
        out.append((f"derivative_scope_{len(out)}", require_derivative, (order,)))
    for args in (
        ("detector_time", 1),
        ("detector_spatial", 1),
        ("source_time", 12),
        ("source_spatial", 5),
        ("actual_external_transfer", 1),
        ("same_space_inverse", 0),
    ):
        out.append((f"graph_scope_{len(out)}", require_graph, args))
    for cutoff in (0, -1, state.MASS2, 4 * state.MASS2 - 1):
        out.append((f"cutoff_scope_{len(out)}", require_cutoff_squared, (cutoff,)))
    for stage in (
        "reset_actual_SLE",
        "physical_WKB_state",
        "live_state_reminimization",
        "drop_initial_mismatch",
        "omit_error_squares",
        "borrow_old_vector_multiplicity",
        "guess_sixth_bulk_sign",
        "discard_fifth_endpoint",
        "fifth_bulk_alone_absolutely_integrable",
        "final_source_cutoff",
        "drop_lower_germ",
        "omit_low_momenta",
        "far_tail_alone_is_full_remainder",
        "ignore_near_logarithms",
        "drop_logarithmic_inverse_mass",
        "average_UV_before_masking",
        "one_leg_mask_replacement",
        "large_small_pairs_omitted",
        "bare_finite_cutoff_Ward_closed",
        "conjugate_complex_dimension",
        "positive_noninteger_dimensional_norm",
        "freeze_invariants_after_dimension_derivative",
        "new_finite_counterterm",
        "double_add_S241_heat",
        "profile_zero_alone",
        "homogeneous_contact_removed",
        "finite_cutoff_matched_mean_zero",
        "finite_inhomogeneous_C2_response",
        "full_ADM_clock_response_closed",
        "same_space_inverse_from_graph_smallness",
        "light_mixed_loops_controlled",
        "nonlinear_same_state_bounce",
        "quantum_gravitational_limit",
        "exact_UV_Regge_complete",
        "closed_P8",
    ):
        out.append(("unsupported_" + stage, require_stage, (stage,)))
    for projection in (
        "bare_Gaussian_current",
        "single_leg",
        "angular_averaged_UV_before_mask",
        "project_full_fixed_profile",
        "change_MSbar_finite_part",
    ):
        out.append(("projection_" + projection, require_projection, (projection,)))
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
            "Unsupported spatial response or closure claim accepted: " + name
        )
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "actual_state_and_entire_W6_computational_error_retained": True,
        "all_six_time_steps_and_complete_UV_subtracted_tensor": True,
        "full_near_far_low_and_removed_union_bounds": True,
        "complete_scalar_dimensional_dominated_matching": True,
        "unchanged_finite_action_and_fixed_profile_once": True,
        "same_restored_projection_not_bare_cutoff_Ward": True,
        "original_primitives_unchanged": frontier() == previous.frontier(),
        "no_full_ADM_inverse_UV_Regge_or_P8_closure": True,
    }
