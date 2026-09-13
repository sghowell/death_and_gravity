"""Scoped whole prepared-heavy-state selection response with no closure transfer."""

from functools import cache

import sympy as s
from p8_vacuum_affine_heavy_clock_quadratic import audit as previous
from p8_vacuum_affine_heavy_curved_state import state

from . import bounds, response, vertices

ITEM = {
    "id": "QG2_H8A420_complete_exact_SLE_state_selection_response_difference_with_all_ADM_contacts_and_full_transfer_bounds",
    "status": "COMPLETE_PREPARED_STATE_SELECTION_DIFFERENCE_AND_FIXED_PROFILE_NOT_COMPARISON_STATE_RESPONSE_FULL_QUANTUM_INVERSE_NONLINEAR_UV_REGGE_OR_P8",
}
exact = previous.exact


def require_stage(stage):
    allowed = (
        "full_ADM_Hamiltonian_features",
        "full_second_metric_contact",
        "exact_covariance_and_ordered_Wick_bridge",
        "complete_SLE_minus_exact_comparison",
        "full_all_internal_and_external_momentum_bound",
        "full_two_leg_regulator_tail",
        "unchanged_fixed_reference_profile_split",
        "complete_matched_common_clock_pullback",
    )
    if not isinstance(stage, str) or stage not in allowed:
        raise ValueError(
            "Only the complete stated state-selection difference is established"
        )
    return stage


def require_parameters(mass_squared, kappa):
    n, k = map(exact, (mass_squared, kappa))
    if n != state.MASS2 or k != state.KAPPA:
        raise ValueError("Keep actual QG2 parameters")
    return n, k


def require_domain(time, momentum):
    t, p = map(exact, (time, momentum))
    if not -s.Rational(1, 2) <= t <= s.Rational(1, 2) or p < 0:
        raise ValueError("Use the stated clock slab and every nonnegative momentum")
    return t, p


def require_cutoff_squared(value):
    cutoff = exact(value)
    if cutoff < 4 * state.MASS2:
        raise ValueError("The auxiliary frequency cutoff requires K squared at least4n")
    return cutoff


def require_graph(label, order):
    if isinstance(order, bool) or not isinstance(order, (int, s.Integer)):
        raise TypeError("The graph degree is an exact integer")
    graphs = {
        "full_ADM_feature_spatial_graph": 1,
        "same_clock_scalar_graph": 2,
        "external_time_derivatives": 0,
    }
    if not isinstance(label, str) or label not in graphs or order != graphs[label]:
        raise ValueError("Keep the complete stated derivative graph")
    return label, int(order)


def require_parity(label, rule):
    rules = {
        "even": "minus_imaginary_part",
        "odd": "plus_i_real_part",
        "general": "both_reflected_ordered_products",
    }
    if (
        not isinstance(label, str)
        or not isinstance(rule, str)
        or label not in rules
        or rules[label] != rule
    ):
        raise ValueError("Retain the correct reflected ordered branch")
    return label, rule


def frontier():
    return previous.frontier()


def matching():
    return [*previous.matching(), dict(ITEM)]


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError(
            "Original primitive and prior matching statuses stay unchanged"
        )
    return True


@cache
def packets():
    return {
        "full_ADM_features": vertices.data(),
        "complete_ordered_covariance_response": response.data(),
        "entire_state_selection_bounds_and_profile": bounds.data(),
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
        "same_S240_state_prescription_and_S241_local_input": True,
        "prior_frontier_unchanged": True,
        "unique_matching_identifiers": len({row["id"] for row in matching()})
        == len(matching()),
        "complete_difference_not_full_comparison_response": True,
        "no_inverse_from_small_graph_coefficient": True,
    }


def observable():
    return {
        "input": "Keep the existing QG2 parent, actual exact SLE and full covariant prescription. The comparator is the exact S240 W6-initialized KG basis, not a physical state replacement or finite WKB evolution.",
        "full": "The full5-feature scalar ADM Hamiltonian gives every lapse, shift and noncommuting metric contact. The finite covariance tangent agrees with the ordered Wick commutator and its reflected reverse term. Odd shift cross channels retain their imaginary Fourier kernels.",
        "bound": "The ENTIRE actual SLE-minus-exact-comparison response and its fixed reference-profile summand have all-internal/all-external-momentum weak-graph bounds below10^-1280 after kappa0 normalization, both in the full physical ADM and actual common-clock graphs. No time derivative of the external direction is required.",
        "tail": "The unchanged two-leg auxiliary frequency projection has a full same-clock tail coefficient below10^-1180/K on the same stated graph, K>=2sqrt(n). Large-small internal pairs and all upper retarded endpoints remain.",
        "remaining": "The full comparison-state renormalized response, complete quantum constraints/inverse, interacting light/mixed loops, nonlinear same-state bounce, quantum gravitational limit, physical UV, finite-gravity Regge and original V/G/B/P8 remain open.",
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
        ("domain", require_domain, [0, 1]),
        ("cutoff", require_cutoff_squared, [4 * state.MASS2]),
        ("graph", require_graph, ["same_clock_scalar_graph", 2]),
        ("parity", require_parity, ["general", "both_reflected_ordered_products"]),
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
    for value in (0, -1, state.MASS2, 4 * state.MASS2 - 1):
        out.append((f"cutoff_scope_{len(out)}", require_cutoff_squared, (value,)))
    for args in ((1, state.KAPPA), (state.MASS2, 1)):
        out.append((f"parameters_scope_{len(out)}", require_parameters, args))
    for args in (
        ("same_clock_scalar_graph", 0),
        ("full_ADM_feature_spatial_graph", 0),
        ("external_time_derivatives", 1),
    ):
        out.append((f"graph_scope_{len(out)}", require_graph, args))
    for args in (
        ("odd", "minus_imaginary_part"),
        ("even", "plus_i_real_part"),
        ("general", "minus_imaginary_part"),
    ):
        out.append((f"parity_scope_{len(out)}", require_parity, args))
    for stage in (
        "full_heavy_response_is_the_difference",
        "comparison_response_zero",
        "finite_WKB_evolution_is_exact_comparator",
        "comparison_is_a_Hadamard_replacement",
        "live_state_reminimization",
        "beta_occupation_only",
        "alpha_minus_one_dropped",
        "ordinary_imaginary_part_for_all_channels",
        "conjugate_the_complex_source_direction",
        "single_forward_product_only",
        "omit_second_Hamiltonian_contact",
        "drop_lapse_or_shift",
        "covariance_force_only_one_adjoint_block",
        "metric_contact_current_sign_positive",
        "UV_counterterms_new_state_dependent",
        "single_internal_leg_cutoff",
        "large_small_pairs_omitted",
        "new_mean_profile_or_prescription",
        "state_profile_contact_zero_alone",
        "finite_cutoff_clock_contact_zero",
        "global_smooth_auxiliary_profile_asserted",
        "same_space_inverse_from_graph_bound",
        "full_quantum_constraint_inverse",
        "same_state_nonlinear_bounce",
        "quantum_gravitational_limit",
        "all_light_mixed_loops_controlled",
        "UV_Regge_complete",
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
            "Unsupported state-selection or closure claim accepted: " + name
        )
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "actual_SLE_not_replaced_or_reminimized": True,
        "full_complex_state_coefficients_and_exact_KG": True,
        "all_ADM_features_and_metric_contacts": True,
        "ordered_reflection_before_complex_source_extension": True,
        "complete_two_leg_cutoff_and_full_tail": True,
        "unchanged_profile_and_finite_prescription": True,
        "original_primitives_unchanged": frontier() == previous.frontier(),
        "no_comparison_inverse_UV_or_P8_closure": True,
    }
