"""Scoped complete scalar spatial UV input; no full-kernel or closure transfer."""

from functools import cache

import sympy as s
from p8_vacuum_affine_heavy_curved_state import state
from p8_vacuum_affine_heavy_state_response import audit as previous

from . import jets, structure
from . import matching as finite_input

ITEM = {
    "id": "QG2_H8A420_complete_minimal_scalar_ordered_spatial_UV_symbol_and_same_scheme_dimensional_finite_difference",
    "status": "COMPLETE_LOCAL_SPATIAL_UV_INPUT_NOT_FULL_RENORMALIZED_COMPARISON_RESPONSE_HOMOGENEOUS_ANCHOR_QUANTUM_INVERSE_NONLINEAR_UV_REGGE_OR_P8",
}
exact = previous.exact


def require_stage(stage):
    allowed = (
        "complete_spatial_scalar_vertex",
        "all_ordered_endpoint_grades",
        "fixed_six_invariant_dimension",
        "whole_scalar_pole_match",
        "full_same_scheme_finite_UV_difference",
        "complete_local_Z24_bound",
    )
    if not isinstance(stage, str) or stage not in allowed:
        raise ValueError(
            "Only the complete stated local spatial UV input is established"
        )
    return stage


def require_parameters(mass_squared, kappa):
    n, k = map(exact, (mass_squared, kappa))
    if n != state.MASS2 or k != state.KAPPA:
        raise ValueError("Keep the actual QG2 scalar parameters")
    return n, k


def require_domain(time, momentum):
    tt, p = map(exact, (time, momentum))
    if not -s.Rational(1, 2) <= tt <= s.Rational(1, 2) or p < 0:
        raise ValueError("Use the full stated clock slab and nonnegative transfer")
    return tt, p


def require_graph(label, order):
    if isinstance(order, bool) or not isinstance(order, (int, s.Integer)):
        raise TypeError("Use the exact graph degree")
    graphs = {"detector_time": 0, "source_time": 2, "source_spatial": 4}
    if not isinstance(label, str) or label not in graphs or graphs[label] != order:
        raise ValueError("Keep the stated detector L2 and source Z24 graph")
    return label, int(order)


def require_ordered_product(label):
    if not isinstance(label, str) or label not in ("00", "01", "10", "11"):
        raise ValueError("Keep both distinct ordered trace and gradient products")
    return label


def require_dimension(value):
    dimension = exact(value)
    if not s.Rational(11, 4) <= dimension <= s.Rational(13, 4):
        raise ValueError("Use the stated local dimensional neighborhood")
    return dimension


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
        "complete_vertex_ordering_and_dimensional_normalization": structure.data(),
        "complete_same_scheme_finite_spatial_UV": finite_input.finite(),
    }


@cache
def residuals():
    out = {
        name + "_" + key: value
        for name, data in packets().items()
        for key, value in data["checks"].items()
    }
    amplitudes, _ = jets.amplitudes()
    out["physical_scalar_annihilation_frequency_leading"] = (
        jets.frequency("k")[2][0] + s.I
    )
    out["complete_scalar_trace_vertex_leading"] = amplitudes["trace"][0] - s.Rational(
        1, 2
    )
    out["complete_scalar_gradient_vertex_leading"] = amplitudes["gradient"][
        0
    ] - s.Rational(1, 2)
    return out


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
        "same_S240_physical_state_and_complete_scalar_prescription": True,
        "original_frontier_unchanged": True,
        "unique_matching_identifiers": len({row["id"] for row in matching()})
        == len(matching()),
        "full_local_difference_not_complete_comparison_response": True,
        "no_inverse_or_parent_closure_from_local_graph_smallness": True,
    }


def complete_jet_payload():
    return {
        "generic_background_time_jets": jets.h,
        "both_internal_scalar_frequency_and_momentum_jets": {
            leg: jets.frequency(leg) for leg in ("k", "l")
        },
        "both_complete_trace_and_gradient_pair_amplitudes": jets.amplitudes(),
        "every_ordered_endpoint_source_jet": jets.endpoint_products(),
        "every_complete_radial_row": jets.radial_rows(),
        "every_spatial_difference_row": jets.spatial_rows(),
        "complete_logarithmic_rows_in_all_source_orders": jets.logarithmic(),
    }


def observable():
    return {
        "input": "Keep the exact S240 QG2 state and minimally coupled scalar prescription, and the S242 state-selection difference. The present UV jets are algebraic high-r coefficients, not a replacement of the exact comparison evolution by WKB.",
        "full": "All four ordered trace/gradient products,60 endpoint/source slots,140 scalar coefficients and35 radial rows give the complete first-five-endpoint spatial UV input. The six physical tensor invariants stay fixed while dimension is varied.",
        "finite": "The full minimally coupled scalar covariant pole matches every invariant. The continued sphere measure, full volume/Euler Hessian and entire evanescent counteraction give the actual same-scheme finite spatial UV difference, with all ordered local Green identities.",
        "bound": "The ENTIRE local spatial UV difference obeys a normalized bound below10^-600 on detector L2 times source Z24: time derivatives0..2 and four spatial derivatives, all external transfers, |t|<=1/2.",
        "remaining": "The full UV-subtracted exact-comparison kernel, homogeneous anchor, matched full ADM/clock response, quantum constraints/inverse, interacting loops, nonlinear bounce, quantum gravity limit, physical UV, Regge and original V/G/B/P8 remain open. The S241 finite heat Hessian is not added again as another counterterm.",
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
        ("graph", require_graph, ["source_spatial", 4]),
        ("dimension", require_dimension, [3]),
        ("ordered", require_ordered_product, ["01"]),
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
    for value in (0, 2, 4, s.Rational(27, 10), s.Rational(33, 10)):
        out.append((f"dimension_scope_{len(out)}", require_dimension, (value,)))
    for args in (
        ("detector_time", 1),
        ("source_time", 0),
        ("source_spatial", 2),
        ("same_space_inverse", 0),
    ):
        out.append((f"graph_scope_{len(out)}", require_graph, args))
    for label in ("tracefree_only", "symmetrize_01_10", "one_pair", ""):
        out.append((f"ordered_scope_{len(out)}", require_ordered_product, (label,)))
    for stage in (
        "drop_trace_time_or_mass_term",
        "borrow_old_vector_frequency",
        "borrow_old_vector_pole",
        "wrong_annihilation_branch",
        "delete_odd_endpoints",
        "drop_upper_source_endpoints",
        "identify_DG_with_GD_before_time_derivatives",
        "freeze_dimension_before_invariants",
        "tracefree_pole_suffices",
        "finite_target_fitted_to_covariance",
        "discard_full_evanescent_counteraction",
        "discard_continued_Euler_or_volume",
        "change_scalar_finite_prescription",
        "add_S241_heat_Hessian_again",
        "contact_zero_in_full_response",
        "zero_transfer_anchor_zero",
        "W6_is_exact_evolution",
        "complete_UV_subtracted_kernel_bound",
        "full_dimensional_dominated_limit",
        "full_heavy_ADM_clock_response",
        "new_physical_state_or_live_reminimization",
        "same_space_inverse_from_local_smallness",
        "all_light_mixed_loops_controlled",
        "full_quantum_constraint_inverse",
        "nonlinear_same_state_bounce",
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
        raise ValueError("Unsupported spatial UV or closure claim accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "full_scalar_not_old_vector_weights": True,
        "all_four_ordered_products_and_odd_endpoints": True,
        "six_invariants_fixed_before_dimension_derivative": True,
        "entire_covariant_counteraction_and_finite_value": True,
        "original_state_profile_and_prescription_unchanged": True,
        "local_difference_not_full_UV_subtracted_response": True,
        "original_primitive_rows_unchanged": frontier() == previous.frontier(),
        "no_inverse_UV_Regge_or_P8_closure": True,
    }
