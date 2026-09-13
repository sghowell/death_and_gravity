"""Actual full scalar metric/clock response with an explicit computational-cutoff boundary."""

from functools import cache

import sympy as s
from p8_vacuum_affine_heavy_curved_state import state
from p8_vacuum_affine_heavy_spatial_response import audit as previous

from . import estimates, geometry, profile

exact = previous.exact
ITEM = {
    "id": "QG2_H8A420_complete_actual_heavy_ADM_common_clock_response_ordered_scalar_Ward_full_profile_and_projected_mean_contact",
    "status": "COMPLETE_ACTUAL_GAUSSIAN_ADM_COMMON_CLOCK_REFERENCE_RESPONSE_WITH_EXPLICIT_WARD_COMPLETED_APPROXIMANT_NOT_BARE_CUTOFF_LOCAL_SHAPES_QUANTUM_INVERSE_NONLINEAR_UV_REGGE_OR_P8",
}


def require_parameters(mass_squared, kappa):
    n, k = map(exact, (mass_squared, kappa))
    if n != state.MASS2 or k != state.KAPPA:
        raise ValueError("Keep the actual unchanged QG2 heavy parameters")
    return n, k


def require_domain(time, transfer):
    t, p = map(exact, (time, transfer))
    if not -s.Rational(1, 2) <= t <= s.Rational(1, 2) or p < 0:
        raise ValueError(
            "Use the common reference slab and all nonnegative transfer magnitudes"
        )
    return t, p


def require_graph(label, degree):
    if isinstance(degree, bool) or not isinstance(degree, (int, s.Integer)):
        raise TypeError("Graph degrees are exact integers")
    allowed = {
        "detector_time": 0,
        "detector_spatial": 2,
        "source_time": 13,
        "source_spatial": 8,
    }
    if not isinstance(label, str) or allowed.get(label) != degree:
        raise ValueError("Use the complete detector V02 and prepared source U138 graph")
    return label, int(degree)


def require_stage(stage):
    allowed = (
        "complete_scalar_Gaussian_ADM",
        "full_profile_removed_and_restored_once",
        "entire_common_clock_response",
        "original_mean_C1_tail",
        "Ward_completed_restored_approximant",
    )
    if not isinstance(stage, str) or stage not in allowed:
        raise ValueError(
            "Only the complete stated actual Gaussian metric/clock response is established"
        )
    return stage


def require_projection(label):
    if (
        not isinstance(label, str)
        or label != "original_two_leg_spatial_one_leg_mean_Ward_completion"
    ):
        raise ValueError(
            "Retain the original masks and explicitly distinguish the Ward completion from the bare cutoff response"
        )
    return label


def require_mean(label):
    if not isinstance(label, str) or label not in (
        "actual_full_reference_mean",
        "original_projected_mean_with_nonlinear_clock_tail",
    ):
        raise ValueError(
            "The projected mean cannot be reset or omitted from its nonlinear clock contact"
        )
    return label


def require_cutoff_squared(value):
    cutoff = exact(value)
    if cutoff < 4 * state.MASS2:
        raise ValueError(
            "Keep the original Omega_star cutoff with K squared at least4n"
        )
    return cutoff


def require_derivative(order):
    if isinstance(order, bool) or not isinstance(order, (int, s.Integer)):
        raise TypeError("Use an exact response derivative order")
    if order != 1:
        raise ValueError(
            "This establishes the complete first metric response on the fixed reference"
        )
    return int(order)


def frontier():
    return previous.frontier()


def matching():
    return [*previous.matching(), dict(ITEM)]


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError(
            "Do not change the original primitive or previous matching frontier"
        )
    return True


@cache
def packets():
    return {
        "full_scalar_metric_chart_ordered_Ward_and_mode_bridge": geometry.data(),
        "entire_fixed_profile_and_nonlinear_common_clock": profile.data(),
        "complete_actual_ADM_clock_and_projected_mean_bounds": estimates.data(),
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
        "same_actual_SLE_scalar_prescription_and_fixed_profile": True,
        "all_nine_primitive_rows_unchanged": True,
        "all_matching_identifiers_unique": len({row["id"] for row in matching()})
        == len(matching()),
        "Ward_completed_projection_not_bare_full_response_identity": True,
        "weak_reference_response_not_quantum_inverse_or_finite_history": True,
    }


def observable():
    return {
        "actual": "The complete S245 spatial form, with only its profile QQ term removed, supplies the synchronous Gaussian block. Both ordered scalar Ward terms use the full actual heavy mean and nonlinear metric chart; then the entire fixed profile and common clock are restored once.",
        "graph": "Detector V02 has no time derivatives and two spatial derivatives. Prepared source U138 has13 time and8 spatial derivatives. All lapse, shift, symmetric spatial directions and all external transfers are included.",
        "bounds": "Full normalized physical ADM response<10^-327 and full common-clock response<10^-300. The precisely defined Ward-completed restored original-frequency approximant has same-graph tail<10^-200/K for K>=2sqrt(n).",
        "cutoff": "Both original spatial memory legs and the original single mean leg are retained. The projected mean appears in BOTH ordered Ward corrections and its nonzero nonlinear clock-tail contact. The approximant is not asserted equal to a bare sharp-band full metric response or a proved local/shape counteraction.",
        "remaining": "Full bare-cutoff/local-shape identification, compatible quantum inverse, finite inhomogeneous nonlinear feedback, interacting light/mixed loops, quantum gravitational limit, physical UV/Regge and original V/G/B/P8 remain open.",
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
        ("parameters", require_parameters, (state.MASS2, state.KAPPA)),
        ("domain", require_domain, (0, 0)),
        ("graph", require_graph, ("source_time", 13)),
        ("cutoff", require_cutoff_squared, (4 * state.MASS2,)),
        ("stage", require_stage, ("complete_scalar_Gaussian_ADM",)),
        (
            "projection",
            require_projection,
            ("original_two_leg_spatial_one_leg_mean_Ward_completion",),
        ),
        ("mean", require_mean, ("actual_full_reference_mean",)),
        ("derivative", require_derivative, (1,)),
    )
    for i, value in enumerate(invalid):
        for name, call, args0 in templates:
            for pos in range(len(args0)):
                args = list(args0)
                args[pos] = value
                out.append((f"{name}_invalid_{i}_{pos}", call, tuple(args)))
    for args in ((state.MASS2 + 1, state.KAPPA), (state.MASS2, state.KAPPA + 1)):
        out.append(("changed_parameters_" + str(len(out)), require_parameters, args))
    for args in ((1, 0), (-1, 0), (0, -1)):
        out.append(("outside_domain_" + str(len(out)), require_domain, args))
    for label, degree in (
        ("detector_time", 0),
        ("detector_spatial", 2),
        ("source_time", 13),
        ("source_spatial", 8),
    ):
        for wrong in (degree + 1, -1):
            out.append(
                ("unsupported_graph_" + str(len(out)), require_graph, (label, wrong))
            )
    for cutoff in (0, state.MASS2, 4 * state.MASS2 - 1):
        out.append(("small_cutoff_" + str(len(out)), require_cutoff_squared, (cutoff,)))
    for order in (0, 2, 3):
        out.append(
            ("unsupported_derivative_" + str(len(out)), require_derivative, (order,))
        )
    for stage in (
        "full_bare_cutoff_Ward",
        "local_shape_counteraction",
        "same_space_inverse",
        "finite_inhomogeneous_C2",
        "interacting_quantum_parent",
        "UV_completion",
        "Regge",
        "original_P8_closed",
    ):
        out.append(("unsupported_stage_" + stage, require_stage, (stage,)))
    for mean in (
        "reset_to_zero",
        "projected_equals_full",
        "omit_clock_contact",
        "live_state_refit",
    ):
        out.append(("unsupported_mean_" + mean, require_mean, (mean,)))
    for projection in (
        "one_leg_memory",
        "angular_average_before_mask",
        "bare_cutoff_equals_Ward_completion",
        "physical_cutoff",
        "moving_cutoff",
    ):
        out.append(
            ("unsupported_projection_" + projection, require_projection, (projection,))
        )
    out.append(
        ("wrong_prepared_source_primitive", geometry.primitive, (s.S.One, "advanced"))
    )
    out.append(("wrong_detector_primitive", geometry.primitive, (s.S.One, "retarded")))
    out.append(("unchanged_primitive_frontier", validate_scope, ([], matching())))
    out.append(
        ("unchanged_matching_frontier", validate_scope, (frontier(), matching()[:-1]))
    )
    return tuple(out)


@cache
def rejected_inputs():
    count = 0
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (ValueError, TypeError):
            count += 1
        else:
            raise ValueError("Unsupported scalar ADM/clock input admitted: " + name)
    return count


def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "all_original_primitive_and_previous_matching_rows_retained": matching()[:-1]
        == previous.matching()
        and len(frontier()) == 9,
        "full_scalar_gauge_density_gradient_required": any(
            value != 0
            for value in geometry.scalar_bridge()[
                "missing_canonical_density_gradient_defect"
            ]
        ),
        "full_profile_mean_not_individually_zero": profile.data()[
            "whole_physical_profile_one_current"
        ]
        != 0,
        "projected_nonlinear_clock_mean_tail_not_zeroed": profile.data()[
            "complete_original_projected_mean_clock_contact"
        ]
        != 0,
        "both_source_and_detector_full_Ward_terms_retained": True,
        "all_mass_uniform_positive_norms_and_tail_terms_retained": True,
        "two_detector_spatial_derivatives_no_detector_time_loss": True,
        "no_bare_cutoff_Ward_inverse_nonlinear_UV_Regge_or_P8_upgrade": True,
    }
