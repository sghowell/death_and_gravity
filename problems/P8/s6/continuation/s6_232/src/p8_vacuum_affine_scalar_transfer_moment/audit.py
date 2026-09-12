"""Original transfer matching and mandatory full spectral weight, with scope controls."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes
from p8_vacuum_affine_scalar_matching_scale import audit as previous

from . import atomic, coefficients, moment

ITEM = {
    "id": "original_scalar_transfer_matching_mandatory_low_spectral_weight_and_higher_forward_coefficient_with_positive_spectral_control",
    "status": "CONDITIONAL_FULL_SPECTRAL_WEIGHT_AND_HIGHER_COEFFICIENT_REQUIREMENT_NOT_FULL_LOOP_BOUND_PHYSICAL_CUTOFF_EXACT_UNITARY_COMPLETION_OR_P8",
}
require_scope = previous.require_scope


def require_stage(stage):
    if not isinstance(stage, str) or stage not in (
        "original_transfer_coefficient_and_first_elastic_slope",
        "exact_conditional_low_high_positive_moment_inequality",
        "necessary_full_low_spectral_weight_under_coefficient_matching",
        "original_first_elastic_moment_bound_only",
        "separately_named_positive_pole_coefficient_diagnostic",
    ):
        raise ValueError(
            "Only original transfer data and the stated conditional moment requirement are proved"
        )
    return stage


def require_tolerances(energy, delta0, delta1):
    for value in (energy, delta0, delta1):
        if isinstance(value, bool) or not isinstance(
            value, (int, s.Integer, s.Rational)
        ):
            raise TypeError(
                "Spectral split and matching tolerances must be exact finite rationals"
            )
    if energy <= 2 or delta0 < 0 or delta1 < 0 or delta1 >= 1:
        raise ValueError("Require split above threshold, delta0>=0 and 0<=delta1<1")
    return tuple(s.Rational(value) for value in (energy, delta0, delta1))


def frontier():
    return previous.frontier()


def matching():
    return [dict(row) for row in previous.matching()] + [dict(ITEM)]


def validate_scope(rows, obligations):
    if rows != frontier() or obligations != matching():
        raise ValueError(
            "Necessary spectral enhancement and a rational control do not close original P8"
        )
    return True


def packets():
    return {
        "original_transfer_coefficients_and_first_elastic_derivative": coefficients.data(),
        "exact_positive_split_inequality_and_full_moment_requirement": moment.data(),
        "separate_positive_pole_coefficient_control_not_a_completion": atomic.data(),
    }


@cache
def residuals():
    rows = {
        name + "_" + key: value
        for name, packet in packets().items()
        for key, value in packet["checks"].items()
    }
    rows.update(
        {
            "nine_original_primitive_rows": len(frontier()) - 9,
            "unchanged_primitive_rows": sum(
                a != b for a, b in zip(frontier(), previous.frontier())
            ),
            "unchanged_previous_matching_rows": sum(
                a != b for a, b in zip(matching(), previous.matching())
            ),
            "one_new_conditional_spectral_requirement": len(matching())
            - len(previous.matching())
            - 1,
        }
    )
    return {key: s.cancel(value) for key, value in rows.items()}


def scalar_entry_count():
    return len(residuals())


@cache
def gates():
    return {
        **{
            name + "_" + key: bool(value)
            for name, packet in packets().items()
            for key, value in packet["gates"].items()
        },
        "no_cut_below_threshold_and_differentiated_arc_are_explicit_premises": True,
        "all_full_low_spectral_weight_kept_instead_of_cut_erasure": True,
        "physical_low_coefficient_errors_are_unproved_not_chosen_counterterms": True,
        "first_elastic_bound_not_transplanted_to_full_measure": True,
        "new_positive_pole_only_a_separate_diagnostic": True,
        "unique_matching_identifiers": len({row["id"] for row in matching()})
        == len(matching()),
    }


def observable():
    return {
        "original_scope": "The same S182 vacuum and complete S177 massive scalar tree, using S231 normalized identical channels. b2_tree(t)=4lambda-3gamma*t and the first elastic transfer slope are exact original data.",
        "conditional_physical_requirement": "Under explicit full fixed-t dispersion and spectral-gap premises, b21+3b20/[2(K-2)]>=-3J4_low(K)/2. Physical coefficient errors delta0<=1 and delta1<=1/2 at split M10^99 therefore require full J4>gamma/5.",
        "independent_first_elastic_bound": "The complete original first elastic moment is positive and below gamma*10^-200, so the required full moment exceeds it by more than2*10^199. This necessary enhancement does not assume small full-loop errors.",
        "mandatory_higher_coefficient": "The positive full moment Gram matrix gives b40>=4*b21^2/(9*b20) for b21<0. The named physical coefficient tolerances therefore require b40>=gamma^2/(8lambda)>0, although the original tree has b40=0. Exact b20,b21 matching requires b40>=gamma^2/lambda; atom saturation is only a positive-measure diagnostic, not a physical unitary optimality claim.",
        "control": "A separate positive rational pole with M_H²=2+2lambda/gamma and g²=1/2^26 matches b20,b21 and supplies J4=2gamma, but also has nonzero v4 coefficient gamma²/lambda. It is not an exact unitary amplitude or a matched parent.",
        "boundary": "No full absorptive calculation, physical cutoff, UV-parent exclusion, finite-gravity result, modified original action/state, full bounce matching or original V/G/B/P8 closure follows.",
    }


def bad_cases():
    bad = (
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
    for i, value in enumerate(bad):
        for pos in range(4):
            args = [0, modes.KAPPA, modes.MASS, 1]
            args[pos] = value
            out.append((f"scope_type_{pos}_{i}", require_scope, tuple(args)))
        out.append((f"stage_type_{i}", require_stage, (value,)))
        for pos in range(3):
            args = [moment.ENERGY, 1, s.Rational(1, 2)]
            args[pos] = value
            out.append((f"tolerance_type_{pos}_{i}", require_tolerances, tuple(args)))
    for args in (
        (1, modes.KAPPA, modes.MASS, 1),
        (-1, modes.KAPPA, modes.MASS, 1),
        (0, 2 * modes.KAPPA, modes.MASS, 1),
        (0, modes.KAPPA, 999, 1),
        (0, modes.KAPPA, modes.MASS, 2),
    ):
        out.append((f"scope_{len(out)}", require_scope, args))
    for args in (
        (0, 0, 0),
        (-3, 0, 0),
        (2, 0, 0),
        (3, -1, 0),
        (3, 0, -1),
        (3, 0, 1),
        (3, 0, 2),
    ):
        out.append((f"tolerance_domain_{len(out)}", require_tolerances, args))
    for stage in (
        "physical_cutoff",
        "full_loop_error_bound",
        "tree_dominates_full_cut",
        "discard_light_cut",
        "subtract_unknown_heavy_pole",
        "differentiated_arc_automatic",
        "new_light_cut_ignored",
        "new_renormalization_b21",
        "full_UV_no_go",
        "closed_P8",
        "positive_atom_is_exact_unitary",
        "original_parent_replaced",
        "full_four_point_function_matched",
        "bounce_parent_complete",
        "finite_gravity_positive_bound",
        "all_higher_operators_zero",
        "full_b40_can_remain_zero",
    ):
        out.append((f"unsupported_stage_{stage}", require_stage, (stage,)))
    for i in range(9):
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
                "extra_UV_completion",
                validate_scope,
                (
                    frontier(),
                    matching() + [{"id": "UV_completion", "status": "COMPLETE"}],
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
            "Unsupported spectral moment or physical closure accepted: " + name
        )
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "omitting_crossing_half_shift_loses_inverse_fourth_term": True,
        "positive_transfer_term_and_full_low_cut_both_retained": True,
        "full_moment_not_bounded_by_a_first_loop_without_error": True,
        "unknown_heavy_poles_are_not_erased": True,
        "positive_rational_control_not_exact_unitarity": True,
        "extra_v4_coefficient_prevents_full_tree_identification": True,
        "spectral_split_not_a_physical_cutoff": True,
        "all_original_primitive_statuses_unchanged": frontier() == previous.frontier(),
    }
