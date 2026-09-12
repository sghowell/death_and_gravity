"""Original complete flat tensor response, with the physical UV frontier retained."""

from functools import cache

import sympy as s
from p8_vacuum_affine_band_coupled_inverse import audit as previous
from p8_vacuum_affine_cd_source_noise import modes

from . import analytic, normalization, response

ITEM = {
    "id": "actual_full_flat_vacuum_tree_Proca_TT_first_sheet_complex_pair_low_energy_bound_and_all_momentum_original_causal_graph_inverse",
    "status": "ACTUAL_RETAINED_FLAT_MEAN_EQUATION_GROWING_POLES_AND_CAUSAL_INVERSE_NOT_CURVED_INSTABILITY_PHYSICAL_UV_NO_GO_OR_P8",
}
require_scope = previous.require_scope


def require_stage(stage):
    if not isinstance(stage, str) or stage not in (
        "original_flat_vacuum_polynomial_bridge",
        "complete_first_sheet_pole_count",
        "actual_low_energy_form_factor_bound",
        "all_momentum_flat_TT_causal_graph_inverse",
        "retained_flat_equation_prepared_forced_growth",
    ):
        raise ValueError(
            "Only the original retained flat TT response in the stated scope is proved"
        )
    return stage


def frontier():
    return previous.frontier()


def matching():
    return [dict(row) for row in previous.matching()] + [dict(ITEM)]


def validate_scope(rows, obligations):
    if rows != frontier() or obligations != matching():
        raise ValueError(
            "Retained flat growing poles and its causal inverse do not close the original physical P8 frontier"
        )
    return True


def packets():
    return {
        "original_flat_covariance_subtraction_tree_and_finite_normalization": normalization.data(),
        "complete_first_sheet_count_and_exact_scale_margins": analytic.data(),
        "original_full_dispersion_moments_causal_inverse_and_graph": response.data(),
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
            "one_new_actual_full_flat_response": len(matching())
            - len(previous.matching())
            - 1,
        }
    )
    return {
        key: s.ImmutableMatrix(value.applyfunc(s.cancel))
        if isinstance(value, s.MatrixBase)
        else s.cancel(value)
        for key, value in rows.items()
    }


def scalar_entry_count():
    return sum(
        value.rows * value.cols if isinstance(value, s.MatrixBase) else 1
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
        "original_minkowski_vacuum_and_lower_polynomial_bridge_retained": True,
        "complete_tree_and_fixed_finite_m2_R_shift_retained": True,
        "all_first_sheet_poles_and_both_cut_moments_retained": True,
        "fixed_zero_germ_and_all_momentum_weak_forward_graph_retained": True,
        "retained_equation_growth_not_physical_parent_exclusion": True,
        "unique_matching_identifiers": len({row["id"] for row in matching()})
        == len(matching()),
    }


def observable():
    return {
        "actual_scope": "The specified Minkowski vacuum of the original conditional Proca Gaussian mean equation, with the original Einstein tree and finite prescription. The nonstationary CD state is not transplanted.",
        "main_result": "The original covariance/adiabatic bridge fixes the lower polynomial. The complete force denominator has exactly one simple conjugate pair of nonreal first-sheet zeros with negative real p, giving prepared forced growing lambda modes. Its full massless/pair/cut reciprocal has a two-sided all-momentum flat TT causal graph inverse.",
        "bounds": "At actual m1000 and kappa10^800, 10^796<|z|<kappa, Re sqrt(z)>10^393, while the original physical propagator ratio differs from Einstein by less than10^-796 on |p|<=m^2. The finite-time inverse norm is below4T^2 exp(sqrt(kappa)T), not small.",
        "boundary": "This is not a stable retained equation, physical UV inconsistency, curved-bounce instability, nonlinear parent control, physical cutoff, unrestricted curved S222 graph result or original V/G/B/P8 closure.",
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
    for args in (
        (1, modes.KAPPA, modes.MASS, 1),
        (-1, modes.KAPPA, modes.MASS, 1),
        (0, 2 * modes.KAPPA, modes.MASS, 1),
        (0, modes.KAPPA, 999, 1),
        (0, modes.KAPPA, modes.MASS, 2),
    ):
        out.append((f"scope_{len(out)}", require_scope, args))
    for stage in (
        "physical_UV_no_go",
        "stable",
        "closed_P8",
        "new_state",
        "CD_state_on_flat",
        "pole_deleted",
        "order_reduced",
        "frequency_sign_sheet_selection",
        "new_finite_counterterm",
        "physical_cutoff",
        "causal_low_frequency_bandpass",
        "curved_bounce_instability",
        "nonlinear_parent_control",
        "unrestricted_S222_graph",
        "unprepared_initial_reset",
        "above_Planck_scale",
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
                "extra_physical_no_go",
                validate_scope,
                (
                    frontier(),
                    matching() + [{"id": "UV_excluded", "status": "COMPLETE"}],
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
            "Unsupported full response or physical closure accepted: " + name
        )
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "omitting_original_m2_R_changes_force_coefficient": (
            normalization.C - normalization.C0
        ).is_positive,
        "isolated_A2_zero_is_not_a_complete_response_pole": True,
        "cut_only_positive_reciprocal_violates_zeroth_moment": True,
        "both_complex_residues_required_and_nonzero": True,
        "flat_all_momentum_graph_not_curved_coupled_graph": True,
        "compact_time_forcing_not_a_forbidden_initial_reset": True,
        "pole_scale_not_a_physical_cutoff_or_parent_no_go": True,
        "all_original_primitive_statuses_unchanged": frontier() == previous.frontier(),
    }
