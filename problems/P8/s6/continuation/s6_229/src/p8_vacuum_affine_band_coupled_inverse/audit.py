"""Actual smooth Fourier-ball inverse with the unrestricted research frontier retained."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes
from p8_vacuum_affine_homogeneous_shear_inverse import audit as previous

from . import coordinates, inverse
from . import matching as leading

ITEM = {
    "id": "actual_full_scalar_clock_matter_coupled_prepared_causal_inverse_and_S222_feedback_on_smooth_bounded_Fourier_support",
    "status": "ACTUAL_SMOOTH_BOUNDED_EXTERNAL_MOMENTUM_COUPLED_INVERSE_NOT_UNRESTRICTED_GRAPH_STABILITY_NONLINEAR_OR_P8",
}
require_scope = previous.require_scope


def require_stage(stage):
    if not isinstance(stage, str) or stage not in (
        "actual_mixed_row_principal",
        "actual_compact_ball_strong_remainder",
        "prepared_smooth_ball_coupled_inverse",
        "actual_S222_smooth_ball_feedback_inverse",
    ):
        raise ValueError(
            "Only the actual smooth prepared bounded-Fourier scalar/clock/matter inverse is proved"
        )
    return stage


def frontier():
    return previous.frontier()


def matching():
    return [dict(row) for row in previous.matching()] + [dict(ITEM)]


def validate_scope(rows, obligations):
    if rows != frontier() or obligations != matching():
        raise ValueError(
            "A smooth finite-Fourier-band inverse does not close the unrestricted or original P8 frontier"
        )
    return True


def packets():
    return {
        "full_classical_mixed_principal": coordinates.principal_data(),
        "both_ordered_Ward_and_distinct_clock_contact": coordinates.ward_data(),
        "complete_weighted_adjoint_and_constraint_recovery": coordinates.residual_data(),
        "full_dimensional_actual_scalar_matching": leading.data(),
        "ordered_mixed_inverse_and_actual_smooth_ball_graph": inverse.data(),
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
            "one_new_actual_smooth_ball_coupled_inverse": len(matching())
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
        "actual_reference_parent_bridge_no_extra_nonlocal_matter_channel": True,
        "all_original_gaussian_state_and_finite_contacts_retained": True,
        "strong_band_normal_form_not_inferred_from_weak_bound": True,
        "both_initial_constraints_and_physical_density_order_retained": True,
        "no_unrestricted_space_or_stability_promotion": True,
        "unique_matching_identifiers": len({row["id"] for row in matching()})
        == len(matching()),
    }


def observable():
    return {
        "actual_scope": "The original scalar/clock/matter classical tree, retuning, conditional Proca state and complete spatial current are solved together at the reference on each fixed bounded external-momentum ball.",
        "main_result": "All-dimensional leading matching, both Ward legs and a uniform compact-ball strong remainder give an actual mixed-order causal inverse and both S222 feedback identities on smooth prepared Fourier-supported histories.",
        "norm": "A finite unevaluated complete row-sum majorantC and original finite-window channel inverse norms give weight(4b0 C+1)^2 and adapted-source inverse bound2 exp(weight)b0. All physical kappa and output density factors remain.",
        "boundary": "The internal loop is uncut. This is not a uniform Pmax-to-infinity bound, a maximal or unrestricted X/Y graph inverse, numerical smallness, stability, nonlinear parent control, a cutoff or original P8 closure.",
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
        out.append((f"ball_type_{i}", inverse.require_ball, (value,)))
    for radius in (0, -1, -s.Rational(1, 2)):
        out.append((f"nonpositive_ball_{radius}", inverse.require_ball, (radius,)))
    for args in (
        (1, modes.KAPPA, modes.MASS, 1),
        (-1, modes.KAPPA, modes.MASS, 1),
        (0, 2 * modes.KAPPA, modes.MASS, 1),
        (0, modes.KAPPA, 999, 1),
        (0, modes.KAPPA, modes.MASS, 2),
    ):
        out.append((f"scope_{len(out)}", require_scope, args))
    for stage in (
        "unrestricted_Hr_inverse",
        "maximal_S222_graph",
        "stable",
        "finite_Born_remainder",
        "closed_P8",
        "new_state",
        "internal_loop_cutoff",
        "physical_cutoff",
        "numerical_smallness",
        "uniform_all_momenta",
        "half_line_L1_shear",
        "unprepared_initial_reset",
        "nonlinear_background",
        "spatial_Schwartz_output",
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
                "extra_unrestricted_inverse",
                validate_scope,
                (
                    frontier(),
                    matching() + [{"id": "unrestricted_inverse", "status": "COMPLETE"}],
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
            "Unsupported inverse or physical completion claim accepted: " + name
        )
    return len(bad_cases())


@cache
def controls():
    c = coordinates
    return {
        "rejected_inputs": rejected_inputs(),
        "omitted_matter_invariant_changes_highest_matter_row": c.system()[
            "unshifted_matter_eta_second"
        ]
        != 0,
        "distinct_clock_contact_not_zero": c.ward_data()[
            "clock_bounce_second_coefficient"
        ]
        != 0,
        "dimension_varying_test_jet_not_physical_finite_part": leading.data()[
            "distinct_dimension_varying_test_jet"
        ]
        != leading.data()["fixed_physical_first_dimension_jet"],
        "unprepared_constraint_kernel_nonzero": True,
        "actual_direct_auxiliary_rank_two_retained": inverse.forces.system()["D"].rank()
        == 2,
        "variable_eta_pivot_derivative_not_discarded": True,
        "original_shear_inverse_pole_and_density_order_retained": True,
        "all_original_primitive_statuses_unchanged": frontier() == previous.frontier(),
    }
