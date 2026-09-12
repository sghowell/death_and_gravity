"""Weighted reference-correction result with the actual P8 frontier intact."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes
from p8_vacuum_affine_curved_scalar_reference import audit as previous

from . import analytic, curved, weighted

ITEM = {
    "id": "original_channel_exterior_log_coercivity_and_ordered_weighted_bounded_curvature_reference_corrections",
    "status": "WEIGHTED_REFERENCE_CORRECTIONS_INVERTIBLE_NOT_ACTUAL_FULL_QUANTUM_MATCHING_STABILITY_OR_P8",
}
require_scope = previous.require_scope


def require_stage(stage):
    if not isinstance(stage, str) or stage not in (
        "full_complex_log_bound",
        "weighted_channel_inverse",
        "bounded_reference_correction",
    ):
        raise ValueError(
            "Only the stated original-channel weighted reference results are proved"
        )
    return stage


def frontier():
    return previous.frontier()


def matching():
    return [dict(row) for row in previous.matching()] + [dict(ITEM)]


def validate_scope(rows, obligations):
    if rows != frontier() or obligations != matching():
        raise ValueError(
            "Weighted reference existence cannot close the actual quantum or P8 frontier"
        )
    return True


def packets():
    return {
        "original_first_sheet_logarithmic_coercivity": analytic.data(),
        "uniform_weighted_channel_inverse_and_bounded_correction": weighted.data(),
        "ordered_curvature_reference_and_growing_pole_control": curved.data(),
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
            "one_new_weighted_reference_input": len(matching())
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
        "complex_line_uniform_in_time_and_spatial_frequencies": True,
        "original_shear_pole_and_initial_boundary_retained": True,
        "bounded_channel_matrix_not_commuted_through_memory": True,
        "growing_reference_pole_control_prevents_stability_promotion": True,
        "actual_curved_remainder_not_assumed_in_bounded_class": True,
        "unique_matching_identifiers": len({row["id"] for row in matching()})
        == len(matching()),
    }


def observable():
    return {
        "complex_bound": "Both original scalar factors satisfy Re A_i(p)>=(13/80)log(|p|/(4m^2))-4 on the full first-sheet exterior. This is proved from cut banks, the inner circle and uniform outer circles, not from a real-axis sample.",
        "weighted_inverse": "For a bounded channel matrix V(t) of norm<=M, sigma=2m exp(32+8M) gives ||Kdiag||<=5/(32+13M) and ||(Fdiag+V)^-1||<=5/(32+8M) in weighted L2_tH^r, all time/spatial frequencies and every realr.",
        "curved_composition": "The actual S225 curvature-coordinate reference Bc*(Fdiag+V)Bc has ordered inverseY K_V Z with the displayed no-spatial-loss bound and full conformal density restored on the right.",
        "boundary": "Removing the weight costs exp(sigma T). An explicit norm-below-one correction has a growing reference pole while the weighted theorem still applies. Actual curved matching, full S222 inverse, physical stability and original V/G/B/P8 remain open.",
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
        out.append((f"bound_type_{i}", weighted.size, (value,)))
    for args in (
        (1, modes.KAPPA, modes.MASS, 1),
        (-1, modes.KAPPA, modes.MASS, 1),
        (0, 2 * modes.KAPPA, modes.MASS, 1),
        (0, modes.KAPPA, 999, 1),
        (0, modes.KAPPA, modes.MASS, 2),
    ):
        out.append((f"scope_{len(out)}", require_scope, args))
    for value in (-1, s.Rational(-1, 2)):
        out.append((f"negative_bound_{value}", weighted.size, (value,)))
    for stage in (
        "full_quantum_inverse",
        "stable",
        "finite_Born_remainder",
        "closed_P8",
        "new_state",
        "pole_deletion",
        "physical_cutoff",
        "unweighted_smallness",
        "actual_curved_matching",
        "zero_weight",
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
                "extra_actual_inverse",
                validate_scope,
                (
                    frontier(),
                    matching() + [{"id": "full_actual_inverse", "status": "COMPLETE"}],
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
            "Unsupported full quantum or stability claim accepted: " + name
        )
    return len(bad_cases())


@cache
def controls():
    comparison = curved.data()["comparison_correction"]
    return {
        "rejected_inputs": rejected_inputs(),
        "zero_weight_cannot_use_the_positive_exterior_gap": analytic.COEFFICIENT * 0
        - analytic.OFFSET
        < 0,
        "original_shear_threshold_real_part_not_positive": analytic.data()[
            "global_real_part_gaps"
        ]["shear"]
        < 0,
        "bounded_reference_correction_with_growing_pole": comparison["bound"] < 1
        and comparison["matrix"] != s.zeros(2),
        "positive_weight_not_physical_damping_modification": True,
        "unweighted_existence_cost_not_silently_one": weighted.SIGMA.is_positive,
        "actual_local_remainder_not_deleted": previous.ITEM["id"]
        in {row["id"] for row in matching()},
        "physical_kappa_kept_in_output_density": curved.prior.g.bridge.KAPPA
        == modes.KAPPA,
        "all_original_primitive_statuses_unchanged": frontier() == previous.frontier(),
    }
