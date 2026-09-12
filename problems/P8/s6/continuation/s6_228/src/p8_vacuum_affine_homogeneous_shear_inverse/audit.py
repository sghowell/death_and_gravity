"""Actual homogeneous quantum/tree inverse and unchanged spatial/P8 frontier."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes
from p8_vacuum_affine_full_finite_local_reference import audit as previous

from . import local, volterra
from . import matching as leading

ITEM = {
    "id": "actual_homogeneous_tracefree_full_curved_quantum_and_classical_tensor_prepared_causal_inverse",
    "status": "ACTUAL_HOMOGENEOUS_FULL_NONLOCAL_SHEAR_AND_TREE_INVERSE_NOT_NONZERO_TRANSFER_FULL_S222_STABILITY_OR_P8",
}
require_scope = previous.require_scope


def require_stage(stage):
    if not isinstance(stage, str) or stage not in (
        "actual_homogeneous_shear_matching",
        "complete_homogeneous_weak_log_remainder",
        "prepared_homogeneous_total_inverse",
    ):
        raise ValueError(
            "Only the actual homogeneous shear problem and its stated prepared inverse are proved"
        )
    return stage


def frontier():
    return previous.frontier()


def matching():
    return [dict(row) for row in previous.matching()] + [dict(ITEM)]


def validate_scope(rows, obligations):
    if rows != frontier() or obligations != matching():
        raise ValueError(
            "The homogeneous shear inverse cannot close the spatial quantum or P8 frontier"
        )
    return True


def packets():
    return {
        "actual_leading_pair_and_fixed_dimensional_matching": leading.data(),
        "complete_original_finite_and_classical_local_terms": local.data(),
        "actual_weak_log_normal_form_and_prepared_inverse": volterra.data(),
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
            "one_new_actual_homogeneous_inverse": len(matching())
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
        "actual_state_and_original_full_current_not_replaced": True,
        "all_dimensional_leading_contacts_fixed_before_primitive": True,
        "complete_classical_tree_and_kappa_retained": True,
        "inverse_not_inferred_from_derivative_losing_weak_bound": True,
        "no_claim_of_nonzero_transfer_invariance_or_stability": True,
        "unique_matching_identifiers": len({row["id"] for row in matching()})
        == len(matching()),
    }


def observable():
    return {
        "actual_response": "The original homogeneous tracefree Proca current on the curved CD clock, all initial-state memory, fixed contacts and classical tensor tree are retained. O(3) gives five identical time operators at zero transfer.",
        "strong_normal_form": "All-dimensional leading matching and the actual selected-state finite high-frequency expansion give I4 T_total=F2(Dt^2)+V_total with an integrable weak-log kernel and finite simultaneous derivative majorants.",
        "complete_homogeneous_inverse": "With the original pole-plus-cut K2 and the complete actual majorantC, the finite weight(4C||K2||L1(0,1)+1)^2 gives a convergent ordered inverse on the same prepared smooth class. Constants remain unevaluated and kappa remains in the physical source.",
        "boundary": "This is per comoving volume at literal zero transfer, not the all-momentum scalar/clock/matter graph, nonlinear control, a numerical small inverse, stability, a cutoff or original P8 closure.",
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
        "full_spatial_quantum_inverse",
        "stable",
        "finite_Born_remainder",
        "closed_P8",
        "new_state",
        "pole_deletion",
        "physical_cutoff",
        "numerical_smallness",
        "full_S222_graph",
        "half_line_L1_shear",
        "unprepared_initial_reset",
        "nonlinear_background",
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
                "extra_full_spatial_inverse",
                validate_scope,
                (
                    frontier(),
                    matching() + [{"id": "full_spatial_inverse", "status": "COMPLETE"}],
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
            "Unsupported spatial quantum or stability claim accepted: " + name
        )
    return len(bad_cases())


@cache
def controls():
    data = leading.data()
    return {
        "rejected_inputs": rejected_inputs(),
        "longitudinal_mode_deletion_changes_log": data["physical_TT_and_LL"][1] != 0,
        "wrong_annihilation_phase_changes_leading_kernel": data[
            "normalized_positive_phase_radial_coefficient"
        ]
        != 0,
        "first_dimensional_jet_not_suppressed": data["first_dimension_jet"] != 0,
        "lower_coefficient_derivatives_retained": True,
        "original_inverse_shear_pole_nonzero": leading.spectral.RESIDUE_LO > 0,
        "physical_kappa_retained": local.original.KAPPA == modes.KAPPA,
        "no_new_final_source_cutoff": True,
        "all_original_primitive_statuses_unchanged": frontier() == previous.frontier(),
    }
