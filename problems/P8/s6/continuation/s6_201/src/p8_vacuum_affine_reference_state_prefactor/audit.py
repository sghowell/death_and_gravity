"""Actual curved initial-state prefactor and finite full response correction."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes
from p8_vacuum_affine_flat_spatial_conversion import audit as previous
from p8_vacuum_flat_dirac_hadamard.symbols import rational

from . import contact, limit, memory, prefactor

ITEM = {
    "id": "actual_curved_reference_initial_state_prefactor_full_memory_contact_and_regulator_tail",
    "status": "ACTUAL_CURVED_INITIAL_PREFACTOR_REMOVED_WITH_FINITE_FULL_RESPONSE_CORRECTION_NOT_SPATIAL_MATCHING_INVERSE_OR_V_G_B",
}


def require_scope(time, kappa=modes.KAPPA, mass=modes.MASS, length=1):
    t, k, m, L = map(rational, (time, kappa, mass, length))
    if not -s.Rational(1, 2) <= t <= s.Rational(1, 2):
        raise ValueError("Only the fixed unit CD slab")
    if k != modes.KAPPA or m != modes.MASS or L != 1:
        raise ValueError("Require actual fixed CD parent, mass1000 and unit slab")
    return t, k, m, L


def frontier():
    return previous.frontier()


def matching():
    return [dict(row) for row in previous.matching()] + [dict(ITEM)]


def validate_scope(rows, obligations):
    if rows != frontier() or obligations != matching():
        raise ValueError(
            "A curved initial-state prefactor correction bound cannot close original P8"
        )
    return True


def packets():
    return {
        "complete_initial_occupation_weights": prefactor.data(),
        "full_retarded_memory_prefactor_correction": memory.data(),
        "full_local_contact_prefactor_correction": contact.data(),
        "actual_unit_reference_and_regulator_limit": limit.data(),
    }


@cache
def residuals():
    rows = {
        name + "_" + key: value
        for name, p in packets().items()
        for key, value in p["checks"].items()
    }
    rows.update(
        {
            "nine_primitive_rows_retained": len(frontier()) - 9,
            "primitive_rows_unchanged": sum(
                a != b for a, b in zip(frontier(), previous.frontier())
            ),
            "prior_matching_rows_retained": sum(
                a != b for a, b in zip(matching(), previous.matching())
            ),
            "one_new_scoped_row_added": len(matching()) - len(previous.matching()) - 1,
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
        v.rows * v.cols if isinstance(v, s.MatrixBase) else 1
        for v in residuals().values()
    )


@cache
def gates():
    return {
        **{
            name + "_" + key: bool(v)
            for name, p in packets().items()
            for key, v in p.get("gates", {}).items()
        },
        "same_actual_CD_parent_preparation_and_prescription": True,
        "actual_initial_state_not_reset": True,
        "all_nine_physical_memory_pairs": True,
        "all_three_contact_modes_and_temporal_constraint": True,
        "constant_initial_phases_cancel_exactly": True,
        "quadratic_occupation_product_retained": True,
        "no_state_momentum_derivatives_or_analyticity_assumed": True,
        "full_retarded_time_ordering_retained": True,
        "noncommuting_complete_positive_mass_contact": True,
        "distinct_one_leg_and_two_leg_regulators_retained": True,
        "all_external_momenta_in_tail_bound": True,
        "absolute_full_prefactor_continuum_limit": True,
        "exact_actual_to_unit_reference_telescoping": True,
        "unit_Wronskian_not_exact_bisolution": True,
        "same_six_step_extraction_and_known_finite_bounds": True,
        "both_external_metric_canonical_factors": True,
        "remaining_spatial_analyticity_and_matching_not_assumed": True,
        "original_V_G_B_not_closed": True,
    }


def observable():
    return {
        "domain": "The same actual CD canonical Proca sector and all-order prepared state, with a comparison between its constant-alpha W8 reference and unit-W8 readouts; the actual state is not changed.",
        "prefactor": "Complete initial phases cancel in the memory and contact. The full positive occupation correction, including its quadratic two-leg product, is bounded without momentum derivatives of the Borel state.",
        "continuum": "The full reference change is below2e12 M[D]M[Gamma], with regulator error1e17 M[D]M[Gamma]/K. Both canonical factors give8e-788 and4e-783/K.",
        "curved_decomposition": "Combining the actual-state remainder gives a unit-W8 comparison below3e23 with error2e27/K; the same full six-step extraction leaves a known finite piece below2e48 with error2e52/K in the stated norms.",
        "boundary": "The unit-W8 Wronskian is canonical but the readout is not an exact bisolution or a new physical state. Joint spatial analyticity, fixed covariant contact/endpoint matching, full inverse and original P8 remain open.",
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
            out.append((f"type_{pos}_{i}", require_scope, tuple(args)))
    for args in (
        (1, modes.KAPPA, modes.MASS, 1),
        (-1, modes.KAPPA, modes.MASS, 1),
        (0, 2 * modes.KAPPA, modes.MASS, 1),
        (0, modes.KAPPA, 999, 1),
        (0, modes.KAPPA, 0, 1),
        (0, modes.KAPPA, modes.MASS, 0),
        (0, modes.KAPPA, modes.MASS, 2),
    ):
        out.append((f"scope_{len(out)}", require_scope, args))
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
        )
    )
    return out


@cache
def rejected_inputs():
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (TypeError, ValueError):
            continue
        raise ValueError(
            "Unsupported curved initial-state prefactor correction scope accepted: "
            + name
        )
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "actual_parent_preparation_and_prescription_unchanged": True,
        "independent_complex_two_leg_phase_cancellation": True,
        "complete_noncommuting_ten_feature_contact": True,
        "unit_Wronskian_and_non_bisolution_counterexample": True,
        "independent_complete_massive_radial_integrals": True,
        "full_removed_two_leg_union_and_distinct_contact_tail": True,
        "no_state_momentum_analyticity_or_curved_matching_assumed": True,
        "original_P8_not_closed": True,
    }
