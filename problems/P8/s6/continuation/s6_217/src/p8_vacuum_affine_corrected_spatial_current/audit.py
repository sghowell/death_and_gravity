"""Corrected full tracefree current and restored conditional Ward input."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes
from p8_vacuum_affine_ordered_scalar_symbol import audit as previous
from p8_vacuum_flat_dirac_hadamard.symbols import rational

from . import assembly, response, subtraction

CORRECTIONS = {
    "complete_matched_tracefree_spatial_Gaussian_current_with_anchored_original_regulator": "HISTORICAL_PHASE_ERROR_SUPERSEDED_BY_CORRECTED_S217_FULL_TRACEFREE_CURRENT_AND_ANCHORED_REGULATOR",
    "prepared_ordered_Ward_reconstruction_full_chart_correction_and_three_remaining_scalar_kernels": "PREPARED_WARD_RECONSTRUCTION_WITH_CORRECTED_S217_KNOWN_TRACEFREE_INPUT_THREE_SCALAR_KERNELS_REMAIN_OPEN",
}
ITEM = {
    "id": "corrected_full_retarded_tracefree_spatial_current_original_regulator_and_prepared_Ward_input",
    "status": "CORRECTED_COMPLETE_TRACEFREE_SPATIAL_WEAK_CURRENT_AND_KNOWN_WARD_INPUT_NOT_FULL_SCALAR_REDUCED_INVERSE_BACKGROUND_OR_V_G_B",
}


def require_scope(time, kappa=modes.KAPPA, mass=modes.MASS, length=1):
    t, k, m, L = map(rational, (time, kappa, mass, length))
    if not -s.Rational(1, 2) <= t <= s.Rational(1, 2):
        raise ValueError("Only the original fixed unit CD slab")
    if k != modes.KAPPA or m != modes.MASS or L != 1:
        raise ValueError("Require actual kappa,mass1000 and unit slab")
    return t, k, m, L


def require_cutoff(cutoff, mass=modes.MASS):
    K, m = map(rational, (cutoff, mass))
    if m != modes.MASS or K < 2 * m:
        raise ValueError("Require the original auxiliary band K>=2m, mass1000")
    return K, m


def frontier():
    return previous.frontier()


def matching():
    rows = [dict(row) for row in previous.matching()]
    if set(CORRECTIONS) - {row["id"] for row in rows}:
        raise ValueError("A named historical restoration target is missing")
    for row in rows:
        if row["id"] in CORRECTIONS:
            row["status"] = CORRECTIONS[row["id"]]
    return rows + [dict(ITEM)]


def validate_scope(rows, obligations):
    if rows != frontier() or obligations != matching():
        raise ValueError("A corrected tracefree current cannot close original P8")
    return True


def packets():
    return {
        "independent_unequal_mode_covariance_and_CCR_phase_bridge": response.canonical_data(),
        "complete_corrected_endpoint_bulk_and_known_decomposition": response.decomposition_data(),
        "physical_UV_and_original_nondecaying_subtraction": subtraction.physical_data(),
        "both_analytic_phase_branches_and_dimension_limit": subtraction.dimension_data(),
        "corrected_full_current_homogeneous_anchor_and_regulator_error": assembly.data(),
        "corrected_prepared_Ward_known_input_and_initial_only_domain": assembly.ward_data(),
    }


@cache
def residuals():
    rows = {
        name + "_" + key: value
        for name, p in packets().items()
        for key, value in p["checks"].items()
    }
    before, after = previous.matching(), matching()
    rows.update(
        {
            "nine_original_primitive_rows": len(frontier()) - 9,
            "primitive_rows_unchanged": sum(
                a != b for a, b in zip(frontier(), previous.frontier())
            ),
            "exact_two_historical_recovery_status_changes": sum(
                a != b for a, b in zip(before, after)
            )
            - 2,
            "every_other_historical_row_unchanged": sum(
                a != b for a, b in zip(before, after) if a["id"] not in CORRECTIONS
            ),
            "one_new_corrected_current_checkpoint": len(after) - len(before) - 1,
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
            name + "_" + key: bool(value)
            for name, p in packets().items()
            for key, value in p["gates"].items()
        },
        "full_physical_current_not_only_finite_UV_adjustment": True,
        "independent_entire_retarded_integral_and_sixth_bulk_checks": True,
        "actual_unaveraged_odd_log_products_and_source_jets_retained": True,
        "original_masks_parent_state_and_finite_prescription": True,
        "two_precise_successor_status_recoveries_not_historical_formula_endorsement": True,
        "no_frozen_scientific_or_test_bytes_edited": True,
        "prepared_known_Ward_input_not_missing_scalar_operator": True,
        "original_V_G_B_and_P8_open": True,
        "unique_matching_identifiers": len({row["id"] for row in matching()})
        == len(matching()),
    }


def observable():
    return {
        "current": "The correctly oriented complete tracefree spatial Gaussian current obeys2e95 M Z136, with original homogeneous-anchored regulator error3e54 M Z136/K.",
        "canonical": "Both external factors give8e-705 and12e-746/K; these are weak derivative-losing tracefree bounds, not a reduced scalar inverse or stable background.",
        "Ward_input": "The bounded known tracefree-plus-Ward functional now uses this corrected current and remains below1e118 V04 U138. The conditional full reconstruction still has three missing ordered scalar kernels.",
        "historical_boundary": "S212's old finite value and S213's old physical formula are not re-endorsed. They are superseded by S216's corrected finite term and this complete corrected decomposition.",
        "remaining": "Full scalar anchor/state/time/contact bounds, reduced inverse, nonlinear quantum background, remaining parent/loop/heavy/cutoff and canonical vacuum/finite-gravity obligations remain open.",
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
            (
                "old_withdrawn_statuses_not_a_corrected_report",
                validate_scope,
                (frontier(), previous.matching() + [dict(ITEM)]),
            ),
        )
    )
    for value in bad:
        out.append((f"cutoff_type_{len(out)}", require_cutoff, (value,)))
        out.append((f"cutoff_mass_type_{len(out)}", require_cutoff, (2000, value)))
    for args in ((0, 1000), (1999, 1000), (2000, 999), (2000, 0), (2000, -1)):
        out.append((f"cutoff_scope_{len(out)}", require_cutoff, args))
    return out


@cache
def rejected_inputs():
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (ValueError, TypeError):
            continue
        raise ValueError("Unsupported corrected current input accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "literal_Fock_real_covariance_and_unequal_time_dependent_modes": True,
        "all_six_endpoints_and_nondeleted_sixth_bulk": True,
        "wrong_phase_omitted_bulk_and_frozen_finite_negative_controls": True,
        "unaveraged_odd_log_shell_and_original_cutoff": True,
        "both_complex_dimension_phase_branches_without_parameter_conjugation": True,
        "actual_homogeneous_anchor_and_full_chart_Ward_scope": True,
        "two_explicit_successor_status_recoveries_and_frozen_history": True,
        "original_P8_not_closed": True,
    }
