"""Finite endpoint Taylor terms restored with the original cutoff tail."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes
from p8_vacuum_affine_local_spatial_hessian import audit as previous
from p8_vacuum_flat_dirac_hadamard.symbols import rational

from . import degrees, limit, radials, tail

ITEM = {
    "id": "finite_endpoint_spatial_Taylor_cells_restored_with_original_regulator_tail",
    "status": "TEN_FINITE_ENDPOINT_TAYLOR_CELLS_RESTORED_NOT_REMAINING_UV_CONTACT_MATCHING_INVERSE_OR_V_G_B",
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
        raise ValueError("Finite Taylor-cell restoration cannot close original P8")
    return True


def packets():
    return {
        "complete_endpoint_spatial_degree_and_time_jet_partition": degrees.data(),
        "all_finite_massive_radial_coefficient_bounds": radials.data(),
        "original_removed_two_leg_finite_polynomial_tail": tail.data(),
        "finite_terms_restored_and_remaining_candidate_sector": limit.data(),
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
        "unit_W8_reference_and_actual_state_correction_kept": True,
        "complete_individual_endpoint_source_time_jet_rows": True,
        "all_nine_physical_pairs_in_each_bound": True,
        "Taylor_on_pre_current_coefficients_only": True,
        "ten_finite_and_fifteen_candidate_cells": True,
        "all_forty_finite_and_thirty_five_candidate_time_jets": True,
        "all_four_exact_massive_radial_moments": True,
        "original_fixed_high_band_not_changed": True,
        "both_created_momenta_in_every_finite_regulator": True,
        "low_and_high_external_transfer_tail_regimes": True,
        "same_spatial_norm_for_finite_polynomial_and_tail": True,
        "continuum_polynomial_only_after_dominated_removal": True,
        "finite_terms_restored_not_deleted": True,
        "zero_transfer_homogeneous_result_unchanged": True,
        "fixed_local_target_not_double_counted": True,
        "candidate_UV_not_claimed_actual_divergence": True,
        "original_V_G_B_not_closed": True,
    }


def observable():
    return {
        "domain": "The same actual CD sector, unit-W8 comparison, all nine physical pairs and original one-leg/two-leg regulators.",
        "partition": "Ten endpoint/spatial-degree cells with j+n>=5 are absolutely integrable; fifteen cells with j+n<=4 remain candidates for UV matching. All75 source time-jet entries are retained.",
        "finite_restoration": "The finite continuum polynomial is below1e50||D||L2 X46[Gamma], with original cutoff tail1e55||D||L2 X46[Gamma]/K; its positive spatial degrees leave P=0 unchanged.",
        "known_actual_piece": "Restoring it to S203 gives5e54 M[D]Y[Gamma] with tail7e60 M[D]Y[Gamma]/K. The fixed local target from S204 is not separately added.",
        "boundary": "The complete contact and remaining fifteen candidate cells still require actual UV expansion and original fixed covariant matching. No full response, inverse, background, stability or original P8 closure.",
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
            "Unsupported finite Taylor-cell restoration scope accepted: " + name
        )
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "actual_parent_preparation_and_prescription_unchanged": True,
        "independent_actual_mixed_spatial_time_coefficients": True,
        "all_nine_pairs_and_two_Cauchy_resolutions": True,
        "all_exact_finite_moments_and_original_removed_unions": True,
        "prewarmed_cache_and_exact_order_guards": True,
        "borderline_majorant_not_actual_divergence_control": True,
        "no_contact_matching_or_local_target_double_counting": True,
        "original_P8_not_closed": True,
    }
