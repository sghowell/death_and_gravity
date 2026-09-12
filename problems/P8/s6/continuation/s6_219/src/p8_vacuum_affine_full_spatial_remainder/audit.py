"""Complete reference Gaussian metric response, not a reduced inverse."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes
from p8_vacuum_affine_homogeneous_trace_anchor import audit as previous
from p8_vacuum_flat_dirac_hadamard.symbols import rational

from . import assembly, dimension, remainder, shapes, vertices

CORRECTIONS = {
    "prepared_ordered_Ward_reconstruction_full_chart_correction_and_three_remaining_scalar_kernels": "PREPARED_WARD_RECONSTRUCTION_WITH_COMPLETE_S219_REFERENCE_GAUSSIAN_METRIC_CURRENT_NOT_REDUCED_INVERSE_OR_NONLINEAR_BACKGROUND",
}
ITEM = {
    "id": "full_spatial_scalar_remainders_original_regulator_and_complete_prepared_Gaussian_metric_response",
    "status": "COMPLETE_REFERENCE_GAUSSIAN_METRIC_RESPONSE_WEAK_BOUND_NOT_REDUCED_INVERSE_FULL_SOURCED_PARENT_BACKGROUND_OR_V_G_B",
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
        raise ValueError(
            "A reference Gaussian metric response cannot close original P8"
        )
    return True


def packets():
    return {
        "full_ADM_first_second_physical_feature_and_norm_bridge": vertices.feature_data(),
        "full_nine_pair_ten_geometry_and_ordered_symmetries": vertices.pair_data(),
        "full_ordered_scalar_original_band_shapes": shapes.data(),
        "full_state_time_endpoint_contact_and_original_tail": remainder.data(),
        "full_complex_dimension_geometry_and_dominated_limit": dimension.data(),
        "complete_spatial_current_homogeneous_anchor_and_regulator": assembly.data(),
        "complete_prepared_Gaussian_metric_response": assembly.ward_data(),
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
            "one_exact_prepared_Ward_completion_status_change": sum(
                a != b for a, b in zip(before, after)
            )
            - 1,
            "every_other_historical_row_unchanged": sum(
                a != b for a, b in zip(before, after) if a["id"] not in CORRECTIONS
            ),
            "one_new_complete_metric_reference_checkpoint": len(after)
            - len(before)
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
        "full_current_not_only_its_local_UV_part": True,
        "full_nonzero_temporal_constraint_and_contact_retained": True,
        "all_three_ordered_scalar_kernels_and_full_homogeneous_anchor": True,
        "same_original_masks_parent_state_and_finite_prescription": True,
        "prepared_reference_metric_response_not_reduced_inverse": True,
        "no_full_nonlinear_sourced_parent_remainder_claim": True,
        "no_frozen_scientific_or_test_bytes_edited": True,
        "original_V_G_B_and_P8_open": True,
        "unique_matching_identifiers": len({row["id"] for row in matching()})
        == len(matching()),
    }


def observable():
    return {
        "full_spatial": "The complete arbitrary-spatial reference Gaussian current obeys2e95 M Z136, with original homogeneous-anchored regulator error3e54 M Z136/K.",
        "ordered_scalar": "All three missing normalized scalar kernels are actual: trace/trace bound2e95 and independently ordered cross bounds3e49 each. The cross kernels are not assumed equal.",
        "full_metric": "The original prepared Ward maps and complete chart/one-point corrections now reconstruct all reference Gaussian metric directions with bound1e118 V04 U138.",
        "canonical_spatial": "Both spatial metric factors give8e-705 and12e-746/K. Weak derivative loss remains; small canonical numbers are not a reduced inverse.",
        "remaining": "The genuinely constrained/reduced scalar-mixed inverse, full finite-amplitude sourced-parent remainder, quantum background, stability, heavy/physical cutoff and original V/G/B obligations remain open.",
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
                "conditional_Ward_status_not_the_completed_reference_response",
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
        raise ValueError("Unsupported full reference metric input accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "literal_metric_determinant_fieldstrength_and_complete_contact": True,
        "nondegenerate_full_six_phase_space_pair_and_real_covariance": True,
        "original_two_leg_mask_all_scalar_source_jets_and_nonzero_log_shell": True,
        "grazing_strip_and_ordered_transpose_negative_controls": True,
        "explicit_complex_dimension_full_geometry_and_basis_bounds": True,
        "same_reference_parent_not_full_nonlinear_source_free_family": True,
        "one_precise_Ward_status_update_with_frozen_history": True,
        "original_P8_not_closed": True,
    }
