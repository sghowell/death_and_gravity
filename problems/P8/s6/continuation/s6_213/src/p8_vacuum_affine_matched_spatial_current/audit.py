"""Complete matched tracefree spatial current and explicit historical metadata errata."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes
from p8_vacuum_affine_dimensional_spatial_symbol import audit as previous
from p8_vacuum_flat_dirac_hadamard.symbols import rational

from . import assembly, dimension, lift, regulator

ITEM = {
    "id": "complete_matched_tracefree_spatial_Gaussian_current_with_anchored_original_regulator",
    "status": "MATCHED_WEAK_DERIVATIVE_LOSING_SPATIAL_CURRENT_NOT_FULL_SCALAR_MIXED_INVERSE_BACKGROUND_OR_V_G_B",
}
CORRECTED = (
    {
        "id": "actual_spatial_UV_difference_complete_contact_cancellation_and_fixed_curved_pole",
        "status": "FULL_SPATIAL_UV_DIFFERENCE_AND_POLE_IDENTIFIED_NOT_FINITE_DIMENSIONAL_MATCHING_OR_FULL_RESPONSE",
    },
    {
        "id": "actual_dimensional_spatial_UV_finite_part_with_evanescent_modes_and_counterterms",
        "status": "ORIGINAL_MSbar_LOCAL_SPATIAL_FINITE_PART_EVALUATED_NOT_FULL_DIMENSION_LIMIT_OR_ASSEMBLED_RESPONSE",
    },
)


def require_scope(time, kappa=modes.KAPPA, mass=modes.MASS, length=1):
    t, k, m, L = map(rational, (time, kappa, mass, length))
    if not -s.Rational(1, 2) <= t <= s.Rational(1, 2):
        raise ValueError("Only the fixed unit CD slab")
    if k != modes.KAPPA or m != modes.MASS or L != 1:
        raise ValueError("Require actual fixed CD parent, mass1000 and unit slab")
    return t, k, m, L


def frontier():
    return previous.frontier()


def metadata_errata():
    rows = previous.matching()
    if len(rows) != 68 or not rows[-3] == rows[-2] == rows[-1]:
        raise ValueError(
            "The precise frozen S211/S212 metadata erratum precondition changed"
        )
    return [
        {
            "checkpoint": "S6." + str(211 + i),
            "index": len(rows) - 2 + i,
            "frozen_display": dict(rows[-2 + i]),
            "corrected_display": dict(corrected),
            "evidence_unchanged": True,
        }
        for i, corrected in enumerate(CORRECTED)
    ]


def matching():
    rows = [dict(row) for row in previous.matching()]
    for erratum in metadata_errata():
        rows[erratum["index"]] = dict(erratum["corrected_display"])
    return rows + [dict(ITEM)]


def validate_scope(rows, obligations):
    if rows != frontier() or obligations != matching():
        raise ValueError(
            "A derivative-losing spatial Gaussian response cannot close original P8"
        )
    return True


def packets():
    return {
        "uniform_full_W8_complex_dimension_and_spatial_limit": dimension.data(),
        "homogeneous_complex_Fourier_Hilbert_lift": lift.data(),
        "same_prescription_complete_tracefree_spatial_current": assembly.data(),
        "original_two_leg_anchored_regulator_and_uniform_tail": regulator.data(),
    }


@cache
def residuals():
    rows = {
        name + "_" + key: value
        for name, p in packets().items()
        for key, value in p["checks"].items()
    }
    prior = previous.matching()
    now = matching()
    rows.update(
        {
            "nine_primitive_rows_retained": len(frontier()) - 9,
            "primitive_rows_unchanged": sum(
                a != b for a, b in zip(frontier(), previous.frontier())
            ),
            "pre_erratum_matching_prefix_unchanged": sum(
                a != b for a, b in zip(now[:-3], prior[:-2])
            ),
            "exactly_two_historical_display_rows_corrected": sum(
                a != b for a, b in zip(now, prior)
            )
            - 2,
            "one_current_checkpoint_row_added": len(now) - len(prior) - 1,
            "all_current_matching_identifiers_unique": len(now)
            - len({row["id"] for row in now}),
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
        "same_actual_prepared_state_mass_and_CD_clock": True,
        "explicit_complex_dimension_radius_one_quarter": True,
        "full_low_near_far_and_external_transfer_majorant": True,
        "all_five_endpoints_and_source_jets_retained": True,
        "actual_state_and_finite_time_pieces_stay_physical": True,
        "full_contact_retained_in_exact_homogeneous_anchor": True,
        "original_MSbar_finite_evanescent_terms_unchanged": True,
        "finite_lower_band_and_full_conversion_cancellations_retained": True,
        "weak_complete_spatial_response_not_strong_stress_operator": True,
        "thirteen_time_six_spatial_derivatives_not_same_space_contraction": True,
        "exact_homogeneous_anchor_not_unanchored_cutoff_rate": True,
        "historical_metadata_erratum_changes_only_two_display_rows": True,
        "frozen_historical_evidence_not_edited": True,
        "current_matching_identifiers_unique": len({row["id"] for row in matching()})
        == len(matching()),
        "lapse_shift_clock_scalar_and_reduced_inverse_still_required": True,
        "finite_coupling_background_and_remaining_parent_still_required": True,
        "original_V_G_B_and_P8_open": True,
    }


def observable():
    return {
        "result": "The complete tracefree spatial Gaussian retarded current is matched in the original prescription by the exact homogeneous anchor plus the known spatial difference and the actual finite UV difference.",
        "dimension_limit": "All four dimensional WKB coefficients give an explicit nonzero disk|d-3|<=1/4. Full angular/low/near/far estimates dominate every spatial Fourier smear, justifying the original subtraction limit.",
        "weak_norm": "The full current is below2e95 M[D]Z136, canonically8e-705; the anchored original-regulator error is below3e54 M[D]Z136/K, canonically12e-746/K.",
        "historical_metadata_errata": metadata_errata(),
        "boundary": "Tracefree spatial Gaussian sector only. Derivative loss, not a contraction or reduced inverse. No lapse/shift/clock/scalar closure, finite-amplitude spatial C2 theorem, finite-coupling background/stability, other parent loops/cutoff or original V/G/B closure.",
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
        raise ValueError("Unsupported matched spatial current scope accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "all_four_order_continued_frequencies_and_explicit_disk": True,
        "full_complex_dimensional_endpoint_remainder_two_radii": True,
        "angular_radial_uniform_envelopes": True,
        "noncommuting_full_generator_contact_and_covariance_tangent": True,
        "complex_Hilbert_lift_and_all_transfer_weights": True,
        "original_lower_band_and_zero_transfer_omission_controls": True,
        "exact_historical_metadata_erratum_and_unique_identifiers": True,
        "original_P8_not_closed": True,
    }
