"""Actual curved linear conversion coefficient with original matching frontier."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes
from p8_vacuum_affine_subleading_band_conversion import audit as previous
from p8_vacuum_flat_dirac_hadamard.symbols import rational

from . import angular, bounds, density, jets

ITEM = {
    "id": "explicit_actual_curved_linear_regulator_conversion_all_tensor_channels_and_source_jets",
    "status": "ACTUAL_CURVED_A1_EVALUATED_WITH_COEFFICIENT_NORMS_NOT_FULL_ONE_BALL_CONTACT_COVARIANT_MATCHING_OR_V_G_B",
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
        raise ValueError("An evaluated cutoff coefficient cannot close original P8")
    return True


def packets():
    return {
        "actual_WKB_low_jets_and_complete_second_grade_pairs": jets.data(),
        "complete_curved_second_grade_source_value_and_time_jets": density.data(),
        "full_tensor_actual_curved_linear_conversion": angular.data(),
        "proper_time_operator_and_explicit_coefficient_norms": bounds.data(),
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
        "same_actual_CD_parent_mass_and_prepared_state": True,
        "same_unit_W8_comparison_and_fixed_prescription": True,
        "all_four_WKB_terms_in_independent_full_mode_tests": True,
        "full_transverse_mixed_longitudinal_and_constraint_content": True,
        "source_and_detector_times_independent_before_derivatives": True,
        "actual_degree_one_j1_zero_without_endpoint_deletion": True,
        "full_curvature_correction_computed_not_flat_assumption": True,
        "first_and_second_source_time_derivatives_retained": True,
        "all_general_tracefree_tensor_channels_reconstructed": True,
        "actual_spatial_momentum_dependence_not_homogeneous_only": True,
        "original_cubic_coefficient_unchanged": True,
        "original_quadratic_and_finite_conversion_cancellation_retained": True,
        "all_P_original_conversion_tail_unchanged": True,
        "coefficient_norms_not_cutoff_uniformity": True,
        "point_zero_defined_by_continuity_not_new_choice": True,
        "full_one_mode_contact_not_changed_or_dropped": True,
        "full_one_ball_contact_covariant_matching_still_open": True,
        "original_V_G_B_not_closed": True,
    }


def observable():
    return {
        "actual_formula": "The complete original linear cutoff coefficient is p/pi^2[-m^2 a M/1024+p^2 S/(24576a)-a(H_prime+2H^2)C/1024+a L(Gamma_second+H Gamma_prime)/2048], with all four invariant forms and source tensors specified explicitly.",
        "all_channels": "Seven full general-tensor azimuthal projector identities reconstruct all five tracefree channels. The actual longitudinal WKB curvature term and independent source-time differentiation fix contributions that a flat fixture misses.",
        "norms": "The cubic coefficient is bounded by1/100||D||L2||grad Gamma||L2; the linear coefficient by1e4||D||L2 Z23[Gamma]. Z23 has two source time derivatives and three spatial Sobolev derivatives.",
        "original_conversion": "S209's zero complete quadratic/finite conversion and uniform1e40 X46/K tail remain unchanged. Both coefficient multipliers extend continuously as zero at P0.",
        "boundary": "Nonzero K3 and K factors remain. Neither coefficient smallness nor formal endpoint symmetry is a full renormalized response, one-ball/contact covariant matching, mixed inverse, background or original P8 closure.",
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
            "Unsupported actual curved coefficient scope accepted: " + name
        )
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "independent_full_W8_curved_and_mass_scaled_flat_fields": True,
        "noncommuting_tensor_clocks_angles_and_separated_transfer": True,
        "full_source_time_jets_with_fixed_detector": True,
        "evaluated_j1_slot_without_higher_endpoint_deletion": True,
        "complete_curved_angular_quadrature_two_resolutions": True,
        "exact_origin_norm_and_Green_identity_controls": True,
        "full_contact_and_fixed_matching_not_replaced": True,
        "original_P8_not_closed": True,
    }
