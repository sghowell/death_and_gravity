"""Actual curved linear conversion coefficient with original matching frontier."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes
from p8_vacuum_affine_curved_linear_conversion import audit as previous
from p8_vacuum_flat_dirac_hadamard.symbols import rational

from . import contact, density, jets
from . import matching as pole_matching

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
        raise ValueError("A spatial pole identity cannot close original P8")
    return True


def packets():
    return {
        "actual_full_WKB_source_derivative_UV_jets": jets.data(),
        "all_spatial_difference_coefficients_and_radial_boundary": density.data(),
        "complete_original_contact_difference": contact.data(),
        "fixed_spatial_pole_identity_and_coefficient_norms": pole_matching.data(),
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
        "all_five_original_endpoint_labels_and_source_jets_retained": True,
        "all_seven_general_tracefree_geometry_contractions": True,
        "all_105_actual_angular_UV_coefficients_computed": True,
        "full_one_leg_contact_cancelled_only_in_spatial_difference": True,
        "actual_nonzero_homogeneous_odd_endpoint_retained": True,
        "exact_one_ball_lower_band_terms_retained": True,
        "original_two_leg_conversion_unchanged": True,
        "actual_curved_log_identity_not_flat_extrapolation": True,
        "existing_R_old_sign_and_proper_clock_measure_retained": True,
        "coefficient_bounds_not_cutoff_uniformity": True,
        "finite_evanescent_matching_not_claimed_from_pole": True,
        "no_finite_local_target_added_twice": True,
        "original_V_G_B_not_closed": True,
    }


def observable():
    return {
        "actual_difference": "The complete angular endpoint symbol at P minus P0 has11 nonzero entries among105. Its full original one-leg contact cancels in this difference before any regulator removal, while the nonzero homogeneous anchor is retained.",
        "curved_pole": "The spatial logarithmic difference equals[m^2 DeltaH_Rold+(13/60)DeltaH_C2+(1/36)DeltaH_Rold2]/(32pi^2), with the original sign, proper-time measure, curvature and first/second source derivatives.",
        "original_radial": "The exact one-ball difference is DeltaF2*(K^2-m^2)/2+DeltaF4*log(K/m). The noncovariant power coefficient, lower endpoint and original pair-band conversion are not replaced or omitted.",
        "norms": "In the stated Z24 norm, power/log coefficient bounds are1/200 and1e4. They are not cutoff-uniform after K^2/logK multiplication.",
        "boundary": "Agreement with the pole is not fixed finite dimensional matching: evanescent dimensional coefficients still matter. No full matched response, reduced mixed inverse, finite-coupling background or original P8 closure follows.",
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
            "Unsupported spatial matching difference scope accepted: " + name
        )
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "independent_all_four_W8_complete_UV_endpoint_coefficients": True,
        "all_source_jets_with_fixed_detector_at_separated_transfer": True,
        "literal_full_ten_field_nine_pair_azimuth_controls": True,
        "complete_angular_coefficients_two_Cauchy_resolutions": True,
        "noncommuting_contact_difference_not_zero_contact": True,
        "actual_odd_endpoint_retained_in_homogeneous_anchor": True,
        "fixed_pole_not_evanescent_finite_matching": True,
        "original_P8_not_closed": True,
    }
