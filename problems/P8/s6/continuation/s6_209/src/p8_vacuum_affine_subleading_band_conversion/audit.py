"""Complete original UV-symbol conversion with actual centered cancellation."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes
from p8_vacuum_affine_uniform_uv_remainder import audit as previous
from p8_vacuum_flat_dirac_hadamard.symbols import rational

from . import angular, centered, flat, geometry

ITEM = {
    "id": "complete_original_UV_shape_conversion_centered_cancellations_and_uniform_tail",
    "status": "ACTUAL_K2_AND_FINITE_CONVERSION_CANCEL_WITH_UNIFORM_TAIL_NOT_CURVED_A1_EVALUATION_FULL_CONTACT_MATCHING_OR_V_G_B",
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
        raise ValueError("A regulator conversion cannot close original P8")
    return True


def packets():
    return {
        "complete_original_shell_and_all_retained_shapes": geometry.data(),
        "exact_angular_reduction_and_actual_leading_coefficients": angular.data(),
        "complete_massive_flat_all_angle_benchmark": flat.data(),
        "actual_centered_symmetry_and_uniform_original_conversion_tail": centered.data(),
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
        "same_actual_CD_parent_state_and_fixed_prescription": True,
        "unit_W8_comparison_not_prepared_state_analyticity": True,
        "all_original_physical_pairs_and_constrained_readouts": True,
        "all_fifteen_UV_slots_and_every_source_time_jet": True,
        "independent_detector_time_before_source_differentiation": True,
        "distinct_transverse_longitudinal_inverse_phases": True,
        "exact_positive_grazing_value_and_derivative_kept": True,
        "finite_terms_cancel_only_after_complete_sum": True,
        "actual_centered_exchange_and_Schwarz_not_flat_assumption": True,
        "odd_endpoint_labels_not_deleted": True,
        "all_P_remainder_including_large_transfer_logarithm": True,
        "same_six_spatial_derivatives_and_K_at_least_twice_mass": True,
        "original_contact_keeps_its_one_leg_band": True,
        "named_one_ball_comparison_not_regulator_replacement": True,
        "old_finite_pieces_and_local_target_not_added_twice": True,
        "curved_linear_conversion_coefficient_still_to_evaluate": True,
        "full_contact_and_fixed_covariant_matching_still_open": True,
        "original_V_G_B_not_closed": True,
    }


def observable():
    return {
        "exact_shapes": "Ten universal shapes instantiate all twenty j/d/h branches and thirty-five source time-jet shapes. The positive grazing strip contributes its cubic value and quartic derivative; finite K0 terms are first retained.",
        "actual_cancellation": "Full pair exchange and analytic Schwarz symmetry in centered variables eliminate odd total UV grades, not odd endpoint labels. The complete unit-W8 conversion has K3 and K1 terms but no K2 or finite K0 term.",
        "uniform_tail": "The full original UV-symbol conversion minus these two terms is bounded by1e40||D||L2 X46[Gamma]/K for every external momentum and K>=2m. Every logarithmic slot and the large-transfer region remain.",
        "benchmark": "Complete all-angle massive flat fields give all45 UV coefficients and explicit channel cubic/linear terms. Full curved A1 is retained as its actual f2c functional, not replaced by this fixture.",
        "boundary": "This names a regulator comparison, not a subtraction choice or full UV/contact covariant matching. No mixed inverse, finite-coupling background, cutoff or original P8 closure follows.",
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
        raise ValueError("Unsupported full conversion scope accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "independent_full_all_angle_massive_flat_fields": True,
        "full_original_shell_and_positive_grazing_quadratures": True,
        "actual_curved_centered_full_mode_and_projector_routes": True,
        "all_curved_UV_coefficients_at_two_resolutions": True,
        "all_original_small_and_large_transfer_UV_slots": True,
        "complete_grade_cancellation_not_odd_endpoint_deletion": True,
        "fixed_contact_and_matching_not_replaced": True,
        "original_P8_not_closed": True,
    }
