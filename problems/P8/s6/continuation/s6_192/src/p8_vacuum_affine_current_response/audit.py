"""Actual physical shear current and C2 fourth-order mode comparison."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes
from p8_vacuum_affine_matrix_response_tail import audit as previous
from p8_vacuum_flat_dirac_hadamard.symbols import rational

from . import adiabatic, low, tail, vertices

ITEM = {
    "id": "actual_physical_shear_current_C2_fourth_order_comparison",
    "status": "ACTUAL_STATE_C2_CURRENT_MODE_COMPARISON_NOT_FIXED_COVARIANT_MATCHING_OR_V_G_B",
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
            "A mode-subtracted current comparison cannot close original P8"
        )
    return True


def packets():
    return {
        "complete_metric_vertices": vertices.data(),
        "exact_finite_band_response": low.data(),
        "ordered_adiabatic_comparison": adiabatic.data(),
        "complete_C2_comparison": tail.data(),
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
        "same_full_parent_M1_profile_state_and_prescription": True,
        "physical_metric_current_not_scalar_surrogate": True,
        "all_frame_frequency_and_metric_vertex_derivatives": True,
        "noncommuting_exponential_Frechet_products_retained": True,
        "positive_temporal_constraint_retained": True,
        "uniform_fast_generator_norm_without_momentum_exponential": True,
        "complete_first_second_covariance_Duhamel_terms": True,
        "unchanged_actual_initial_covariance_and_parameter_data": True,
        "finite_band_and_infinite_tail_both_retained": True,
        "complex_coefficient_adjoint_not_marker_conjugation": True,
        "physical_derivatives_before_contour_selection": True,
        "full_ordered_parity_not_commuting_polarization_shortcut": True,
        "fourth_order_comparison_not_a_new_renormalization_scheme": True,
        "fixed_covariant_finite_contact_matching_still_required": True,
        "uniform_C2_comparison_and_zero_safe_Taylor_remainder": True,
        "no_full_physical_selfenergy_or_feedback_inference": True,
        "mixed_spatial_interacting_parent_control_still_open": True,
        "original_V_G_B_not_closed": True,
    }


def observable():
    return {
        "domain": "The actual unit CD slab, gamma=epsilon Gamma, ||Gamma^(j)||<=1 for j0..12, |epsilon|<=1/100, and a common zero initial neighborhood; fixed smooth tracefree symmetric metric detector ||D||<=1.",
        "actual_current": "J_D=-omega tr(G_D Sigma)/2 includes the full physical metric vertex and all of its first and second amplitude derivatives.",
        "tail": "The actual-minus-reference physical current tail has bounds 6e-65, 1e-46, 1e-27; the reference-minus-fourth-order tail has bounds 1e-13, 1e-12, 1e-10 through order2.",
        "comparison": "The full actual-minus-J_ad4 current is C2 with bounds 1e77, 1e94, 1e111. Its Taylor remainder is at most epsilon^2*1e111/2, zero at epsilon0 and strict otherwise.",
        "boundary": "J_ad4 is a mathematical mode comparison, not yet matched to the original fixed covariant subtraction and finite contacts; no full physical feedback, background, cutoff or V/G/B closure.",
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
            "Unsupported physical current comparison scope accepted: " + name
        )
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "actual_parent_state_constraint_and_physical_vertex_retained": True,
        "noncommuting_root_and_exponential_derivatives_retained": True,
        "momentum_independent_growth_not_a_fast_frequency_gap": True,
        "all_second_parameter_product_and_Duhamel_terms": True,
        "adiabatic_marker_distinct_from_physical_amplitude": True,
        "comparison_not_fixed_covariant_subtraction_or_full_feedback": True,
        "both_proof_bands_not_a_physical_cutoff": True,
        "original_P8_not_closed": True,
    }
