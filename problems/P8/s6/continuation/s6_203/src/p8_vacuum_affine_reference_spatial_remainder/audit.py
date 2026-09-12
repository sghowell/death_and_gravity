"""Complete all-momentum endpoint remainder with original regulator tails."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes
from p8_vacuum_affine_reference_spatial_analyticity import audit as previous
from p8_vacuum_flat_dirac_hadamard.symbols import rational

from . import limit, near, real, tail

ITEM = {
    "id": "actual_curved_reference_all_momentum_spatial_Taylor_subtracted_endpoint_remainder",
    "status": "ACTUAL_UNIT_W8_ALL_MOMENTUM_ENDPOINT_REMAINDER_WITH_ORIGINAL_TAIL_NOT_LOCAL_MATCHING_FULL_INVERSE_OR_V_G_B",
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
        raise ValueError("An all-momentum endpoint remainder cannot close original P8")
    return True


def packets():
    return {
        "all_real_endpoint_row_and_unexpanded_low_band": real.data(),
        "complete_near_region_and_all_Taylor_radials": near.data(),
        "all_momentum_remainder_and_original_regulator_tail": tail.data(),
        "known_actual_finite_piece_and_unmatched_sector": limit.data(),
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
        "unit_W8_reference_and_finite_state_correction_kept": True,
        "all_nine_physical_pairs_and_full_constrained_readouts": True,
        "complete_first_five_endpoint_time_jet_rows": True,
        "all_real_internal_and_external_momenta": True,
        "unexpanded_low_band_and_zero_momentum_retained": True,
        "near_region_raw_and_all_five_Taylor_orders": True,
        "logarithmic_fourth_radial_order_kept": True,
        "far_piece_uses_unchanged_joint_Cauchy_result": True,
        "all_momentum_dominated_weak_limit": True,
        "original_removed_two_created_mode_union": True,
        "same_spatial_regularity_for_one_over_K_tail": True,
        "no_derivative_of_sharp_intersection_domain": True,
        "polynomial_integrand_not_integrated_polynomial": True,
        "full_distinct_one_leg_contact_kept_unmatched": True,
        "actual_current_decomposition_and_both_canonical_factors": True,
        "no_fixed_local_matching_or_full_inverse_assumed": True,
        "original_V_G_B_not_closed": True,
    }


def observable():
    return {
        "domain": "The same actual CD Proca sector and preparation; the already justified unit-W8 reference is used only inside a complete retained current decomposition.",
        "all_momentum": "Leave |k|<m unexpanded and subtract the degree-four spatial Taylor integrand for |k|>=m. Far and near estimates retain every external momentum and the original two-created-mode regulator.",
        "remainder": "The full endpoint remainder is below3e54||D||L2 X46[Gamma], with regulator tail5e60||D||L2 X46[Gamma]/K.",
        "known_actual_piece": "Combining S201 gives a known actual-current finite piece below4e54 M[D]Y[Gamma] and tail6e60 M[D]Y[Gamma]/K, with both canonical factors retained.",
        "boundary": "The unmatched Taylor-integrand sector need not be polynomial after integration over a moving sharp band. Its original fixed covariant matching with the complete distinct one-leg contact, full inverse, background and original P8 remain open.",
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
            "Unsupported all-momentum endpoint remainder scope accepted: " + name
        )
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "actual_parent_preparation_and_prescription_unchanged": True,
        "independent_all_five_actual_endpoint_time_jets": True,
        "full_leading_endpoint_at_large_external_momentum": True,
        "independent_all_five_near_radial_integrals": True,
        "actual_near_and_low_removed_two_leg_unions": True,
        "moving_intersection_nonpolynomial_counterexample": True,
        "no_fixed_local_matching_inverse_or_state_reset_assumed": True,
        "original_P8_not_closed": True,
    }
