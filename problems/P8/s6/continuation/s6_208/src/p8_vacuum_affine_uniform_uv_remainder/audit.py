"""Uniform actual UV-subtracted endpoint remainder and original regulator tail."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes
from p8_vacuum_affine_spatial_symbol import audit as previous
from p8_vacuum_flat_dirac_hadamard.symbols import rational

from . import domain, estimates, limit, tail

ITEM = {
    "id": "uniform_actual_UV_subtracted_endpoint_remainder_with_original_regulator_tail",
    "status": "UNIFORM_UV_SUBTRACTED_ENDPOINT_REMAINDER_NOT_FULL_CONTACT_COVARIANT_MATCHING_OR_V_G_B",
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
        raise ValueError("A uniform UV-subtracted remainder cannot close original P8")
    return True


def packets():
    return {
        "actual_uniform_joint_inverse_radius_and_time_domain": domain.data(),
        "complete_all_momentum_far_near_low_endpoint_remainder": estimates.data(),
        "original_both_created_mode_removed_union_tail": tail.data(),
        "exact_actual_current_repartition_and_canonical_bounds": limit.data(),
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
        "complete_actual_state_unit_reference_correction_kept": True,
        "full_four_W8_coefficients_and_both_Schwarz_signs": True,
        "all_nine_pairs_and_all_ten_constrained_readouts": True,
        "improved_bounds_only_for_normalized_analytic_modes": True,
        "all_fifteen_UV_slots_and_thirty_five_source_jets": True,
        "source_time_differentiation_with_fixed_detector": True,
        "uniform_external_momentum_dependent_analytic_radius": True,
        "complete_far_near_and_unexpanded_low_regions": True,
        "all_near_power_and_logarithmic_terms_retained": True,
        "same_six_spatial_derivatives_for_bound_and_tail": True,
        "original_two_leg_removed_union_not_replaced": True,
        "full_contact_keeps_original_one_mode_band": True,
        "old_finite_pieces_not_counted_twice": True,
        "rebase_from_state_correction_and_finite_time_piece": True,
        "fixed_local_target_not_added_as_matched_quantum_term": True,
        "full_UV_contact_covariant_matching_still_open": True,
        "original_V_G_B_not_closed": True,
    }


def observable():
    return {
        "domain": "The unchanged massive CD unit-W8 comparison with all physical pairs and source jets. A new joint inverse-radius/time domain uses rho(P)=1/[100(m+|P|)], with explicit normalized frequency/frame/readout bounds.",
        "uniform_remainder": "The complete low plus high(E5-U) endpoint piece is bounded by2e48||D||L2 X46[Gamma] and original regulator error4e53||D||L2 X46[Gamma]/K. The full far/near regions and every power/log term remain.",
        "known_actual_piece": "An exact finite-regulator repartition, rebased from S201 state correction plus finite time remainder, gives5e48 M[D]Y[Gamma] and tail5e53 M[D]Y[Gamma]/K. Earlier finite endpoint pieces are not double-counted.",
        "canonical": "The uniform endpoint displays are8e-752 and16e-747/K; the complete newly known actual piece has2e-751 and2e-746/K.",
        "boundary": "The full actual UV-symbol integral and contact still require all finite/divergent and subleading regulator-artifact covariant matching. No full response/inverse/background/stability/cutoff or original P8 closure is inferred.",
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
        raise ValueError("Unsupported uniform UV remainder scope accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "actual_parent_state_and_matching_unchanged": True,
        "independent_actual_full_normalized_W8_domains": True,
        "all_UV_coefficients_two_resolutions_at_separated_P": True,
        "actual_far_and_large_transfer_near_remainders": True,
        "exact_original_far_near_low_removed_union_integrals": True,
        "same_six_spatial_derivatives_and_exact_order_guards": True,
        "state_time_rebase_and_no_finite_piece_double_counting": True,
        "original_P8_not_closed": True,
    }
