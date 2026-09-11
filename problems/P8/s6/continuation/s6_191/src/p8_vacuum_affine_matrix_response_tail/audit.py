"""Actual C2 covariance-tail response and finite-amplitude remainder."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes
from p8_vacuum_affine_matrix_adiabatic import audit as previous
from p8_vacuum_flat_dirac_hadamard.symbols import rational

from . import covariance, mixed, reference, variation

ITEM = {
    "id": "actual_matrix_covariance_tail_C2_response_and_Taylor_remainder",
    "status": "ACTUAL_UNCHANGED_STATE_C2_HIGH_BAND_COVARIANCE_TAIL_NOT_COMPLETE_COVARIANT_RESPONSE_OR_V_G_B",
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
        raise ValueError("A C2 covariance-tail response cannot close original P8")
    return True


def packets():
    return {
        "mixed_frame_jets": mixed.data(),
        "ordered_reference_derivatives": reference.data(),
        "actual_graph_response": variation.data(),
        "complete_weighted_C2_tail": covariance.data(),
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
        "same_full_parent_M1_profile_and_prescription": True,
        "same_actual_all_order_state_and_nonzero_initial_mismatch": True,
        "only_initial_parameter_mismatches_zero": True,
        "full_three_mode_constraint_and_mixing_retained": True,
        "noncommuting_mixed_parameter_time_jets_retained": True,
        "finite_C12_direction_norm_and_amplitude_domain_stated": True,
        "no_positive_root_condition_number_assumption": True,
        "same_finite_tenth_order_reference": True,
        "complete_first_second_actual_error_equations": True,
        "quadratic_first_response_error_retained": True,
        "all_covariance_adjoint_and_product_derivatives": True,
        "varying_energy_weight_differentiated": True,
        "fixed_analysis_band_not_physical_cutoff": True,
        "uniform_infinite_tail_majorants_through_second_order": True,
        "dominated_C2_limit_and_actual_Taylor_remainder": True,
        "complete_covariant_subtraction_contacts_still_required": True,
        "full_low_band_mixed_nonlinear_parent_control_still_open": True,
        "original_V_G_B_not_closed": True,
    }


def observable():
    return {
        "domain": "The actual unit CD slab, gamma=epsilon Gamma, ||Gamma^(j)||<=1 for j0..12, |epsilon|<=1/100, and a common zero neighborhood of the initial Cauchy surface.",
        "state": "The actual all-order state and its nonzero finite-reference mismatch remain; only initial parameter derivatives of that mismatch vanish.",
        "bound": "The full balanced covariance error derivatives have coefficients2e31,4e33,2e36 and inverse-frequency powers10,9,8 in the fixed nu_minus>=1e16 high band.",
        "integral": "The complete energy-weighted tail and its first two derivatives have infinite integral bounds1e-65,1e-47,1e-28; the tail integral's finite-amplitude Taylor remainder is at most epsilon^2*1e-28/2 (zero at epsilon0 and strict otherwise).",
        "boundary": "Not complete covariant stress subtraction/contacts, low-band or full physical response, gravitational feedback inverse, interacting quantum background, cutoff or V/G/B closure.",
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
        raise ValueError("Unsupported C2 covariance-tail scope accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "actual_all_order_state_and_initial_mismatch_retained": True,
        "temporal_constraint_and_polarization_mixing_retained": True,
        "no_matrix_commutation_or_condition_number_shortcut": True,
        "full_first_and_second_parameter_errors_retained": True,
        "covariance_and_frequency_weight_both_differentiated": True,
        "finite_amplitude_tail_not_full_feedback_remainder": True,
        "proof_partition_not_a_physical_cutoff": True,
        "original_P8_not_closed": True,
    }
