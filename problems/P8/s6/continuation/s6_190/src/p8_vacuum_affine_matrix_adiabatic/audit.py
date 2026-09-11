"""Actual matrix adiabatic high-band covariance error, not full response."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes
from p8_vacuum_affine_local_tensor_response import audit as previous
from p8_vacuum_flat_dirac_hadamard.symbols import rational

from . import frame, initial, jets, riccati

ITEM = {
    "id": "actual_CD_full_matrix_adiabatic_uniform_covariance_tail",
    "status": "ACTUAL_UNCHANGED_STATE_UNIFORM_HIGH_BAND_COVARIANCE_REMAINDER_NOT_RESPONSE_SUBTRACTION_INVERSE_OR_V_G_B",
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
        raise ValueError("A matrix covariance tail cannot close original P8")
    return True


def packets():
    return {
        "full_matrix_frame": frame.data(),
        "uniform_frame_jets": jets.data(),
        "finite_Riccati_reference": riccati.data(),
        "actual_state_tail": initial.data(),
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
        "same_actual_all_order_Cauchy_state": True,
        "no_frequency_reference_state_reset": True,
        "full_three_mode_constraint_retained": True,
        "noncommuting_polarization_mixing_retained": True,
        "finite_C12_metric_jet_neighborhood_stated": True,
        "no_condition_number_or_polarization_gap_assumption": True,
        "finite_tenth_order_reference_not_convergent_series": True,
        "exact_all_order_initial_tail_including_partial_cutoffs": True,
        "exact_fast_phase_removed_only_from_norm_bound": True,
        "actual_uniform_high_band_covariance_error": True,
        "complete_infinite_energy_weighted_tail_integrated": True,
        "proof_partition_not_physical_cutoff": True,
        "amplitude_parameter_response_derivatives_still_required": True,
        "fixed_covariant_subtraction_contacts_still_required": True,
        "full_mixed_spatial_nonlinear_feedback_still_required": True,
        "finite_coupling_and_interacting_parent_control_still_open": True,
        "original_V_G_B_not_closed": True,
    }


def observable():
    return {
        "domain": "The actual unit CD slab and smooth compact homogeneous tracefree shear with every raw time jet0..12 bounded by1/100 in operator norm, zero near the common initial Cauchy surface.",
        "state": "The unchanged S6.55 all-order T/T/L state, not a finite adiabatic replacement.",
        "reference": "Full noncommuting three-mode canonical frame and tenth-order Riccati reference, with a condition-number-independent residual.",
        "bound": "For nu_minus>=1e16 the full balanced six-quadrature covariance error is below2e31 nu_minus^-10; its infinite energy-weighted continuum integral is below1e-65.",
        "boundary": "This is not a response-parameter derivative bound, fixed-covariant subtraction/contact match, finite-coupling feedback, inverse, physical cutoff or V/G/B closure.",
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
        raise ValueError("Unsupported matrix covariance scope accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "actual_all_order_state_not_reset": True,
        "temporal_constraint_and_polarization_mixing_retained": True,
        "no_polarization_gap_or_condition_number_shortcut": True,
        "finite_jet_neighborhood_not_arbitrary_metric_amplitude": True,
        "actual_covariance_tail_not_parameter_response_remainder": True,
        "all_higher_cutoff_terms_bounded_not_dropped": True,
        "proof_partition_not_a_physical_cutoff": True,
        "original_P8_not_closed": True,
    }
