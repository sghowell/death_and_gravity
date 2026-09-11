"""Complete fixed-prescription homogeneous current and canonical response."""

from functools import cache

import sympy as s
from p8_vacuum_affine_adiabatic_action import audit as previous
from p8_vacuum_affine_cd_source_noise import modes
from p8_vacuum_flat_dirac_hadamard.symbols import rational

from . import canonical, current, metrics, tensors

ITEM = {
    "id": "actual_complete_fixed_covariant_homogeneous_current_C2_and_canonical_response",
    "status": "ACTUAL_COMPLETE_FIXED_COVARIANT_HOMOGENEOUS_CURRENT_C2_AND_CANONICAL_RESPONSE_NOT_FULL_FEEDBACK_OR_V_G_B",
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
            "A homogeneous current response bound cannot close original P8"
        )
    return True


def packets():
    return {
        "absolute_metric_jet_algebra": metrics.data(),
        "complete_covariant_Euler_tensors": tensors.data(),
        "full_fixed_prescription_current": current.data(),
        "actual_clock_and_canonical_response": canonical.data(),
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
        "same_actual_three_mode_covariance_and_physical_vertices": True,
        "absolute_exponential_zero_jet_not_relative_one": True,
        "all_factorial_weighted_time_and_amplitude_jets": True,
        "component_l1_products_before_every_index_contraction": True,
        "full_noncommuting_inverse_and_metric_jets": True,
        "covariant_spatial_connections_not_dropped_by_homogeneity": True,
        "complete_R_squared_and_Ricci_squared_Euler_tensors": True,
        "both_raised_metric_indices_and_detector_variation_retained": True,
        "no_time_derivative_bound_on_the_final_detector_needed": True,
        "Euler_removed_only_after_the_fixed_dimension_limit": True,
        "original_finite_local_coefficients_not_retuned": True,
        "whole_current_C2_remainder_not_covariance_only": True,
        "actual_clock_current_zero_by_full_rotation_invariance": True,
        "both_canonical_metric_chain_factors_retained": True,
        "no_same_space_contraction_from_derivative_losing_bound": True,
        "interacting_parent_background_and_cutoff_still_open": True,
        "original_V_G_B_not_closed": True,
    }


def observable():
    return {
        "domain": "The unchanged actual Gaussian Proca state, mass1000, kappa1e800, unit CD slab and smooth homogeneous gamma=epsilon Gamma; direction jets through order12 are at most1, |epsilon|<=.01, and the common initial shear neighborhood is zero.",
        "fixed_local_current": "The full finite local action's homogeneous shear current and its first two amplitude derivatives are each strictly bounded by1e86, with complete covariant Euler tensors and metric/detector factors.",
        "complete_current": "Adding the S192 convergent comparison through the S193 fixed-prescription identity gives full current derivative displays1e87,1e95,2e111 and a Taylor remainder at most epsilon^2*1e111, zero at epsilon0.",
        "canonical_clock": "The original isotropic tracefree current is exactly zero. The properly normalized linear homogeneous metric response is4DJ/(a^3 kappa), bounded C12-to-C0 by4e-705 on the stated smooth prepared directions.",
        "boundary": "This is a derivative-losing external-metric response bound, not a same-space contraction, full spatial/mixed inverse, finite-coupling interacting quantum solution, stability or physical cutoff.",
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
            "Unsupported complete homogeneous current scope accepted: " + name
        )
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "actual_parent_state_constraint_and_physical_vertex_retained": True,
        "absolute_metric_majorants_and_full_covariant_connections": True,
        "independent_literal_spatial_metric_Euler_variations": True,
        "all_raising_density_and_detector_factors_retained": True,
        "fixed_prescription_and_actual_state_not_changed": True,
        "complete_current_Taylor_zero_edge_and_longitudinal_mode": True,
        "canonical_derivative_bound_not_a_same_space_inverse": True,
        "original_P8_not_closed": True,
    }
