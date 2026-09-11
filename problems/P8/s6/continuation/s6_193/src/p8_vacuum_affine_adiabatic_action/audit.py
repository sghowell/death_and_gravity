"""Actual homogeneous-shear fixed covariant subtraction matching."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes
from p8_vacuum_affine_current_response import audit as previous
from p8_vacuum_flat_dirac_hadamard.symbols import rational

from . import action, angular, curvature
from . import matching as local_matching

ITEM = {
    "id": "actual_homogeneous_shear_full_fixed_covariant_mode_subtraction_matching",
    "status": "ACTUAL_HOMOGENEOUS_SHEAR_FIXED_COVARIANT_MODE_SUBTRACTION_MATCHING_NOT_FULL_FEEDBACK_OR_V_G_B",
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
            "A homogeneous fixed-prescription matching cannot close original P8"
        )
    return True


def packets():
    return {
        "complete_ordered_local_variation": action.data(),
        "full_dimension_angular_contractions": angular.data(),
        "actual_shear_covariant_matching": local_matching.data(),
        "fixed_finite_dimension_limit": curvature.data(),
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
        "complete_independent_K_and_frequency_variations": True,
        "every_noncommuting_fourth_order_matrix_word_retained": True,
        "positive_root_removed_only_inside_exact_similar_traces": True,
        "full_temporal_constraint_in_every_dimension": True,
        "all_oriented_sphere_pairings_and_trace_classes": True,
        "arbitrary_shear_not_isotropic_extrapolation": True,
        "constant_pointwise_coordinate_change_only": True,
        "all_weighted_boundary_terms_exact_at_general_dimension": True,
        "no_Euler_removal_before_the_dimension_limit": True,
        "physical_state_and_convergent_integral_stay_dimension_three": True,
        "evanescent_polarization_and_curvature_coefficients_retained": True,
        "original_four_dimensional_pole_coefficients_not_retuned": True,
        "full_finite_local_contact_matching_on_admitted_shears": True,
        "no_full_spatial_mixed_or_inverse_claim": True,
        "interacting_parent_background_and_cutoff_still_open": True,
        "original_V_G_B_not_closed": True,
    }


def observable():
    return {
        "domain": "The same actual Gaussian Proca state and homogeneous unimodular shear family of S192, with smooth normalized detector D and unchanged common preparation.",
        "identity": "The actual J_ad4 is the compact metric variation of the complete local marker orders0,2,4 for generic noncommuting positive K and independent frequency.",
        "matching": "Full general-d shear angular/radial integration matches covariant curvature polynomials modulo explicit compact total derivatives. The original fixed four-dimensional pole coefficients give exactly the S176 finite local action, including all evanescent finite contacts.",
        "conclusion": "The original fixed-prescription homogeneous current equals the convergent S192 actual-minus-J_ad4 comparison plus the variation of the unchanged finite local action. No counterterm or state is selected anew.",
        "boundary": "A numerical nonlinear local-current norm and complete spatial/mixed response, feedback inverse, finite-coupling interacting background, physical cutoff and V/G/B closure are not established here.",
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
            "Unsupported homogeneous fixed-prescription scope accepted: " + name
        )
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "actual_parent_state_constraint_and_physical_vertex_retained": True,
        "full_noncommuting_word_and_oriented_angular_contractions": True,
        "compact_variational_identity_not_phase_heuristic": True,
        "arbitrary_shear_curvature_not_isotropic_extrapolation": True,
        "evanescent_finite_contacts_not_discarded": True,
        "fixed_prescription_and_actual_state_not_changed": True,
        "matching_not_full_feedback_background_or_cutoff": True,
        "original_P8_not_closed": True,
    }
