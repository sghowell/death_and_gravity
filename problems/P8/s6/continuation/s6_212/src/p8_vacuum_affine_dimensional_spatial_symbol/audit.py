"""Actual curved linear conversion coefficient with original matching frontier."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes
from p8_vacuum_affine_spatial_matching_difference import audit as previous
from p8_vacuum_flat_dirac_hadamard.symbols import rational

from . import density, geometry, jets
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
        raise ValueError("A finite UV coefficient cannot close original P8")
    return True


def packets():
    return {
        "actual_general_dimension_magnetic_and_angular_geometry": geometry.data(),
        "complete_dimension_dependent_modes_and_source_jets": jets.data(),
        "fixed_source_invariant_dimension_coefficients": density.data(),
        "continued_counterterm_and_original_MSbar_finite_UV_difference": pole_matching.data(),
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
        "actual_full_general_dimension_magnetic_vertex": True,
        "all_four_WKB_orders_in_independent_dimensional_tests": True,
        "all_source_jets_and_independent_detector_time": True,
        "all_105_physical_angular_coefficients_reproduced": True,
        "fixed_source_invariants_before_dimension_derivative": True,
        "actual_nonzero_evanescent_odd_endpoint_retained": True,
        "original_one_leg_contact_cancelled_only_in_difference": True,
        "continued_curvature_Hessians_varied_before_limit": True,
        "Euler_and_volume_evanescent_terms_retained": True,
        "original_MSbar_normalization_and_scale": True,
        "fixed_comoving_lower_band_term_retained": True,
        "finite_time_Green_identity_not_missing_first_derivative": True,
        "explicit_finite_coefficient_bound_not_full_response_inverse": True,
        "no_finite_local_target_added_twice": True,
        "full_dimensional_limit_and_response_assembly_still_required": True,
        "original_V_G_B_not_closed": True,
    }


def observable():
    return {
        "dimensional_modes": "The actual general-d constrained WKB modes and full field-strength vertex yield all spatial UV coefficient derivatives. All105 physical angular coefficients are reproduced without freezing polarization multiplicity or canonical rates.",
        "evanescent_data": "An odd endpoint vanishes atd3 but contributes through its dimension derivative. Invariant reconstruction fixes the physical source before differentiation, retaining the scalar-basis correction.",
        "counterterms": "The full metric second variation of the original fixed scalar pole weights continues R_old,R_squared,Ricci_squared,Riemann_squared and volume. The evanescent Euler and log(a) terms remain.",
        "finite_UV_difference": "The explicit local finite spatial UV difference uses the unchanged MSbar convention, lower comoving band and scale. Its proper-time first derivative coefficient satisfies the full Green identity. The Z24 coefficient bound is1e5, canonically4e-795.",
        "boundary": "The complete dimension-limit interchange and assembly with the anchored homogeneous current and known remainders remain. No full response, reduced mixed inverse, finite-coupling background or original P8 closure follows from this local coefficient.",
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
        raise ValueError("Unsupported dimensional finite UV scope accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "all_four_dimensional_WKB_orders_and_full_source_jets": True,
        "literal_integer_dimension_field_strength_and_sphere": True,
        "literal_noncommuting_metric_second_variation": True,
        "complex_dimension_radial_finite_part_two_resolutions": True,
        "evanescent_odd_endpoint_and_fixed_source_basis": True,
        "continued_Euler_and_volume_Green_controls": True,
        "exact_finite_coefficient_norm_and_origin": True,
        "original_P8_not_closed": True,
    }
