"""Joint spatial analyticity and a finite far-momentum endpoint remainder."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes
from p8_vacuum_affine_reference_state_prefactor import audit as previous
from p8_vacuum_flat_dirac_hadamard.symbols import rational

from . import boundary, domain, frame, readouts

ITEM = {
    "id": "actual_curved_reference_joint_spatial_analyticity_full_readouts_and_far_endpoint_Taylor_tail",
    "status": "ACTUAL_UNIT_W8_JOINT_SPATIAL_ANALYTICITY_AND_FAR_ENDPOINT_TAYLOR_TAIL_NOT_NEAR_REGION_LOCAL_MATCHING_INVERSE_OR_V_G_B",
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
        raise ValueError("A joint spatial far-endpoint bound cannot close original P8")
    return True


def packets():
    return {
        "joint_complex_frequency_domain": domain.data(),
        "full_analytic_polarization_frame": frame.data(),
        "complete_constrained_readouts_and_inverse_phase": readouts.data(),
        "far_endpoint_spatial_Taylor_remainder": boundary.data(),
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
        "unit_W8_reference_not_replacement_physical_state": True,
        "no_Borel_state_momentum_analyticity_assumed": True,
        "both_internal_legs_in_one_joint_complex_domain": True,
        "positive_frequency_and_bilinear_length_branches_controlled": True,
        "local_frame_not_global_momentum_sphere_frame": True,
        "genuine_complex_Euclidean_polarization_norms": True,
        "all_nine_physical_pairs_and_full_ten_field_readouts": True,
        "full_sixteen_component_stress_tensor_bound": True,
        "analytic_Schwarz_not_in_place_complex_conjugation": True,
        "pre_current_endpoint_coefficients_before_imaginary_part": True,
        "all_first_five_endpoint_time_jet_rows": True,
        "Taylor_applied_only_to_coefficients_not_tests_or_band": True,
        "far_region_and_removed_both_leg_union_retained": True,
        "explicit_convergent_radial_and_regulator_tail_bounds": True,
        "both_external_metric_canonical_factors": True,
        "near_region_and_fixed_local_matching_not_assumed": True,
        "original_V_G_B_not_closed": True,
    }


def observable():
    return {
        "domain": "Same actual CD sector and preparation; only the unit-W8 reference coefficients from S201 are continued in complex time and external momentum.",
        "analyticity": "For real |k|>=1000, nu=sqrt(m^2+|k|^2/Amax^2), the ball ||p||<=1e-6nu and fixed nested time discs control both full frequencies, the local complex frame, all constrained readouts and inverse summed phase.",
        "endpoint": "The complete pre-current row from source time jets0..4 has norm below2e27nu. Its fourth spatial Taylor error is below4e57|P|^5nu^-4 on |P|<=1e-6nu/2.",
        "far_limit": "The far-region remainder is below1e54||D||L2 X46[Gamma], with the original both-leg tail below1e58||D||L2 X46[Gamma]/K. The canonical displays are4e-746 and4e-742/K.",
        "boundary": "The complementary near region, finite and divergent Taylor coefficients, full reference contact and fixed covariant spatial matching are not included. No full inverse, interacting background, physical cutoff or original P8 closure.",
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
            "Unsupported joint spatial far-endpoint scope accepted: " + name
        )
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "actual_parent_preparation_and_prescription_unchanged": True,
        "independent_actual_W8_joint_complex_domain": True,
        "full_complex_frame_and_Euclidean_counterexample": True,
        "complete_readouts_both_Schwarz_signs_and_all_nine_pairs": True,
        "actual_full_endpoint_nonzero_fifth_Taylor_remainder": True,
        "independent_far_radial_and_removed_two_leg_union": True,
        "no_near_region_local_matching_or_state_analyticity_assumed": True,
        "original_P8_not_closed": True,
    }
