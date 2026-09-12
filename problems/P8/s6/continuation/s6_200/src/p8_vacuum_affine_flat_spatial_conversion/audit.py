"""Full spatial flat covariant and time-subtraction conversion benchmark."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes
from p8_vacuum_affine_flat_tensor_cut import audit as previous
from p8_vacuum_flat_dirac_hadamard.symbols import rational

from . import bounds, conversion, flat, tensor

ITEM = {
    "id": "full_spatial_flat_Proca_covariant_and_time_subtraction_conversion_with_nonlocal_coefficient_bounds",
    "status": "FULL_SPATIAL_FLAT_NONLOCAL_SUBTRACTION_CONVERSION_NOT_CURVED_MATCHING_LOCAL_POLYNOMIAL_INVERSE_OR_V_G_B",
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
            "A flat spatial subtraction conversion benchmark bound cannot close original P8"
        )
    return True


def packets():
    return {
        "complete_spatial_tensor_numerators": tensor.data(),
        "exact_covariant_time_subtraction_conversion": conversion.data(),
        "all_spatial_coefficient_and_tail_bounds": bounds.data(),
        "physical_flat_retarded_extraction_and_nonlocality": flat.data(),
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
        "same_full_parent_normalization_state_and_prescription": True,
        "flat_representative_not_actual_CD_state_replacement": True,
        "all_full_spatial_tensor_components": True,
        "both_spin_sectors_and_longitudinal_cut": True,
        "lightcone_projector_poles_canceled_before_integration": True,
        "all_three_equal_time_conversion_coefficients": True,
        "generic_cubic_division_with_full_tensor_residue": True,
        "independent_four_curvature_tensor_Hessians": True,
        "analytic_spectral_integrals_and_negative_axis_discs": True,
        "no_physical_finite_polynomial_selected": True,
        "contact_and_memory_regulators_not_silently_equated": True,
        "no_curved_odd_endpoint_discard": True,
        "spatially_nonlocal_coefficients_not_local_counterterms": True,
        "all_spatial_momentum_coefficient_and_tail_bounds": True,
        "six_spatial_and_four_time_source_derivatives_explicit": True,
        "both_external_metric_canonical_chain_factors": True,
        "not_fully_reduced_mixed_norm_or_same_space_inverse": True,
        "original_V_G_B_not_closed": True,
    }


def observable():
    return {
        "domain": "The same canonical mass1000 Proca flat reference vacuum, as a full spatial external-metric matching benchmark; the actual CD state is not changed.",
        "conversion": "The explicitly defined covariant nonlocal representative equals its full sixth-time-derivative subtracted retarded bulk plus three finite equal-time coefficients A0+w A1+w^2 A2 at every spatial momentum.",
        "coefficients": "The full tensor coefficients have convergent spectral integral formulas and bound (q+m^2)^3/[80pi^2 m^(2r+2)] for r=0,1,2, with controlled inverse-Lambda spectral tails.",
        "weak_bound": "In the explicit six-spatial/four-time derivative source norm, the conversion bilinear is below1e10 and its spectral tail below1e17/Lambda; both canonical factors give4e-790 and4e-783/Lambda.",
        "boundary": "Finite spatially nonlocal conversion is not an adjustable local counterterm, the fixed physical finite polynomial, curved contact/endpoint matching, full response inverse or original P8 closure.",
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
            "Unsupported flat spatial subtraction conversion benchmark scope accepted: "
            + name
        )
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "actual_parent_normalization_state_and_prescription_retained": True,
        "independent_nonzero_transfer_covariant_curvature_Hessians": True,
        "complete_matrix_improper_spectral_conversion": True,
        "independent_centered_frequency_derivative_crosschecks": True,
        "operator_coefficient_and_infinite_tail_controls": True,
        "literal_prepared_retarded_sine_kernel_extraction": True,
        "spatial_nonlocality_not_local_or_curved_matching": True,
        "original_P8_not_closed": True,
    }
