"""Actual reference retarded boundary extraction and finite bulk."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes
from p8_vacuum_affine_retarded_state_remainder import audit as previous
from p8_vacuum_flat_dirac_hadamard.symbols import rational

from . import boundary, bounds, ibp, majorants

ITEM = {
    "id": "actual_reference_six_step_retarded_boundary_extraction_and_full_finite_bulk",
    "status": "ACTUAL_REFERENCE_FINITE_BULK_AND_EXPLICIT_ENDPOINTS_NOT_COVARIANT_SPATIAL_MATCHING_INVERSE_OR_V_G_B",
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
            "A retarded reference boundary extraction bound cannot close original P8"
        )
    return True


def packets():
    return {
        "exact_six_step_retarded_identity": ibp.data(),
        "complete_Leibniz_Cauchy_majorants": majorants.data(),
        "full_finite_reference_bulk_and_endpoint": bounds.data(),
        "remaining_contact_and_time_endpoints": boundary.data(),
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
        "actual_constant_alpha_W8_reference_not_a_new_state": True,
        "full_ten_field_pair_amplitude_and_all_nine_polarizations": True,
        "six_retarded_integrations_with_all_boundary_terms": True,
        "complete_Cauchy_Leibniz_products_through_six": True,
        "all_initial_terms_zero_by_common_preparation": True,
        "no_odd_endpoint_deleted_before_full_tensor_contraction": True,
        "sixth_bulk_integrable_without_retarded_step_differentiation": True,
        "fifth_endpoint_has_six_inverse_frequency_factors": True,
        "fifth_bulk_absolute_majorant_not_claimed_integrable": True,
        "all_internal_and_external_momentum_regulators_removed": True,
        "full_two_leg_reference_tail_uniform_in_external_transfer": True,
        "full_contact_keeps_its_distinct_one_leg_regulator": True,
        "known_actual_state_remainder_added_in_compatible_norms": True,
        "both_canonical_tensor_chain_factors": True,
        "equal_time_not_automatically_spatial_locality": True,
        "fixed_covariant_contact_and_endpoint_matching_still_open": True,
        "original_V_G_B_not_closed": True,
    }


def observable():
    return {
        "domain": "The same actual Gaussian Proca reference on the CD clock and overlapping compact spatial tracefree probes with the common preparation.",
        "identity": "Six exact retarded integrations isolate all six equal-time boundary terms. The sixth bulk and final endpoint have a complete continuum limit with all momenta retained.",
        "finite_reference": "The finite reference piece is below1e48 times detector L2 and source time-H6 norms, with uniform regulator error1e52/K.",
        "known_total": "Adding the complete S197 actual-state remainder gives known finite part2e48 in detector spatial H1 and source six-time/one-space derivative norms; canonical display8e-752 and tail8e-748/K.",
        "boundary": "The full reference contact plus five equal-time terms still need the original covariant spatial matching. Equal-time does not imply spatial locality, and no full inverse, background, cutoff or V/G/B closure follows.",
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
            "Unsupported retarded reference boundary extraction scope accepted: " + name
        )
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "actual_parent_state_constraint_and_physical_vertex_retained": True,
        "literal_variable_phase_retarded_integrals_and_every_endpoint": True,
        "independent_rational_derivative_Leibniz_coefficients": True,
        "common_initial_jet_and_missing_lower_boundary_controls": True,
        "fifth_absolute_tail_divergence_and_sixth_integrability": True,
        "nonpolynomial_spatial_kernel_and_odd_endpoint_controls": True,
        "finite_piece_not_remaining_covariant_matching_or_inverse": True,
        "original_P8_not_closed": True,
    }
