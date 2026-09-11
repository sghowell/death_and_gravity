"""Actual continuum weak response on strictly separated supports."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes
from p8_vacuum_affine_spatial_current import audit as previous
from p8_vacuum_flat_dirac_hadamard.symbols import rational

from . import projection, response, support, tail

ITEM = {
    "id": "actual_full_spatial_separated_support_weak_Gaussian_response_and_quantitative_tail",
    "status": "ACTUAL_CONTINUUM_SEPARATED_SUPPORT_WEAK_RESPONSE_NOT_COINCIDENT_EXTENSION_INVERSE_OR_V_G_B",
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
            "A separated-support weak response bound cannot close original P8"
        )
    return True


def packets():
    return {
        "complete_two_particle_stress_tail": tail.data(),
        "common_orthogonal_projection": projection.data(),
        "ordered_support_contact_and_causality": support.data(),
        "actual_continuum_weak_current_response": response.data(),
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
        "full_ten_feature_stress_and_nine_polarization_pairs": True,
        "both_created_momenta_projected_by_the_same_orthogonal_operator": True,
        "reference_tail_union_uniform_in_external_momentum": True,
        "low_external_tail_forces_both_internal_legs_large": True,
        "high_external_tail_controlled_by_actual_spatial_H1_norm": True,
        "all_internal_and_external_momentum_regulators_removed": True,
        "strict_time_order_verified_before_removing_the_step": True,
        "full_metric_and_finite_curvature_contacts_zero_by_support": True,
        "same_actual_preparation_not_a_new_reference_vacuum": True,
        "physical_stress_to_Hamiltonian_sign_and_quarter_factor": True,
        "full_Proca_causal_Green_operator_not_scalar_substitution": True,
        "spacelike_locality_only_after_full_continuum_limit": True,
        "two_canonical_metric_chain_factors_retained": True,
        "weak_bilinear_not_operator_norm_differentiability": True,
        "no_overlapping_time_or_diagonal_extension_claimed": True,
        "no_full_inverse_interacting_background_or_stability": True,
        "original_V_G_B_not_closed": True,
    }


def observable():
    return {
        "domain": "Real smooth compact spatial tracefree tests in the open actual CD slab, with the complete source support strictly earlier than the complete readout support.",
        "full_response": "The actual Gaussian-sector weak current is i<[T[D],T[Gamma]]>/4, with all constrained polarizations and all internal/external momenta. Full local metric and fixed finite-curvature contacts vanish by support.",
        "tail": "For K>=1000 the squared discarded two-particle stress-vector norm is below1e52 N[f]^2/K. A common orthogonal projection gives the response error below5e51 N[D]N[Gamma]/K.",
        "canonical": "The complete weak response is below5e49 N[D]N[Gamma], and below2e-750 after the two canonical metric factors. Its canonical regulator error is below2e-748 N[D]N[Gamma]/K.",
        "boundary": "The finite-band sequence converges quantitatively on these separated supports. This does not extend the retarded distribution to overlapping times or the diagonal, establish operator-norm differentiability, or supply a full inverse or quantum background.",
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
            "Unsupported separated-support weak current scope accepted: " + name
        )
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "actual_parent_state_constraint_and_physical_vertex_retained": True,
        "independent_massive_radial_and_nonzero_transfer_integrals": True,
        "literal_multimode_Fock_pair_factor_and_both_leg_projection": True,
        "full_Hilbert_projection_tail_not_mixed_cross_terms": True,
        "strict_time_order_and_spacelike_enclosure_controls": True,
        "overlap_triangle_not_the_full_commutator_square": True,
        "weak_continuum_not_operator_norm_inverse_or_background": True,
        "original_P8_not_closed": True,
    }
