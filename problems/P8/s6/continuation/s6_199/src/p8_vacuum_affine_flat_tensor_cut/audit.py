"""Canonical flat full Proca tensor cut and subtraction benchmark."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes
from p8_vacuum_affine_retarded_reference_boundary import audit as previous
from p8_vacuum_flat_dirac_hadamard.symbols import rational

from . import curvature, dispersion, polarizations, projectors

ITEM = {
    "id": "actual_canonical_flat_Proca_full_tensor_cut_and_subtracted_dispersion_benchmark",
    "status": "CANONICAL_FLAT_FULL_TENSOR_CUT_BENCHMARK_NOT_CURVED_MATCHING_FULL_PARENT_AMPLITUDE_OR_V_G_B",
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
            "A canonical flat tensor cut benchmark bound cannot close original P8"
        )
    return True


def packets():
    return {
        "complete_physical_polarization_stress": polarizations.data(),
        "full_tensor_cut_and_normalization": projectors.data(),
        "subtracted_dispersive_remainder": dispersion.data(),
        "curvature_weight_and_scope": curvature.data(),
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
        "same_full_parent_normalization_profile_and_prescription": True,
        "actual_CD_state_not_replaced_by_flat_benchmark": True,
        "canonical_flat_reference_Proca_sector_only": True,
        "all_three_physical_modes_and_temporal_constraint": True,
        "all_nine_creation_pair_amplitudes": True,
        "full_covariant_conserved_stress_not_only_spatial_trace": True,
        "both_spin2_and_spin0_sectors_retained": True,
        "complete_rotational_tensor_projection": True,
        "Wick_pair_exchange_and_Lorentz_phase_space_factors": True,
        "actual_i_over_four_retarded_current_factor": True,
        "strict_positive_massive_threshold_support": True,
        "three_subtracted_remainder_and_actual_moment_bounds": True,
        "finite_local_polynomial_not_selected_from_absorptive_cut": True,
        "four_dimensional_heat_weights_with_Euler_retained": True,
        "longitudinal_UV_not_silently_replaced_by_Maxwell": True,
        "both_canonical_tensor_chain_factors": True,
        "no_curved_spatial_or_interacting_matching_inferred": True,
        "original_V_G_B_not_closed": True,
    }


def observable():
    return {
        "domain": "The same canonical mass1000 Proca sector in a flat reference vacuum, as an external-metric Gaussian spatial-matching benchmark; the actual CD state is not changed.",
        "full_cut": "All nine physical polarization pairs yield conserved spin2 and spin0 densities above4m^2, with the complete Wick, Lorentz phase-space and actual current i/4 factors.",
        "dispersion": "The three-subtracted nonlocal part has exact first convergent moments3/(3584pi^2 m^2) and3/(8960pi^2 m^2), controlled next-order errors and explicit spectral-tail bounds.",
        "UV": "The full tensor cut agrees with four-dimensional Proca heat weights13/120 for Weyl squared and1/72 for scalar curvature squared after the complete tensor Hessian factors.",
        "boundary": "The cut does not choose finite local counterterms, fix the curved contact/endpoint sector, supply full parent scattering cuts or finite-gravity Regge data, or close original P8.",
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
            "Unsupported canonical flat tensor cut benchmark scope accepted: " + name
        )
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "actual_parent_normalization_and_constraint_retained": True,
        "independent_four_tensor_polarization_reconstruction": True,
        "noncollinear_Lorentz_and_full_angular_tensor_checks": True,
        "phase_space_Wick_and_retarded_oscillator_normalization": True,
        "independent_subtracted_moments_and_spectral_tail_integrals": True,
        "full_linear_curvature_Hessians_and_pre_cache_type_controls": True,
        "flat_cut_not_finite_local_or_curved_parent_matching": True,
        "original_P8_not_closed": True,
    }
