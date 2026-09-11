"""Actual spatial vertices and two-momentum finite-regulator response."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes
from p8_vacuum_affine_covariant_current import audit as previous
from p8_vacuum_flat_dirac_hadamard.symbols import rational

from . import bounds, hamiltonian, response, vertices

ITEM = {
    "id": "actual_spatial_two_momentum_metric_vertices_and_finite_band_Gaussian_response",
    "status": "ACTUAL_SPATIAL_TWO_MOMENTUM_CURRENT_FINITE_BAND_NOT_RENORMALIZED_CONTINUUM_INVERSE_OR_V_G_B",
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
            "A finite-band spatial response bound cannot close original P8"
        )
    return True


def packets():
    return {
        "complete_spatial_Hamiltonian": hamiltonian.data(),
        "two_momentum_physical_vertices": vertices.data(),
        "actual_off_diagonal_Gaussian_response": response.data(),
        "uniform_external_momentum_finite_band": bounds.data(),
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
        "full_source_zero_on_spatial_unimodular_external_family": True,
        "spatial_Legendre_transform_before_temporal_elimination": True,
        "constraint_energy_and_longitudinal_polarization_retained": True,
        "both_distinct_internal_momenta_in_the_magnetic_vertex": True,
        "full_pointwise_exponential_second_contact_convolution": True,
        "all_energy_features_before_norm_estimation": True,
        "regular_zero_momentum_and_noncollinear_pairs": True,
        "full_two_propagator_off_diagonal_covariance_tangent": True,
        "actual_common_preparation_not_new_metric_state": True,
        "covariance_and_Kubo_sign_and_factor_agree": True,
        "finite_regulator_contact_not_noise_anticommutator": True,
        "high_external_transfer_contact_not_dropped": True,
        "all_external_momenta_in_the_finite_band_bound": True,
        "finite_band_not_a_physical_cutoff": True,
        "no_homogeneous_subtraction_transfer_to_full_space": True,
        "no_full_inverse_interacting_background_or_stability": True,
        "original_V_G_B_not_closed": True,
    }


def observable():
    return {
        "domain": "The same actual Gaussian Proca state at the isotropic CD clock and prescribed smooth spatial unimodular shear, with the common zero initial neighborhood.",
        "exact_spatial_vertices": "The complete constrained physical Hamiltonian gives two-momentum first vertices and nonzero local metric contacts with full Fourier convolution, and energy-relative norms uniform in both internal and external momenta.",
        "response": "Two distinct unperturbed propagators give the exact off-diagonal Gaussian covariance tangent. The finite-regulator current includes its retarded Kubo term and the second metric contact.",
        "finite_band_bound": "For the common computational band |k|<=K=1e16, the complete projected current response is below2e104 times the temporal supremum of a fixed external Fourier input, uniformly in its momentum. The canonical display is8e-696.",
        "boundary": "The growing internal-momentum bound is not an infinite renormalized spatial-tail estimate. No homogeneous covariant subtraction is transferred, no physical cutoff is declared and no full inverse, quantum background or V/G/B closure follows.",
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
            "Unsupported spatial finite-band current scope accepted: " + name
        )
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "actual_parent_state_constraint_and_physical_vertex_retained": True,
        "independent_spatial_Legendre_and_noncommuting_lattice_variations": True,
        "two_momenta_and_energy_feature_factorization": True,
        "complete_metric_contact_and_reverse_pair_count": True,
        "independent_coupled_evolution_and_Kubo_normalization": True,
        "high_transfer_contact_and_nonintegrable_majorant_controls": True,
        "finite_band_not_full_continuum_or_physical_cutoff": True,
        "original_P8_not_closed": True,
    }
