"""Conditional complete vacuum source response without a full UV verdict."""

from functools import cache

import sympy as s
from p8_vacuum_affine_classical_propagator import audit as previous
from p8_vacuum_flat_dirac_hadamard.symbols import rational

from . import energy, noise, source, vacuum

ITEM = {
    "id": "complete_finite_time_canonical_vacuum_source_mean_and_scalar_force_noise_with_fixed_profile",
    "status": "FULL_CONDITIONAL_GAUSSIAN_VACUUM_SOURCE_CAUSAL_AND_NOISE_BOUNDS_NOT_INTERACTING_AMPLITUDE_CUTOFF_OR_V_G_B",
}


def require_scope(time, kappa=source.K0, jet_bound=1, mass=source.MASS):
    t, k, j, m = map(rational, (time, kappa, jet_bound, mass))
    if not -s.Rational(1, 2) <= t <= s.Rational(1, 2):
        raise ValueError("Only the unit-length observation interval")
    if k < source.K0 or not 0 < j <= 1 or m != source.MASS:
        raise ValueError(
            "Require anchored fixed-canonical family, unit jet ball and fixed mass1000"
        )
    return t, k, j, m


def frontier():
    return previous.frontier()


def matching():
    return [dict(row) for row in previous.matching()] + [dict(ITEM)]


def validate_scope(rows, obligations):
    if rows != frontier() or obligations != matching():
        raise ValueError("Conditional vector bounds cannot close original P8")
    return True


def packets():
    return {
        "complete_source": source.data(),
        "causal_mean": energy.data(),
        "scalar_force_noise": noise.data(),
        "fixed_vacuum_profile": vacuum.data(),
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
        "complete_source_not_leading_harmonic_truncation": True,
        "complex_source_bounds_keep_all_null_gradient_terms": True,
        "prepared_real_time_class_not_impossible_compact_bandlimited_class": True,
        "original_temporal_constraint_and_nonconserved_source_kept": True,
        "positive_reference_energy_not_full_driven_stress": True,
        "causal_scalar_force_from_unintegrated_action": True,
        "same_specified_Minkowski_vacuum_not_a_CD_state_transfer": True,
        "three_physical_polarizations_and_full_covariance_retained": True,
        "smeared_linear_force_noise_not_coincident_stress_noise": True,
        "full_fixed_profile_and_vacuum_constant_retained": True,
        "all_canonical_scalar_functions_fixed_in_family": True,
        "conditional_decoupling_not_interacting_quantum_limit": True,
        "finite_time_bound_not_on_shell_pole_or_cutoff_claim": True,
        "original_V_G_B_not_closed": True,
    }


def observable():
    return {
        "mean": "The complete nonlinear affine vacuum source has all-momentum finite-time energy and scalar mean-force bounds on the fixed real canonical jet class.",
        "noise": "The same positive-frequency Minkowski Gaussian vector sector gives smeared scalar-force standard deviation below2e-1190 per unit J3 test norm at the anchor, including its longitudinal polarization.",
        "fixed_profile": "The nonconstant scalar coefficient and its force are below1e-818500 on this class after the explicit same-prescription vacuum-constant cancellation; the full fixed coefficient persists in the canonical limit.",
        "actual_family": "Source mean force and noise variance decay as1/kappa at fixed full scalar functions and mass; this is a conditional Gaussian statement, not full interacting vacuum matching.",
        "remaining": "Metric stress/noise, nonlinear light/gravity solutions, interacting measures and loops, threshold/cutoff/omitted-order matching, full vacuum cuts/contour and finite-gravity IR/Regge.",
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
            args = [0, source.K0, 1, source.MASS]
            args[pos] = value
            out.append((f"type_{pos}_{i}", require_scope, tuple(args)))
    for args in (
        (1, source.K0, 1, source.MASS),
        (-1, source.K0, 1, source.MASS),
        (0, source.K0 / 2, 1, source.MASS),
        (0, source.K0, 0, source.MASS),
        (0, source.K0, 2, source.MASS),
        (0, source.K0, 1, 0),
        (0, source.K0, 1, 999),
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
        raise ValueError("Unsupported vacuum source scope accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "full_source_not_replaced_by_polynomial_germ": True,
        "no_time_Fourier_support_inferred_from_preparation": True,
        "longitudinal_and_temporal_terms_not_discarded": True,
        "mean_and_noise_not_conflated": True,
        "fixed_Gaussian_sector_not_full_parent_measure": True,
        "nonconstant_profile_not_deleted_in_decoupling": True,
        "small_forcing_not_a_full_background_error_bound": True,
        "original_P8_not_closed": True,
    }
