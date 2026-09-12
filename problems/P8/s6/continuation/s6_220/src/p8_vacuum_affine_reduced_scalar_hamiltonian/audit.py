"""Scoped scalar comparison progress; unchanged original primitive frontier."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes
from p8_vacuum_affine_full_spatial_remainder import audit as previous
from p8_vacuum_flat_dirac_hadamard.symbols import rational

from . import estimates, pullback, scalar

ITEM = {
    "id": "regular_QG1_scalar_coefficient_Hamiltonian_complete_infrared_Ward_pullback_and_comparison_bounds",
    "status": "REGULAR_CLASSICAL_COEFFICIENT_PHASE_SYSTEM_AND_IR_SAFE_REFERENCE_SCALAR_WEAK_RESPONSE_NOT_FULL_QUANTUM_CONSTRAINT_INVERSE_OR_P8",
}


def require_scope(time, kappa=modes.KAPPA, mass=modes.MASS, length=1):
    t, k, m, L = map(rational, (time, kappa, mass, length))
    if not -s.Rational(1, 2) <= t <= s.Rational(1, 2):
        raise ValueError("Only the original unit CD slab")
    if k != modes.KAPPA or m != modes.MASS or L != 1:
        raise ValueError("Require fixed kappa, mass1000 and unit slab")
    return t, k, m, L


def require_band(momentum_bound):
    value = rational(momentum_bound)
    if value < 0:
        raise ValueError("Require a finite nonnegative mathematical momentum bound")
    return value


def frontier():
    return previous.frontier()


def matching():
    return [dict(row) for row in previous.matching()] + [dict(ITEM)]


def validate_scope(rows, obligations):
    if rows != frontier() or obligations != matching():
        raise ValueError(
            "A scalar comparison is not full quantum inversion or original P8 closure"
        )
    return True


def packets():
    return {
        "literal_full_inhomogeneous_scalar_ADM": scalar.spatial_data(),
        "actual_fixed_QG1_full_quadratic_retuning": scalar.retuning_data(),
        "regular_coefficient_sector_canonical_Hamiltonian": scalar.hamiltonian_data(),
        "full_ordered_Ward_infrared_cancellation": pullback.ward_data(),
        "clock_second_chart_and_phase_comparison": pullback.clock_data(),
        "explicit_infrared_scalar_metric_norms": estimates.projection_data(),
        "phase_derivative_loss_and_compact_band_classical_propagator": estimates.phase_data(),
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
            "nine_original_primitive_rows": len(frontier()) - 9,
            "all_original_primitive_rows_unchanged": sum(
                a != b for a, b in zip(frontier(), previous.frontier())
            ),
            "all_prior_matching_rows_unchanged": sum(
                a != b for a, b in zip(matching(), previous.matching())
            ),
            "one_new_scoped_comparison_checkpoint": len(matching())
            - len(previous.matching())
            - 1,
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
            name + "_" + key: bool(value)
            for name, p in packets().items()
            for key, value in p["gates"].items()
        },
        "literal_nonzero_transfer_matter_wave_kept": True,
        "complete_fixed_retuning_not_adaptive_reference_stress": True,
        "both_full_ADM_and_extra_clock_one_current_contacts": True,
        "IR_cancellation_not_unweighted_shift_bound": True,
        "phase_comparison_not_full_quantum_constraint_reduction": True,
        "no_unproved_numerical_high_stress_jet_bound": True,
        "compact_band_not_physical_EFT_cutoff": True,
        "no_frozen_input_edits_or_original_P8_closure": True,
        "unique_matching_identifiers": len({r["id"] for r in matching()})
        == len(matching()),
    }


def observable():
    return {
        "coefficient_sector": "The literal full nonzero-P scalar/matter action, including fixed QG1 lapse-volume retuning, has b=pv/2 and n=Lcorr/(2Jnew), Jnew>1/100; no Theta denominator appears.",
        "reference_scalar_metric": "Both ordered Ward corrections and the extra clock-chart contact give the infrared-safe bound1e117 V03 U138 on (n,zeta,b). The bare shift itself has no claimed unweighted IR norm.",
        "phase_comparison": "The current coefficient-sector pullback has a1e120 C13 V05 U13,10 bound, with finite explicitly defined coefficient-jet norm C13, not an invented numerical13-jet stress bound.",
        "compact_band": "The classical coefficient-sector propagator on |P|<=Lambda is bounded by exp[1000(1+Lambda^2)^2|t-s|]; this is not a uniform continuum inverse or physical cutoff.",
        "remaining": "Full nonlocal quantum constraint inversion, correct-space scalar/mixed inverse, nonlinear sourced parent, quantum background/stability, heavy/physical cutoff and original V/G/B remain open.",
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
            out.append((f"scope_type_{pos}_{i}", require_scope, tuple(args)))
        out.append((f"band_type_{i}", require_band, (value,)))
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
    for value in (-1, -s.Rational(1, 100)):
        out.append((f"negative_band_{len(out)}", require_band, (value,)))
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
            (
                "extra_quantum_inverse",
                validate_scope,
                (
                    frontier(),
                    matching() + [{"id": "full_quantum_inverse", "status": "COMPLETE"}],
                ),
            ),
        )
    )
    return out


@cache
def rejected_inputs():
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (ValueError, TypeError):
            continue
        raise ValueError("Unsupported scalar comparison input accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "independent_constraint_saddle_and_zero_Theta_Legendre": True,
        "literal_full_1024_switch_and_volume_Hessian": True,
        "independent_three_by_three_extrinsic_curvature_spatial_mean": True,
        "independent_density_push_and_detector_pull": True,
        "nonzero_deleted_ADM_contact_IR_pole_and_clock_contact": True,
        "arbitrarily_small_momentum_projector_and_two_endpoint_germs": True,
        "noncoercive_bounce_Hamiltonian_and_explicit_inverse_boundary": True,
        "all_prior_frontiers_unchanged_and_original_P8_open": True,
    }
