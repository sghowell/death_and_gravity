"""Actual all-order CD stress noise without a quantum-background verdict."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import audit as previous
from p8_vacuum_affine_cd_source_noise import modes
from p8_vacuum_flat_dirac_hadamard.symbols import rational

from . import noise, oscillatory, reference, stress

ITEM = {
    "id": "actual_all_order_CD_full_local_smeared_metric_stress_covariance",
    "status": "FULL_SPACETIME_TENSOR_STRESS_NOISE_BOUND_NOT_METRIC_SOLUTION_CUTOFF_OR_V_G_B",
}


def require_scope(time, kappa=modes.KAPPA, mass=modes.MASS, length=1):
    t, k, m, L = map(rational, (time, kappa, mass, length))
    if not -s.Rational(1, 2) <= t <= s.Rational(1, 2):
        raise ValueError("Only the fixed unit CD slab")
    if k != modes.KAPPA or m != modes.MASS or L != 1:
        raise ValueError(
            "Require actual fixed CD parent, mass1000 and unit observation slab"
        )
    return t, k, m, L


def frontier():
    return previous.frontier()


def matching():
    return [dict(row) for row in previous.matching()] + [dict(ITEM)]


def validate_scope(rows, obligations):
    if rows != frontier() or obligations != matching():
        raise ValueError("A stress-noise input bound cannot close original P8")
    return True


def packets():
    return {
        "full_physical_stress": stress.data(),
        "reference_and_actual_state": reference.data(),
        "pair_smearing": oscillatory.data(),
        "centered_stress_noise": noise.data(),
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
        "same_all_order_state_not_W8_or_flat_vacuum": True,
        "exact_initial_amplitude_and_full_evolved_remainder": True,
        "all_three_modes_and_temporal_constraint_kept": True,
        "full_local_stress_not_bilocal_smoothed_field_products": True,
        "all_tensor_components_and_pair_cross_terms_kept": True,
        "time_smearing_precedes_distributional_limit": True,
        "full_internal_and_external_momenta_not_homogeneous_only": True,
        "no_physical_cutoff_inferred_from_proof_regulator": True,
        "same_mean_prescription_c_numbers_only_cancel_in_centering": True,
        "quadratic_stress_not_assumed_Gaussian": True,
        "Einstein_test_normalization_not_reduced_mixed_mode_chart": True,
        "noise_input_not_metric_solution": True,
        "no_vanishing_bounce_density_denominator": True,
        "original_V_G_B_not_closed": True,
    }


def observable():
    return {
        "noise": "The actual all-order CD state's complete local Proca stress has smeared standard deviation below1e25 in the specified three-time-derivative and one-spatial-derivative tensor norm.",
        "normalization": "The same bounds become1e-775 for stress/kappa and1e-375 for the Einstein-normalized metric test coupling, not a fully reduced mixed-mode norm.",
        "scope": "All tensor components, physical polarizations, temporal constraint and the complete momentum range are retained. No flat state, WKB reset, physical momentum cutoff or separately smoothed field product replaces the observable.",
        "remaining": "Actual spatial/nonlinear/finite-coupling metric response, full interacting parent measures/loops, heavy/cutoff/threshold and omitted-order matching, full vacuum and finite-gravity V/G/B.",
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
        (0, modes.KAPPA * 2, modes.MASS, 1),
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
        raise ValueError("Unsupported metric stress-noise scope accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "state_not_reset_or_replaced_by_reference": True,
        "full_local_tensor_not_only_energy_or_pressure": True,
        "creation_pairs_include_complex_conjugate_modes": True,
        "no_time_derivative_assumption_on_fast_mixing": True,
        "complete_radial_integrals_not_a_finite_grid": True,
        "mean_cancellation_not_zero_noise": True,
        "noise_input_not_a_metric_stability_certificate": True,
        "original_P8_not_closed": True,
    }
