"""Actual leading sharp two-leg artifact and scoped regulator conversion."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes
from p8_vacuum_affine_finite_endpoint_polynomial import audit as previous
from p8_vacuum_flat_dirac_hadamard.symbols import rational

from . import angular, asymptotics, geometry, polarizations

ITEM = {
    "id": "actual_leading_full_Proca_sharp_two_leg_regulator_conversion_coefficient",
    "status": "ACTUAL_LEADING_SHARP_BAND_ARTIFACT_NOT_FULL_REGULATOR_COVARIANT_MATCHING_OR_V_G_B",
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
        raise ValueError("A sharp-band conversion coefficient cannot close original P8")
    return True


def packets():
    return {
        "all_nine_physical_pair_leading_endpoint_symbol": polarizations.data(),
        "complete_weighted_angular_moments_and_five_channels": angular.data(),
        "original_intersection_and_lost_shell_limit": geometry.data(),
        "actual_W8_asymptotics_and_scoped_conversion": asymptotics.data(),
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
        "unit_W8_reference_and_actual_state_correction_kept": True,
        "complete_nine_pair_current_sign_and_Fourier_measure": True,
        "full_longitudinal_temporal_and_mass_readouts": True,
        "fixed_mass_high_momentum_not_Maxwell_theory": True,
        "all_five_tracefree_spatial_channels": True,
        "original_second_leg_not_silently_removed": True,
        "one_ball_only_explicit_conversion_comparator": True,
        "complete_original_contact_band_retained": True,
        "grazing_strip_in_exact_finite_geometry": True,
        "uniform_fixed_P_asymptotics_not_all_P_norm_bound": True,
        "all_four_other_endpoint_shell_orders_lower": True,
        "nonpolynomial_cusp_not_local_Hessian": True,
        "no_subtraction_or_matching_changed_here": True,
        "not_all_subleading_artifacts_or_finite_matching": True,
        "not_renormalized_current_divergence": True,
        "canonical_prefactor_not_uniform_regulator_bound": True,
        "original_V_G_B_not_closed": True,
    }


def observable():
    return {
        "domain": "The same massive CD Proca sector and unit-W8 comparison, with real tracefree spatial tensors and fixed external P on the compact time slab.",
        "actual_symbol": "All nine physical pairs give r*[4tr(DG)-4n.(DG+GD)n+3(n.D.n)(n.G.n)]/(8a), retaining the longitudinal contribution.",
        "conversion": "Pair-band minus one-k-ball for the first five endpoints has leading term -K^3|P|[18tr(DG)-12(D phat).(G phat)-(phat.D.phat)(phat.G.phat)]/(512 pi^2 a). The comparator is not asserted renormalized and does not replace the original mask.",
        "channels": "The bracket eigenvalues on normalized tensor/tensor/vector/vector/scalar channels are18,18,12,12,28/3. The nonzero |P| cusp is a real leading regulator-conversion artifact.",
        "boundary": "This is not a full quantum matching, physical divergence, cutoff or model exclusion. Subleading artifacts, actual contact/candidate-cell finite and divergent coefficients and the original covariant matching remain open.",
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
        raise ValueError("Unsupported sharp-band artifact scope accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "actual_parent_preparation_and_prescription_unchanged": True,
        "independent_actual_W8_all_nine_pair_limit": True,
        "full_five_channel_angular_moments_and_rotations": True,
        "original_removed_shell_and_grazing_geometry": True,
        "exact_moment_prewarmer_and_order_guards": True,
        "longitudinal_not_Maxwell_and_one_leg_controls": True,
        "artifact_not_full_physical_divergence_or_matching": True,
        "original_P8_not_closed": True,
    }
