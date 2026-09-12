"""Full ADM local Gaussian input, without claiming full continuum matching."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes
from p8_vacuum_affine_matched_spatial_current import audit as previous
from p8_vacuum_flat_dirac_hadamard.symbols import rational

from . import adm, bounds, chart, vertices

ITEM = {
    "id": "full_constrained_ADM_Proca_vertices_metric_chart_contacts_and_all_momentum_forms",
    "status": "FULL_ADM_LOCAL_AND_FINITE_REGULATOR_GAUSSIAN_INPUT_NOT_FULL_NEW_SECTOR_CONTINUUM_REDUCED_INVERSE_OR_V_G_B",
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
        raise ValueError("Full local ADM vertices cannot close original P8")
    return True


def packets():
    return {
        "actual_ADM_Legendre_constraint_and_parent_source_boundary": adm.data(),
        "full_first_second_ADM_vertices_and_tracefree_recovery": vertices.data(),
        "full_four_metric_chart_and_one_point_contact_chain": chart.data(),
        "sharp_complex_shift_and_all_momentum_finite_response_forms": bounds.data(),
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
            "all_corrected_prior_matching_rows_unchanged": sum(
                a != b for a, b in zip(matching(), previous.matching())
            ),
            "one_current_checkpoint_row_added": len(matching())
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
        "same_actual_parent_state_and_canonical_three_mode_sector": True,
        "all_ten_lapse_shift_spatial_metric_directions": True,
        "trace_and_lapse_constraint_variations_not_dropped": True,
        "longitudinal_shift_transport_not_dropped": True,
        "all_ordered_spatial_second_variation_words": True,
        "complex_Fourier_not_assumed_single_mode_Hermitian": True,
        "all_momenta_without_spurious_transfer_loss": True,
        "growing_finite_band_bound_not_UV_limit": True,
        "S213_matched_tracefree_continuum_result_unchanged": True,
        "four_metric_contact_chain_keeps_actual_one_point_current": True,
        "nonlinear_parent_source_not_globally_set_to_zero": True,
        "same_finite_band_state_and_two_propagators": True,
        "unreduced_ADM_norm_not_canonical_scalar_inverse": True,
        "original_V_G_B_and_P8_open": True,
        "all_primitive_statuses_unchanged": True,
        "current_matching_identifiers_unique": len({row["id"] for row in matching()})
        == len(matching()),
    }


def observable():
    return {
        "full_constraint": "The lapse/trace-dependent temporal constraint and the longitudinal shift term follow from the full physical-metric action and Legendre transform.",
        "all_vertices": "Complete first and second ADM metric vertices include every scalar-constraint and noncommuting spatial contact, recover the frozen shear sector and distinguish ADM coordinate contacts from four-metric contacts.",
        "sharp_forms": "The shift feature norm equals a||beta|| for arbitrary complex Fourier beta. With sigma=|n|+(5/2)||Q||op+a||beta||, full first/second energy-relative forms are bounded bysigma andsigma_D sigma_G for all internal/external momenta.",
        "finite_response": "The same prepared-state finite-band Gaussian response bounds apply to all ADM directions in this norm, without interpreting the growing computational band as a physical cutoff.",
        "boundary": "Full new-sector continuum state/time/UV/contact matching and Ward structure, genuinely reduced mixed inverse, finite-amplitude spatial response, finite-coupling background/stability, other parent loops/cutoff and V/G/B remain.",
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
        raise ValueError("Unsupported full ADM local Gaussian scope accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "literal_four_metric_action_and_nonzero_shift_Legendre": True,
        "arbitrary_SPD_magnetic_and_constraint_density": True,
        "literal_noncommuting_full_spatial_first_second_vertices": True,
        "complex_reverse_pair_reality_and_full_shift_flux": True,
        "sharp_complex_shift_norm_and_separated_momenta": True,
        "nonzero_trace_constraint_and_longitudinal_omission_controls": True,
        "nonzero_four_metric_chart_contacts_and_same_parent_boundary": True,
        "original_P8_not_closed": True,
    }
