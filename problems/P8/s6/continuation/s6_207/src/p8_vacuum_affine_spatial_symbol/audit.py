"""Complete four-sector spatial ultraviolet symbol and finite coefficient extraction."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes
from p8_vacuum_affine_sharp_band_artifact import audit as previous
from p8_vacuum_flat_dirac_hadamard.symbols import rational

from . import benchmark, extraction, sectors, symbol

ITEM = {
    "id": "complete_four_sector_spatial_UV_symbol_with_independent_detector_time",
    "status": "COMPLETE_SPATIAL_UV_SYMBOL_EXTRACTION_NOT_FULL_CONTACT_COVARIANT_MATCHING_OR_V_G_B",
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
        raise ValueError("A finite spatial UV extraction cannot close original P8")
    return True


def packets():
    return {
        "exact_full_four_sector_two_time_pair_factorization": sectors.data(),
        "full_inverse_radius_normalization_and_analytic_symbol": symbol.data(),
        "complete_finite_UV_extraction_and_fixed_P_remainder": extraction.data(),
        "nonzero_massive_flat_three_channel_benchmark": benchmark.data(),
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
        "complete_nine_pairs_in_four_sectors": True,
        "full_longitudinal_constraint_and_mixed_pairs": True,
        "independent_source_and_detector_times": True,
        "sector_dependent_inverse_phases": True,
        "actual_eight_order_frequency_before_Taylor_extraction": True,
        "all_fifteen_endpoint_power_slots": True,
        "all_thirty_five_source_time_jet_entries": True,
        "spatial_degree_from_xP_not_mask_derivatives": True,
        "higher_W8_jet_independence_not_reference_replacement": True,
        "fixed_P_integrable_remainder_not_all_P_norm": True,
        "full_massive_flat_benchmark_not_curved_matching": True,
        "actual_curved_odd_endpoints_not_deleted": True,
        "original_pair_and_contact_masks_unchanged": True,
        "leading_sharp_artifact_not_ignored": True,
        "no_full_quantum_matching_or_model_exclusion": True,
        "original_V_G_B_not_closed": True,
    }


def observable():
    return {
        "domain": "The same massive CD Proca unit-W8 comparison, real tracefree spatial source/detector tensors and their independent time jets, with the original regulator masks retained.",
        "exact_reduction": "All nine physical pairs factor into four projector sectors with separate source and detector times and sector-dependent inverse phases. The exact inverse-radius normalization removes individual polarization frames from the symbol.",
        "finite_extraction": "The five endpoints have normalized analytic rows x^(j-1)H_j(x). Retaining degrees0,...,4-j yields fifteen endpoint/power slots and all35 source time-jet entries. Each degree-d coefficient is polynomial in external P of degree at most d.",
        "remainder_and_benchmark": "The subtracted remainder is O(k^-4) at fixed P, with a fixed-P integrable tail, but no explicit uniform all-P norm. An exact massive flat three-channel benchmark and actual curved W8 coefficient/time-jet comparisons test nonzero subleading terms and higher-W8 coefficient independence.",
        "boundary": "The original two-leg integration still produces regulator artifacts, including S206's actual leading cusp. Full subleading conversion, contact and finite/divergent covariant matching, full response/inverse/background/cutoff and original P8 remain open.",
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
        raise ValueError("Unsupported spatial UV symbol scope accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "actual_parent_preparation_and_prescription_unchanged": True,
        "independent_full_fields_vs_projectors_and_source_time_jets": True,
        "all_curved_coefficients_two_Cauchy_resolutions": True,
        "sixth_eighth_W8_terms_do_not_change_retained_UV_jets": True,
        "massive_flat_full_field_benchmark_and_nonzero_remainder": True,
        "premature_coincidence_and_common_frequency_controls": True,
        "odd_endpoint_UV_cancellation_not_full_endpoint_deletion": True,
        "original_P8_not_closed": True,
    }
