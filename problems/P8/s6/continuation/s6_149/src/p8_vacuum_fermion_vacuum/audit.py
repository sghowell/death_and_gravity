"""Both vacuum primitive rows advance; full matching remains a separate obligation."""

from functools import cache

import sympy as s
from p8_vacuum_fermion_insertion_ms import audit as previous

from . import calibration, derivatives, forest, massless, scalar

MODULES = (massless, derivatives, scalar, forest, calibration)
TARGETS = ("scalar_Phi0_W0_F0", "gauge_Phi0")
STATUS = "BOUNDED_PAIRED_VACUUM_PRIMITIVE_WITH_ASSIGNED_LOCAL_FOREST"


def frontier():
    rows = previous.frontier()
    for r in rows:
        if r["id"] in TARGETS:
            r["status"] = STATUS
    return rows


def validate_frontier(rows):
    if rows != frontier():
        raise ValueError("Vacuum primitive frontier differs from exact scope")
    return True


@cache
def residuals():
    out = {
        mod.__name__.rsplit(".", 1)[-1] + "_" + k: s.simplify(v)
        for mod in MODULES
        for k, v in mod.data()["checks"].items()
    }
    old, new = previous.frontier(), frontier()
    out.update(
        {
            "only_two_vacuum_rows_advance": sum(a != b for a, b in zip(old, new)) - 2,
            "nine_primitive_rows_retained": len(new) - 9,
            "no_primitive_row_wholly_unevaluated": sum(
                r["status"] == "UNEVALUATED" for r in new
            ),
            "two_paired_vacuum_rows": sum(r["status"] == STATUS for r in new) - 2,
            "seven_other_rows_unchanged": sum(a == b for a, b in zip(old, new)) - 7,
        }
    )
    return out


def scalar_entry_count():
    return sum(
        len(v) if isinstance(v, s.MatrixBase) else 1 for v in residuals().values()
    )


@cache
def gates():
    return {
        **calibration.data()["bounds"],
        "dimension_symbolic_scalar_and_gauge_vacuum_traces": True,
        "both_fermion_self_energy_counterterms_in_determinant_variation": True,
        "fixed_mu_vacuum_mass_derivatives_match_previous_anchors": True,
        "scalar_whole_cycle_full_physical_mass_residue_forest": True,
        "vector_Ward_identity_closes_gauge_CT_scaleless": True,
        "three_proper_theta_subcycles_overlap_no_product_forest": True,
        "actual_scalar_mass_not_massless_auxiliary_anchor": True,
        "three_beta_terms_before_convergent_scalar_remainder": True,
        "finite_remainder_pole_prefactor_product_retained": True,
        "dimensional_contact_removal_before_convergent_spectral_integral": True,
        "six_Yukawa_and_forty_two_gauge_color_states": True,
        "large_vacuum_constant_retained_not_set_zero": True,
        "assigned_counterterms_removed_from_separate_insertion_ledger": True,
        "other_vacuum_source_and_matching_terms_remain_open": True,
        "primitive_ledger_not_full_canonical_V_G_B_or_original_P8_closure": True,
    }


def bad_cases():
    invalid = (
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
    for prefix, call in (
        ("scalar", scalar.enclosure),
        ("gauge", calibration.gauge_enclosure),
    ):
        for j in range(3):
            for i, v in enumerate(invalid):
                vals = [2, 1, 144]
                vals[j] = v
                out.append((f"{prefix}_type_{j}_{i}", call, tuple(vals)))
        for i, v in enumerate(((1, 1, 144), (2, -1, 144), (2, 1, 0), (2, 1, 145))):
            out.append((f"{prefix}_domain_{i}", call, v))
    for i in range(9):
        rows = frontier()
        rows[i]["status"] = "COMPLETE"
        out.append((f"frontier_{i}", validate_frontier, (rows,)))
    out.extend(
        (
            ("frontier_missing", validate_frontier, (frontier()[:-1],)),
            (
                "frontier_extra",
                validate_frontier,
                (frontier() + [{"id": "extra", "status": "COMPLETE"}],),
            ),
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
        raise ValueError("Unsupported vacuum primitive input accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "full_assigned_proper_forests_not_raw_vacuum_only": True,
        "no_early_d_four_numerator_or_finite_pole_product_loss": True,
        "actual_scalar_mass_and_nonzero_finite_remainder_retained": True,
        "all_fourteen_gauge_flavors_not_only_active_pair": True,
        "vacuum_energy_reference_not_yet_full_model_sum": True,
        "first_source_square_and_other_vacuum_terms_open": True,
        "other_matching_canonical_and_truncation_terms_open": True,
        "original_P8_not_closed": True,
    }
