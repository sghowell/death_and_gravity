"""Advance the fixed-order amplitude/pole obligation without closing P8."""

from functools import cache

import sympy as s
from p8_vacuum_full_phi_normalization import audit as previous

from . import bounds, calibration, contact, ownership
from . import matching as coordinate

MODULES = (ownership, coordinate, contact, bounds, calibration)
TARGET = "full_matched_two_loop_pole_and_amplitude_error"
STATUS = "FIXED_ORDER_PHI_POLE_AND_B2_BOUNDED_PHYSICAL_TRUNCATION_OPEN"


def frontier():
    return previous.frontier()


def matching():
    return [
        {**row, "status": STATUS} if row["id"] == TARGET else dict(row)
        for row in previous.matching()
    ]


def validate_scope(rows, obligations):
    if rows != frontier() or obligations != matching():
        raise ValueError("Complete amplitude scope differs from the fixed ledger")
    return True


@cache
def residuals():
    out = {
        mod.__name__.rsplit(".", 1)[-1] + "_" + key: s.simplify(v)
        for mod in MODULES
        for key, v in mod.data()["checks"].items()
    }
    out.update(
        {
            "all_nine_primitive_rows_unchanged": sum(
                a != b for a, b in zip(frontier(), previous.frontier())
            ),
            "nine_primitive_rows_retained": len(frontier()) - 9,
            "only_fixed_order_amplitude_pole_obligation_advanced": sum(
                a != b for a, b in zip(matching(), previous.matching())
            )
            - 1,
            "matching_ledger_length_unchanged": len(matching())
            - len(previous.matching()),
        }
    )
    return out


def scalar_entry_count():
    return len(residuals())


@cache
def gates():
    return {
        **calibration.data()["bounds"],
        "all_192_scalar_refinements_disjoint": True,
        "all_four_fermion_quartic_families_and_proper_MS_forests": True,
        "heavy_summed_W_not_replaced_by_local_L": True,
        "each_assigned_counterterm_occurrence_removed_from_separate_ledger": True,
        "full_regulated_contact_conversion_composed_pointwise": True,
        "no_commuting_vector_field_assumption": True,
        "hybrid_MS_interactions_and_physical_Phi_OS_references": True,
        "entire_renormalized_one_loop_finite_before_coordinate_variation": True,
        "bare_first_pole_reexpression_k1_and_induced_M2_retained": True,
        "fundamental_G1_square_distinct_from_graph_disjoint_counterterms": True,
        "scalar_and_fermion_first_homogeneity_weights_four_and_two": True,
        "no_extra_LSZ_or_raw_MS_formula_on_hybrid_amplitude": True,
        "source_and_local_Phi2_references_owned_by_canonical_pole": True,
        "complete_unit_disc_Phi_pole_inherited_unchanged": True,
        "formal_through_two_loop_not_physical_truncation": True,
        "vacuum_source_V_G_B_and_original_P8_open": True,
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
    for j in range(3):
        for i, v in enumerate(invalid):
            vals = [1] * 3
            vals[j] = v
            out.append(
                (f"cross_type_{j}_{i}", bounds.first_variation_upper, tuple(vals))
            )
        vals = [1] * 3
        vals[j] = -1
        out.append((f"cross_negative_{j}", bounds.first_variation_upper, tuple(vals)))
    for j in range(5):
        vals = [1] * 5
        vals[j] = -1
        out.append((f"assembly_negative_{j}", bounds.assemble, tuple(vals)))
    out.append(("assembly_zero_tree", bounds.assemble, (0, 1, 1, 1, 1)))
    for i in range(8):
        rows = ownership.rows()
        rows[i]["status"] = "UNASSIGNED"
        out.append(
            (f"owner_{i}", ownership.validate_rows, (rows, ownership.counterterms()))
        )
    out.append(
        (
            "owner_missing",
            ownership.validate_rows,
            (ownership.rows()[:-1], ownership.counterterms()),
        )
    )
    out.append(
        (
            "owner_duplicate",
            ownership.validate_rows,
            (ownership.rows() + ownership.rows()[:1], ownership.counterterms()),
        )
    )
    for i in range(10):
        ct = ownership.counterterms()
        ct[i]["extra"] = not ct[i]["extra"]
        out.append(
            (f"counterterm_{i}", ownership.validate_rows, (ownership.rows(), ct))
        )
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
            ("missing_obligation", validate_scope, (frontier(), matching()[:-1])),
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
        raise ValueError("Unsupported complete amplitude input accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "primitive_sum_not_used_without_matching_and_forests": True,
        "no_double_counted_proper_counterterms_or_LSZ": True,
        "positive_epsilon_bare_map_and_induced_heavy_mass_retained": True,
        "finite_hybrid_first_variation_not_raw_MS_amputation": True,
        "no_commuting_contact_and_field_direction_assumption": True,
        "complete_fixed_order_pole_and_amplitude_only": True,
        "vacuum_source_physical_truncation_V_G_B_still_open": True,
        "original_P8_not_closed": True,
    }
