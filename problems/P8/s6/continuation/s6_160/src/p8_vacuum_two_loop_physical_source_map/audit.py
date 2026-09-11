"""Preserve broad matching obligations while recording physical-source transport."""

from functools import cache

import sympy as s
from p8_vacuum_full_heavy_source import audit as previous

from . import calibration, topology, transport, vertices, ward

MODULES = (topology, vertices, ward, transport, calibration)
ITEM = {
    "id": "complete_two_loop_derivative_map_with_matched_physical_source",
    "status": "TRANSPORTED_FOR_FULL_ACTION_AND_PHYSICAL_SOURCE_NOT_ORDINARY_PSI_MS_OR_TRUNCATION",
}


def frontier():
    return previous.frontier()


def matching():
    return [dict(r) for r in previous.matching()] + [dict(ITEM)]


def validate_scope(rows, obligations):
    if rows != frontier() or obligations != matching():
        raise ValueError(
            "Physical-source transport scope differs from the fixed ledger"
        )
    return True


@cache
def residuals():
    out = {
        mod.__name__.rsplit(".", 1)[-1] + "_" + k: s.simplify(v)
        for mod in MODULES
        for k, v in mod.data()["checks"].items()
    }
    out.update(
        {
            "all_nine_primitive_rows_unchanged": sum(
                a != b for a, b in zip(frontier(), previous.frontier())
            ),
            "nine_primitive_rows_retained": len(frontier()) - 9,
            "all_parent_matching_rows_unchanged": sum(
                a != b for a, b in zip(matching(), previous.matching())
            ),
            "only_physical_source_transport_item_added": len(matching())
            - len(previous.matching())
            - 1,
        }
    )
    return out


def scalar_entry_count():
    return len(residuals())


@cache
def gates():
    return {
        **calibration.data()["bounds"],
        "literal_fixed_cubic_derivative_map_used": True,
        "full_scalar_action_including_octic_and_higher_kept": True,
        "generated_Yukawa_and_parent_counterterms_retained": True,
        "physical_source_four_point_requires_map_order_three": True,
        "source_degree_and_counterterm_loop_weights_counted": True,
        "local_H_source_and_vacuum_counterterms_kept_separately": True,
        "arbitrary_polynomial_Gaussian_Ward_identity_through_cubic_map_order": True,
        "finite_dimensional_Jacobian_kept_in_diagnostics": True,
        "continuum_DR_ghost_loops_scaleless_before_finite_parts": True,
        "same_regulated_action_source_and_counterterm_prescription": True,
        "heavy_elimination_and_physical_H_reference_unchanged": True,
        "physical_source_correlators_equal_canonical_parent_at_named_order": True,
        "no_extra_LSZ_factor_on_already_physical_source": True,
        "ordinary_Psi_composite_MS_mixing_not_claimed_computed": True,
        "no_global_inverse_or_bare_derivative_truncation_equivalence": True,
        "physical_truncation_V_G_B_and_original_P8_open": True,
    }


def bad_cases():
    invalid = (
        True,
        False,
        1.0,
        s.Integer(2),
        "2",
        None,
        s.oo,
        -s.oo,
        s.zoo,
        s.I,
        s.nan,
        s.Symbol("x"),
    )
    out = []
    for j in range(2):
        for i, v in enumerate(invalid):
            args = [4, 2]
            args[j] = v
            out.append(
                (f"map_degree_type_{j}_{i}", topology.max_map_degree, tuple(args))
            )
    for e, l in ((-2, 2), (1, 2), (3, 2), (6, 2), (4, -1), (4, 3)):
        out.append((f"map_degree_domain_{e}_{l}", topology.max_map_degree, (e, l)))
    for i in range(len(topology.support_rows())):
        rows = topology.support_rows()
        rows[i]["map_degree"] += 1
        out.append((f"operator_support_{i}", topology.validate_support, (rows,)))
    out.extend(
        (
            (
                "missing_operator",
                topology.validate_support,
                (topology.support_rows()[:-1],),
            ),
            (
                "duplicate_operator",
                topology.validate_support,
                (topology.support_rows() + topology.support_rows()[:1],),
            ),
        )
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
        raise ValueError(
            "Unsupported physical-source transport input accepted: " + name
        )
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "octic_omission_detected_at_two_loop_four_point": True,
        "generated_Yukawa_omission_detected_at_one_loop_four_point": True,
        "source_and_Jacobian_omissions_detected_separately": True,
        "map_order_two_not_mistaken_for_sufficient_at_two_loops": True,
        "physical_source_not_an_independent_ordinary_Psi_MS_operator": True,
        "no_uncomputed_ordinary_residue_or_global_norm_bound": True,
        "frozen_parent_Phi_pole_b2_vacuum_H_references_unchanged": True,
        "original_P8_not_closed": True,
    }
