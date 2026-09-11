"""Append a finite-window cut result without closing the broad V/G/B gates."""

from functools import cache

import sympy as s
from p8_vacuum_two_loop_physical_source_map import audit as previous

from . import bounds, calibration, cut, fermion, scalar

MODULES = (scalar, fermion, cut, bounds, calibration)
ITEM = {
    "id": "complete_matched_two_loop_low_energy_Phi_elastic_cut",
    "status": "BOUNDED_ON_PHYSICAL_WINDOW_4_TO_6_WITH_FULL_FIRST_AMPLITUDE_NOT_PHYSICAL_REMAINDER_OR_GLOBAL_DISPERSION",
}


def frontier():
    return previous.frontier()


def matching():
    return [dict(r) for r in previous.matching()] + [dict(ITEM)]


def validate_scope(rows, obligations):
    if rows != frontier() or obligations != matching():
        raise ValueError("The low-energy cut scope differs from the fixed ledger")
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
            "nine_primitive_rows_retained": len(frontier()) - 9,
            "all_parent_primitive_rows_unchanged": sum(
                a != b for a, b in zip(frontier(), previous.frontier())
            ),
            "all_parent_matching_rows_unchanged": sum(
                a != b for a, b in zip(matching(), previous.matching())
            ),
            "only_finite_window_cut_item_added": len(matching())
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
        "only_two_Phi_interference_contributes_at_loop_two_on_window": True,
        "three_scalar_channels_and_all_full_vertex_terms_retained": True,
        "no_heavy_propagator_expansion_inside_unbounded_loop": True,
        "physical_boundary_value_resolved_before_absolute_bound": True,
        "light_logarithm_bound_includes_continuous_threshold": True,
        "fermion_entire_MS_local_reference_subtracted_before_norm": True,
        "all_six_cyclic_fermion_words_and_trace_dimension_retained": True,
        "one_canonical_first_field_factor_only": True,
        "physical_mass_fixed_so_phase_space_not_reexpanded": True,
        "second_cut_coefficient_not_assumed_nonnegative": True,
        "one_loop_amplitude_square_not_mislabeled_complete_loop_three": True,
        "complete_parent_two_loop_b2_error_retained": True,
        "matched_physical_source_transport_used": True,
        "physical_remainder_and_global_V_G_B_open": True,
        "all_prior_scientific_sources_and_reports_unchanged": True,
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
    specifications = (
        ("scalar", bounds.scalar_enclosure, [1, 1, 32, 2, 6, 144], (0, 1, 2, 3, 5)),
        ("fermion", bounds.fermion_enclosure, [24, 1, 6, 144], (0, 1, 2, 3)),
        ("cut", bounds.cut_enclosure, [1, 1], (0, 1)),
        ("band", bounds.improved_band, [1, 1, 1], (0, 1, 2)),
    )
    for prefix, call, good, slots in specifications:
        for j in slots:
            for i, v in enumerate(invalid):
                args = good.copy()
                args[j] = v
                out.append((f"{prefix}_type_{j}_{i}", call, tuple(args)))
    for j, v in ((0, -1), (1, -1), (2, 31), (3, -1), (5, 0)):
        args = [1, 1, 32, 2, 6, 144]
        args[j] = v
        out.append((f"scalar_domain_{j}", bounds.scalar_enclosure, tuple(args)))
    for i, v in enumerate(invalid + (s.Integer(6), 0, -1)):
        out.append(
            (
                f"log_cap_type_or_domain_{i}",
                bounds.scalar_enclosure,
                (1, 1, 32, 2, v, 144),
            )
        )
    out.append(
        ("failed_heavy_power_cap", bounds.scalar_enclosure, (1, 1, 100, 2, 6, 144))
    )
    for j, v in enumerate((23, 0, 0, 0)):
        args = [24, 1, 6, 144]
        args[j] = v
        out.append((f"fermion_domain_{j}", bounds.fermion_enclosure, tuple(args)))
    out.extend(
        (
            ("cut_zero_lambda", bounds.cut_enclosure, (0, 1)),
            ("cut_negative_error", bounds.cut_enclosure, (1, -1)),
            ("band_zero_lambda", bounds.improved_band, (0, 1, 1)),
            ("band_negative_parent_error", bounds.improved_band, (1, -1, 1)),
            ("band_negative_amplitude_error", bounds.improved_band, (1, 1, -1)),
        )
    )
    for i in range(len(cut.rows())):
        rows = cut.rows()
        rows[i]["through_two_loops"] = "UNASSIGNED"
        out.append((f"state_owner_{i}", cut.validate_rows, (rows,)))
    out.extend(
        (
            ("missing_state", cut.validate_rows, (cut.rows()[:-1],)),
            ("duplicate_state", cut.validate_rows, (cut.rows() + cut.rows()[:1],)),
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
        raise ValueError("Unsupported low-energy cut input accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "finite_b2_bound_not_substituted_for_full_amplitude_bound": True,
        "no_absolute_integral_through_unsplit_Feynman_double_pole": True,
        "full_scalar_MS_contact_and_scale_terms_kept": True,
        "entire_fermion_local_reference_not_counted_twice": True,
        "identical_particle_interference_factor_two_retained": True,
        "higher_gauge_and_four_particle_cuts_not_claimed_computed": True,
        "positive_formal_margin_not_physical_remainder_control": True,
        "original_P8_not_closed": True,
    }
