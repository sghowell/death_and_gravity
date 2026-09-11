"""A specified free in/out observable, not the complete interacting B state."""

from functools import cache

import sympy as s
from p8_vacuum_clock_transparent_map import audit as previous

from . import calibration, energy, mode, profile, transition

MODULES = (mode, profile, transition, energy, calibration)
ITEM = {
    "id": "free_flat_clock_Dirac_out_energy_with_complete_transition_series_remainder",
    "status": "BOUNDED_FOR_SPECIFIED_QUADRATIC_MODEWISE_IN_OUT_STATE_NOT_INTERACTING_OR_CURVED_BOUNCE_STRESS",
}


def frontier():
    return previous.frontier()


def matching():
    return [dict(row) for row in previous.matching()] + [dict(ITEM)]


def validate_scope(rows, obligations):
    if rows != frontier() or obligations != matching():
        raise ValueError(
            "The free Dirac observable cannot change the fixed scope ledger"
        )
    return True


@cache
def residuals():
    out = {
        mod.__name__.rsplit(".", 1)[-1] + "_" + key: s.simplify(value)
        for mod in MODULES
        for key, value in mod.data()["checks"].items()
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
            "only_free_flat_in_out_state_item_added": len(matching())
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
        **transition.data()["bounds"],
        **calibration.data()["bounds"],
        "only_specified_quadratic_flat_Dirac_operator": True,
        "interacting_gauge_particle_observable_not_inferred": True,
        "integrable_mass_tails_give_modewise_in_out_limits": True,
        "modewise_CAR_probability_preserved": True,
        "complete_odd_transition_Dyson_tail_retained": True,
        "full_time_dependent_frequency_in_transition_phase": True,
        "both_integration_by_parts_boundary_terms_vanish": True,
        "profile_total_variation_bounds_not_sampling_claims": True,
        "zero_momentum_and_constant_mass_production_vanish": True,
        "spin_particle_antiparticle_and_active_color_flavor_counts": True,
        "full_momentum_energy_integral_finite": True,
        "finite_energy_not_promoted_to_Hadamard_or_local_stress_bound": True,
        "kappa_reference_ratio_not_bounce_density_error": True,
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
    for i, value in enumerate(invalid + (-1, s.Rational(1, 10))):
        out.append((f"mixing_type_or_range_{i}", transition.mixing_tail, (value,)))
    good = [1, 10, s.Rational(1, 100), 1]
    for slot in range(4):
        for i, value in enumerate(invalid):
            args = good.copy()
            args[slot] = value
            out.append(
                (f"amplitude_type_{slot}_{i}", transition.amplitude_terms, tuple(args))
            )
    for slot, value in ((0, -1), (1, 0), (2, -1), (3, 0), (2, s.Rational(1, 10))):
        args = good.copy()
        args[slot] = value
        out.append(
            (f"amplitude_domain_{slot}", transition.amplitude_terms, tuple(args))
        )
    good = [10, s.Rational(1, 100), 1, 6]
    for slot in range(3):
        for i, value in enumerate(invalid):
            args = good.copy()
            args[slot] = value
            out.append((f"energy_type_{slot}_{i}", energy.enclosure, tuple(args)))
    for i, value in enumerate(invalid + (s.Integer(6), 0, -1)):
        out.append(
            (
                f"energy_multiplicity_{i}",
                energy.enclosure,
                (10, s.Rational(1, 100), 1, value),
            )
        )
    for slot, value in ((0, 0), (1, -1), (2, 0), (1, s.Rational(1, 10))):
        args = good.copy()
        args[slot] = value
        out.append((f"energy_domain_{slot}_{value}", energy.enclosure, tuple(args)))
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
        raise ValueError("Unsupported free Dirac production input accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "higher_transition_terms_not_set_to_zero": True,
        "modewise_evolution_not_global_infinite_volume_Fock_unitary": True,
        "free_out_particle_energy_not_transient_renormalized_stress": True,
        "all_transition_orders_not_all_Feynman_loops": True,
        "colored_free_modes_not_interacting_gauge_asymptotic_states": True,
        "reference_scale_ratio_not_bounce_backreaction_bound": True,
        "finite_energy_not_Hadamard_or_curved_state_certificate": True,
        "original_P8_not_closed": True,
    }
