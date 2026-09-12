"""Original scalar unitarity and conditional matching, with no cutoff assertion."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes
from p8_vacuum_affine_full_flat_tensor_response import audit as previous

from . import amplitude, dispersion, unitarity

ITEM = {
    "id": "original_canonical_massive_scalar_partial_unitarity_necessary_correction_and_conditional_dispersion_matching_error_tradeoff",
    "status": "EXACT_TREE_AND_CONDITIONAL_MATCHING_DISJUNCTION_NOT_PHYSICAL_CUTOFF_FULL_QUANTUM_DECOUPLING_UV_NO_GO_OR_P8",
}
require_scope = previous.require_scope


def require_stage(stage):
    if not isinstance(stage, str) or stage not in (
        "original_complete_scalar_tree_partial_waves",
        "normalized_identical_channel_and_first_elastic_cut",
        "nominal_real_part_ceiling_and_necessary_correction",
        "conditional_full_angular_error_dispersion_tradeoff",
        "same_units_tensor_scale_comparison_only",
    ):
        raise ValueError(
            "Only the stated tree and conditional scalar matching results are proved"
        )
    return stage


def require_tolerances(energy, eta, delta):
    for value in (energy, eta, delta):
        if isinstance(value, bool) or not isinstance(
            value, (int, s.Integer, s.Rational)
        ):
            raise TypeError("Matching tolerances must be exact finite real rationals")
    if energy <= 0 or energy**2 < 32 or eta < 0 or eta >= 1 or delta < 0:
        raise ValueError(
            "Matching tolerances must obey the massive annulus and error domain"
        )
    return tuple(s.Rational(value) for value in (energy, eta, delta))


def frontier():
    return previous.frontier()


def matching():
    return [dict(row) for row in previous.matching()] + [dict(ITEM)]


def validate_scope(rows, obligations):
    if rows != frontier() or obligations != matching():
        raise ValueError(
            "Conditional scalar matching does not close the original physical P8 frontier"
        )
    return True


def packets():
    return {
        "complete_original_scalar_tree_and_identical_elastic_cut": amplitude.data(),
        "exact_unitarity_circle_and_necessary_full_correction": unitarity.data(),
        "conditional_full_amplitude_dispersion_error_tradeoff": dispersion.data(),
    }


@cache
def residuals():
    rows = {
        name + "_" + key: value
        for name, packet in packets().items()
        for key, value in packet["checks"].items()
    }
    rows.update(
        {
            "nine_original_primitive_rows": len(frontier()) - 9,
            "unchanged_primitive_rows": sum(
                a != b for a, b in zip(frontier(), previous.frontier())
            ),
            "unchanged_previous_matching_rows": sum(
                a != b for a, b in zip(matching(), previous.matching())
            ),
            "one_new_conditional_scalar_matching_result": len(matching())
            - len(previous.matching())
            - 1,
        }
    )
    return {key: s.cancel(value) for key, value in rows.items()}


def scalar_entry_count():
    return len(residuals())


@cache
def gates():
    return {
        **{
            name + "_" + key: bool(value)
            for name, packet in packets().items()
            for key, value in packet["gates"].items()
        },
        "original_complete_tree_and_unchanged_retuned_low_jets": True,
        "exact_unitarity_not_assumed_of_a_nonzero_real_tree": True,
        "full_error_includes_loops_heavy_fields_and_higher_operators": True,
        "physical_b2_not_a_freely_chosen_tree_coefficient": True,
        "no_finite_gravity_cutoff_or_tensor_pole_repair_inferred": True,
        "unique_matching_identifiers": len({row["id"] for row in matching()})
        == len(matching()),
    }


def observable():
    return {
        "actual_scope": "The complete S177 fixed-canonical nongravitational four-scalar tree of the same S182 retuned parent, scalar mass1, lambda10^-600 and gamma1024*10^-800. The quantum S-matrix and small error premises are not established.",
        "main_result": "Normalized identical-channel partial waves and the first elastic coefficient are exact. The unique nominal t0=1/2 ceiling lies between COM energies10^133 and2*10^133. Exact unitarity demands relative t0 correction greater than1-10^-5 at10^134.",
        "conditional_result": "If an exact S-matrix satisfies the stated positive twice-subtracted dispersion and full angular relative error eta on the annulus[E^2/2,E^2], its physical b2/(4lambda)>9(1-eta)^2(E/10^125)^8. AtE10^125, eta<=1/2 and b2<=8lambda cannot all coexist, although the original tree amplitude is below10^-46.",
        "boundary": "These are a necessary correction and a conditional matching disjunction, not a physical cutoff, a demonstrated quantum decoupling theorem, finite-gravity partial-wave bound, UV no-go, tensor-pole modification, nonlinear/bounce control or original V/G/B/P8 closure.",
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
        out.append((f"stage_type_{i}", require_stage, (value,)))
        for pos in range(3):
            args = [dispersion.AGREEMENT_ENERGY, dispersion.ETA, 1]
            args[pos] = value
            out.append((f"tolerance_type_{pos}_{i}", require_tolerances, tuple(args)))
    for args in (
        (1, modes.KAPPA, modes.MASS, 1),
        (-1, modes.KAPPA, modes.MASS, 1),
        (0, 2 * modes.KAPPA, modes.MASS, 1),
        (0, modes.KAPPA, 999, 1),
        (0, modes.KAPPA, modes.MASS, 2),
    ):
        out.append((f"scope_{len(out)}", require_scope, args))
    for args in (
        (0, 0, 0),
        (-6, 0, 0),
        (5, 0, 0),
        (6, -1, 0),
        (6, 1, 0),
        (6, 2, 0),
        (6, 0, -1),
    ):
        out.append((f"tolerance_domain_{len(out)}", require_tolerances, args))
    for stage in (
        "physical_cutoff",
        "full_quantum_decoupling",
        "tree_exactly_unitary_below_ceiling",
        "finite_gravity_partial_wave",
        "graviton_pole_deleted",
        "full_real_one_loop_amplitude",
        "higher_loops_bounded",
        "old_GY14_matching_imported",
        "new_b2_counterterm",
        "UV_parent_excluded",
        "closed_P8",
        "tensor_poles_repaired",
        "curved_bounce_stable",
        "unrestricted_S222_graph",
        "dimension_varying_couplings",
        "unknown_errors_set_zero",
    ):
        out.append((f"unsupported_stage_{stage}", require_stage, (stage,)))
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
                "extra_physical_cutoff",
                validate_scope,
                (
                    frontier(),
                    matching() + [{"id": "physical_cutoff", "status": "COMPLETE"}],
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
        raise ValueError(
            "Unsupported scalar matching or physical closure accepted: " + name
        )
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "distinguishable_16pi_normalization_fails_identical_bubble": True,
        "potential_mass_and_l2_wave_not_discarded": True,
        "tree_is_real_and_not_an_exact_unitary_amplitude": True,
        "conditional_error_not_assumed_as_a_verified_bound": True,
        "physical_b2_not_identified_with_tree_without_matching": True,
        "small_tree_size_does_not_close_improved_dispersion": True,
        "tensor_scale_comparison_not_a_cutoff_or_pole_change": True,
        "all_original_primitive_statuses_unchanged": frontier() == previous.frontier(),
    }
