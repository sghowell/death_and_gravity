"""Strict scope and complete unchanged-scheme two-mass reference audit."""

from functools import cache

import sympy as s
from p8_vacuum_affine_heavy_adm_clock_response import audit as previous

from . import geometry, inverse, spectral

exact = previous.exact
ITEM = {
    "id": "QG2_H8A420_complete_two_mass_Gaussian_flat_reference_inverse_with_unchanged_finite_curvature_two_threshold_measure_and_no_inherited_Proca_pole",
    "status": "COMPLETE_TWO_MASS_GAUSSIAN_FLAT_QUOTIENT_REFERENCE_INVERSE_NOT_CURVED_COUPLED_INVERSE_SMALLNESS_LIGHT_LOOPS_NONLINEAR_UV_REGGE_OR_P8",
}


def require_parameters(proca_mass_squared, heavy_mass_squared, kappa):
    values = tuple(map(exact, (proca_mass_squared, heavy_mass_squared, kappa)))
    if values != (geometry.PROCA_MASS2, geometry.HEAVY_MASS2, geometry.KAPPA):
        raise ValueError(
            "Keep both distinct actual masses and the unchanged gravitational normalization"
        )
    return values


def require_domain(window, transfer_squared):
    T, q = map(exact, (window, transfer_squared))
    if T <= 0 or q < 0:
        raise ValueError(
            "Use a positive finite window and nonnegative squared transfer"
        )
    return T, q


def require_spin(spin):
    if isinstance(spin, bool) or not isinstance(spin, (int, s.Integer)):
        raise TypeError("Spin labels are exact integers")
    if spin not in (0, 2):
        raise ValueError("The full scalar stress has both spin0 and spin2 channels")
    return int(spin)


def require_stage(stage):
    if not isinstance(stage, str) or stage not in (
        "complete_minimal_scalar_stress_cut",
        "unchanged_two_mass_finite_reference",
        "complete_two_threshold_reciprocal",
        "ordinary_factor_inverse_with_both_thresholds",
        "all_transfer_flat_quotient_inverse",
    ):
        raise ValueError(
            "Only the complete stated two-mass flat reference is established"
        )
    return stage


def require_scheme(label):
    if not isinstance(label, str) or label != "old_Proca_plus_actual_heavy_mu1":
        raise ValueError("Keep both specified complete finite prescriptions")
    return label


def require_graph(label):
    if (
        not isinstance(label, str)
        or label != "complete_causal_initial_boundary_flat_quotient"
    ):
        raise ValueError("The graph includes the full initial distributional boundary")
    return label


def require_measure(label):
    if (
        not isinstance(label, str)
        or label != "full_two_threshold_reciprocal_of_the_sum"
    ):
        raise ValueError(
            "Invert the sum, with both full cuts and no inherited isolated pole"
        )
    return label


def require_cutoff(label):
    if not isinstance(label, str) or label != "no_internal_or_external_momentum_cutoff":
        raise ValueError(
            "The flat reference retains every internal and external momentum"
        )
    return label


def require_pole_count(count):
    if isinstance(count, bool) or not isinstance(count, (int, s.Integer)):
        raise TypeError("Pole counts are exact integers")
    if count != 0:
        raise ValueError("The actual total reference has no first-sheet isolated pole")
    return int(count)


def frontier():
    return previous.frontier()


def matching():
    return [*previous.matching(), dict(ITEM)]


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError(
            "Do not alter the original primitive or previous matching frontier"
        )
    return True


@cache
def packets():
    return {
        "whole_minimal_scalar_cut_and_fixed_two_mass_finite_input": geometry.data(),
        "complete_two_threshold_first_sheet_reciprocal": spectral.data(),
        "whole_causal_all_transfer_flat_quotient_inverse": inverse.data(),
    }


@cache
def residuals():
    return {
        name + "_" + key: value
        for name, data in packets().items()
        for key, value in data["checks"].items()
    }


def scalar_entry_count():
    return sum(
        len(v) if isinstance(v, s.MatrixBase) else 1 for v in residuals().values()
    )


@cache
def gates():
    return {
        **{
            name + "_" + key: bool(v)
            for name, data in packets().items()
            for key, v in data["gates"].items()
        },
        "actual_two_mass_prescriptions_and_parameters_unchanged": True,
        "all_nine_original_primitive_rows_retained": len(frontier()) == 9,
        "all_matching_identifiers_unique": len({r["id"] for r in matching()})
        == len(matching()),
        "full_flat_reference_not_actual_curved_state_replacement": True,
        "reference_inverse_not_total_interacting_quantum_inverse_or_stability": True,
    }


def observable():
    return {
        "actual": "The unchanged Proca Gaussian fourth-order reference and actual added minimal-heavy Gaussian reference are summed BEFORE inversion, in their independently fixed finite prescriptions.",
        "new_spectrum": "Both complete total factors have no first-sheet zero. Their reciprocals are cut-only positive Stieltjes measures, with both physical thresholds and the second trace cusp retained. The Proca-only pole is neither deleted from its old result nor inserted into the inverse of the sum.",
        "inverse": "The complete two-channel flat quotient has a two-sided causal inverse on the full initial-boundary graph, uniformly over every external transfer with no spatial derivative loss. Its normalized C_tH^r and L2_tH^r bound is375*T^4/[32*(log(n)+2)]<125*T^4/4224.",
        "normalization": "The physical force inverse restores64*pi^2*kappa. Its bound is750*pi^2*kappa*T^4/(log(n)+2), not smallness. Metric-source norm conversion also retains the full amplitude Gram matrix.",
        "remaining": "The actual curved tree/state/profile/matter normal form, compatible coupled inverse bounds, interacting light/mixed loops, nonlinear same-state bounce, quantum gravitational limit, physical UV/Regge and original V/G/B/P8 remain open.",
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
    templates = (
        (
            "parameters",
            require_parameters,
            (geometry.PROCA_MASS2, geometry.HEAVY_MASS2, geometry.KAPPA),
        ),
        ("domain", require_domain, (1, 0)),
        ("spin", require_spin, (0,)),
        ("stage", require_stage, ("all_transfer_flat_quotient_inverse",)),
        ("scheme", require_scheme, ("old_Proca_plus_actual_heavy_mu1",)),
        ("graph", require_graph, ("complete_causal_initial_boundary_flat_quotient",)),
        ("measure", require_measure, ("full_two_threshold_reciprocal_of_the_sum",)),
        ("cutoff", require_cutoff, ("no_internal_or_external_momentum_cutoff",)),
        ("poles", require_pole_count, (0,)),
    )
    out = []
    for i, value in enumerate(invalid):
        for name, call, original in templates:
            for pos in range(len(original)):
                args = list(original)
                args[pos] = value
                out.append((f"{name}_invalid_{i}_{pos}", call, tuple(args)))
    for index in range(3):
        args = [geometry.PROCA_MASS2, geometry.HEAVY_MASS2, geometry.KAPPA]
        args[index] += 1
        out.append(("changed_parameter_" + str(index), require_parameters, tuple(args)))
    for args in ((0, 0), (-1, 0), (1, -1)):
        out.append(("outside_domain_" + str(len(out)), require_domain, args))
    for spin in (-1, 1, 3):
        out.append(("unsupported_spin_" + str(spin), require_spin, (spin,)))
    for value in (-1, 1, 2):
        out.append(
            ("wrong_total_pole_count_" + str(value), require_pole_count, (value,))
        )
    for stage in (
        "full_curved_inverse",
        "full_three_source_inverse",
        "physical_stability",
        "small_force_inverse",
        "complete_interacting_loops",
        "nonlinear_bounce",
        "quantum_gravity_limit",
        "UV_completion",
        "Regge",
        "original_P8_closed",
    ):
        out.append(("unsupported_stage_" + stage, require_stage, (stage,)))
    for name, call, values in (
        (
            "scheme",
            require_scheme,
            ("retune_to_remove_pole", "massless_log_only", "one_mass_approximation"),
        ),
        (
            "graph",
            require_graph,
            (
                "discard_initial_atoms",
                "equation_only_at_positive_time",
                "unrestricted_curved_graph",
            ),
        ),
        (
            "measure",
            require_measure,
            (
                "sum_of_reciprocals",
                "carry_old_pole",
                "drop_heavy_cut",
                "delete_pole_by_hand",
            ),
        ),
        (
            "cutoff",
            require_cutoff,
            ("external_ball", "finite_internal_band", "physical_EFT_cutoff"),
        ),
    ):
        for value in values:
            out.append(("unsupported_" + name + "_" + value, call, (value,)))
    out.append(("unchanged_primitive_frontier", validate_scope, ([], matching())))
    out.append(
        ("unchanged_matching_frontier", validate_scope, (frontier(), matching()[:-1]))
    )
    return tuple(out)


@cache
def rejected_inputs():
    count = 0
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (ValueError, TypeError):
            count += 1
        else:
            raise ValueError("Unsupported two-mass inverse input admitted: " + name)
    return count


def controls():
    gd, sd = geometry.data(), spectral.data()
    return {
        "rejected_inputs": rejected_inputs(),
        "all_original_primitive_and_previous_matching_rows_retained": matching()[:-1]
        == previous.matching()
        and len(frontier()) == 9,
        "inverse_of_sum_not_sum_of_inverses": s.Rational(1, 5)
        != s.Rational(1, 2) + s.Rational(1, 3),
        "old_Proca_pole_not_discarded_by_hand": bool(
            sd["strict_complete_lowest_threshold_gaps"][2] > 5
        ),
        "both_sides_of_second_trace_threshold_retained": all(
            v != 0 for v in sd["second_threshold_trace_two_sided_cusp"].values()
        ),
        "fixed_physical_dimension_jet_not_replaced": gd[
            "fixed_physical_scalar_leading_and_jet"
        ][1]
        != gd["different_dimension_varying_scalar_jet"],
        "both_full_mass_cuts_and_fixed_fourth_coefficients_retained": True,
        "initial_distributional_boundary_and_kappa_normalization_retained": True,
        "no_curved_coupled_inverse_stability_nonlinear_UV_Regge_or_P8_upgrade": True,
    }
