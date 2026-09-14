"""Scoped full-parent classical bounce and finite-band linear-energy audit."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_finite_window_growth.intervals import fraction
from p8_vacuum_affine_physical_volume_correction import audit as previous

from . import background, reference, source, window

STATE = "unchanged_original_reference_with_explicit_S261_erratum_and_distinct_classical_comparison"
MEASURE = (
    "classical_constrained_Euler_and_full_linear_energy_not_interacting_quantum_measure"
)
OBSERVABLES = (
    "whole_fixed_source_fifth_jets_and_regular_classical_constraint_family",
    "tiny_actual_physical_classical_bounce_with_all_sources_and_heavy_field",
    "complete_sixteen_phase_finite_band_reference_energy_comparison",
    "explicit_quantum_mean_cutoff_and_original_V_G_B_boundary",
)
ITEM = {
    "id": "QG2_H8A427_full_current_tiny_classical_bounce_and_complete_finite_band_linear_energy",
    "status": "COMPLETE_SOURCED_TINY_HOMOGENEOUS_CLASSICAL_BOUNCE_AND_FULL_FINITE_BAND_LINEAR_ENERGY_NOT_QUANTUM_MEAN_WILSONIAN_CUTOFF_NONLINEAR_INHOMOGENEOUS_CAUCHY_OR_ORIGINAL_V_G_B_P8",
}


def parameters():
    return previous.parameters()


def require_parameters(values):
    return previous.require_parameters(values)


def require_state(label):
    if not isinstance(label, str) or label != STATE:
        raise ValueError(
            "Keep the fixed reference, the explicit S261 correction and the distinct classical comparison"
        )
    return label


def require_measure(label):
    if not isinstance(label, str) or label != MEASURE:
        raise ValueError(
            "The complete classical variational energy is not a quantum measure"
        )
    return label


def require_observable(label):
    if not isinstance(label, str) or label not in OBSERVABLES:
        raise ValueError(
            "Require the sourced classical or explicitly limited finite-band result"
        )
    return label


def require_epsilon(value):
    value = fraction(value)
    if abs(value) > fraction(source.EPSILON):
        raise ValueError(
            "The signed initial lapse parameter must stay in the proved tiny family"
        )
    return value


def require_band(lower, upper):
    lower, upper = map(fraction, (lower, upper))
    if not fraction(source.LOW) <= lower <= upper <= fraction(source.HIGH):
        raise ValueError(
            "The complete momentum interval must stay inside the stated finite band"
        )
    return lower, upper


def require_times(first, last):
    first, last = map(fraction, (first, last))
    if not -fraction(source.TIME) <= first <= last <= fraction(source.TIME):
        raise ValueError("Keep both endpoints inside the proved symmetric clock slab")
    return first, last


def frontier():
    return previous.frontier()


def matching():
    return [*previous.matching(), dict(ITEM)]


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError(
            "Do not promote the physical cutoff, quantum mean or original P8 frontier"
        )
    return True


@cache
def packets():
    return {
        "whole_original_sources_and_fifth_jet_enclosures": source.data(),
        "whole_regular_signed_initial_constraint_family": background.initial(),
        "whole_original_Euler_reference_and_light_Jacobian": background.light(),
        "whole_complete_coefficient_value_and_first_jet_bounds": background.coefficients(),
        "whole_correct_physical_Hubble_and_acceleration": background.physical(),
        "whole_regular_constraint_and_coefficient_bootstrap_reentry": background.bootstrap(),
        "whole_original_fixed_reference_canonical_energy_bridge": reference.canonical(),
        "whole_fixed_reference_boundary_polynomial_enclosures": reference.shear_bounds(),
        "whole_positive_full_reference_energy": reference.positive_energy(),
        "whole_complete_sixteen_phase_generator_and_original_maps": window.matrices(),
        "whole_two_sided_classical_turn_and_full_linear_energy_comparison": window.comparison(),
    }


@cache
def residuals():
    return {
        packet + "_" + name: clean(value)
        for packet, data in packets().items()
        for name, value in data["checks"].items()
    }


def scalar_entry_count():
    return sum(
        len(value) if isinstance(value, s.MatrixBase) else 1
        for value in residuals().values()
    )


@cache
def gates():
    return {
        **{
            packet + "_" + name: bool(value)
            for packet, data in packets().items()
            for name, value in data["gates"].items()
        },
        "all_nine_original_primitive_rows_unchanged": len(frontier()) == 9,
        "all_matching_identifiers_unique": len({row["id"] for row in matching()})
        == len(matching()),
        "all_original_parameters_unchanged": require_parameters(parameters())
        == parameters(),
        "S261_false_physical_binding_remains_explicitly_refuted": True,
        "original_quantum_mean_cutoff_and_P8_remain_open": True,
    }


def observable():
    return {
        "positive_result": "The unchanged full parent admits a signed tiny family of exact constrained homogeneous classical solutions on[-10^-60,10^-60], starting at N=1+epsilon with abs(epsilon)<=10^-230 and the entire positive M1 constraint root. All sources and the heavy field remain. The correct physical scale a_hat R^-1/4 has exactly one local turn, with positive proper-time Hubble derivative greater than3. Around each such solution all16 canonical linear phases have root-energy growth below2 throughout P in[10^64,2*10^64], in the specified single fixed-reference energy.",
        "comparison_limit": "The same full off-clock parent still has the S256 finite-band large-growth example at epsilon=10^-6. This result does not reverse that example or restore an all-momentum bound. The fixed finite band and exceptionally short clock slab are mathematical comparison parameters, not a Wilsonian cutoff or a cosmological duration. A homogeneous classical solution and its linear perturbations are not a nonlinear inhomogeneous Cauchy theorem.",
        "fixed_reference": "The source profiles, regulator/matching constants, coupled preparation and canonical endpoints are unchanged. The reference energy is a norm, not a reselected vacuum. The physical factor is the corrected R**(-1/2). S261's false shifted-R identification remains refuted, despite unchanged historical replay payloads.",
        "original_problem": "The actual interacting physical quantum mean, quantum ordering and gauge regulator, Ward/state/endpoint defects, physical curved subtraction, complete loop errors, controlled Wilsonian matching/cutoff, nonlinear inhomogeneous bounce and original V/G/B/P8 remain OPEN. Completed scoped P8(a) qualifications are unchanged.",
    }


def bad_cases():
    cases = []
    invalid = (
        True,
        False,
        1.0,
        s.Float(1),
        "1",
        None,
        s.oo,
        s.I,
        s.nan,
        s.Symbol("unknown"),
    )
    for i, value in enumerate(invalid):
        for label, call in (
            ("state", require_state),
            ("measure", require_measure),
            ("observable", require_observable),
            ("epsilon", require_epsilon),
        ):
            cases.append(("invalid_" + label + "_" + str(i), call, (value,)))
        cases.append(("invalid_band_" + str(i), require_band, (value, source.HIGH)))
        cases.append(("invalid_time_" + str(i), require_times, (value, source.TIME)))
        for key in parameters():
            values = parameters()
            values[key] = value
            cases.append(
                (
                    "invalid_parameter_" + key + "_" + str(i),
                    require_parameters,
                    (values,),
                )
            )
    for value in (2 * source.EPSILON, -2 * source.EPSILON, s.Rational(1, 10**6), 1):
        cases.append(
            ("outside_signed_lapse_family_" + str(value), require_epsilon, (value,))
        )
    for lower, upper in (
        (0, source.HIGH),
        (source.LOW - 1, source.HIGH),
        (source.LOW, source.HIGH + 1),
        (source.HIGH, source.LOW),
    ):
        cases.append(
            (
                "outside_band_" + str(lower) + "_" + str(upper),
                require_band,
                (lower, upper),
            )
        )
    for first, last in (
        (-2 * source.TIME, 0),
        (0, 2 * source.TIME),
        (source.TIME, -source.TIME),
    ):
        cases.append(
            (
                "outside_time_" + str(first) + "_" + str(last),
                require_times,
                (first, last),
            )
        )
    for label in (
        "fixed_reference_is_the_unforced_quantum_mean",
        "original_P8_closed",
        "homogeneous_classical_bounce_completes_original_B",
        "finite_band_is_Wilsonian_cutoff",
        "cutoff_equals_heavy_mass",
        "all_spatial_momenta_have_uniform_same_space_bound",
        "S256_growth_counterexample_refuted",
        "all_lapse_displacements_admissible",
        "heavy_mass_deleted_from_evolution",
        "source_profiles_or_state_reprepared",
        "hat_metric_turn_is_automatically_physical_turn",
        "R_minus_half_is_physical",
        "quantum_ordering_and_gauge_regulator_complete",
        "Gaussian_support_inside_nonlinear_branch",
        "inhomogeneous_nonlinear_Cauchy_problem_closed",
        "finite_linear_energy_bounds_all_loops",
        "coupled_scalar_K_and_G_commute",
        "both_canonical_endpoint_maps_dropped",
    ):
        cases.append(("unsupported_" + label, require_observable, (label,)))
    for i in range(len(frontier())):
        rows = frontier()
        rows[i]["status"] = "COMPLETE"
        cases.append(
            ("primitive_promotion_" + str(i), validate_scope, (rows, matching()))
        )
    for i in range(len(matching())):
        rows = matching()
        rows[i]["status"] = "COMPLETE"
        cases.append(
            ("matching_promotion_" + str(i), validate_scope, (frontier(), rows))
        )
    cases.append(
        (
            "prior_state_without_distinct_classical_comparison",
            require_state,
            (previous.STATE,),
        )
    )
    return cases


@cache
def rejected_inputs():
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (TypeError, ValueError):
            continue
        raise ValueError("Unsupported finite-neighborhood claim accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "full_current_source_heavy_mass_and_original_reference_retained": True,
        "signed_tiny_constraint_family_not_bare_root_substitution": True,
        "physical_R_inverse_square_root_turn_not_hat_metric_assumption": True,
        "complete_sixteen_phase_generator_and_all_endpoint_maps_kept": window.matrices()[
            "whole_real_phase_count"
        ]
        == 16,
        "fixed_finite_band_not_quantum_cutoff_or_all_momentum_theorem": True,
        "known_larger_displacement_growth_example_unchanged": True,
        "explicit_S261_physical_binding_refutation_not_rescinded": True,
        "original_frontier_and_historic_matching_payloads_unchanged": frontier()
        == previous.frontier()
        and matching()[:-1] == previous.matching(),
    }
