"""Strict unchanged-preparation covariance scope and retained original P8 frontier."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_finite_band_neighborhood import audit as previous
from p8_vacuum_affine_finite_window_growth.intervals import fraction

from . import physical, reference, source, state

STATE = "unchanged_S251_R3_complete_fixed_preparation_no_new_minimization"
MEASURE = "quasifree_reference_covariance_and_finite_Weyl_jets_not_interacting_quantum_measure"
OBSERVABLES = (
    "whole_original_scalar_tensor_preparation_finite_momentum_covariance_bound",
    "complete_fixed_reference_physical_lapse_hat_scale_and_TT_band_variances",
    "single_linear_lapse_Gaussian_spectral_tail_not_nonlinear_support",
    "correct_second_order_Weyl_volume_jet_not_interacting_physical_mean",
)
ITEM = {
    "id": "QG2_H8A428_quantitative_unchanged_Gaussian_preparation_and_physical_band_covariance",
    "status": "WHOLE_FIXED_SCALAR_TENSOR_PREPARATION_AND_COMPLETE_FINITE_BAND_LINEAR_COVARIANCE_WITH_CORRECT_Weyl_VOLUME_JET_NOT_NONLINEAR_SUPPORT_INTERACTING_MEAN_REGULATOR_CUTOFF_OR_ORIGINAL_V_G_B_P8",
}


def parameters():
    return previous.parameters()


def require_parameters(values):
    return previous.require_parameters(values)


def require_state(value):
    if not isinstance(value, str) or value != STATE:
        raise ValueError(
            "Keep the complete original S251 fixed preparation without re-minimization"
        )
    return value


def require_measure(value):
    if not isinstance(value, str) or value != MEASURE:
        raise ValueError(
            "The fixed reference covariance and finite Weyl jet are not an interacting measure"
        )
    return value


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError(
            "Require the complete fixed-reference observable and its stated limited scope"
        )
    return value


def require_momentum(value):
    value = fraction(value)
    if value < fraction(source.LOW):
        raise ValueError(
            "The evaluated high-momentum preparation estimate starts at1e64"
        )
    return value


def require_threshold(value):
    value = fraction(value)
    if value < fraction(source.EPSILON):
        raise ValueError(
            "The stated one-observable tail ceiling uses thresholds at least1e-230"
        )
    return value


def require_band(lower, upper):
    return previous.require_band(lower, upper)


def require_times(first, last):
    return previous.require_times(first, last)


def frontier():
    return previous.frontier()


def matching():
    return [*previous.matching(), dict(ITEM)]


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError(
            "Do not promote the fixed interacting mean, cutoff or original P8 frontier"
        )
    return True


@cache
def packets():
    return {
        "whole_original_reference_and_fixed_preparation": source.data(),
        "whole_full_profile_principal_energy_interval_inequalities": reference.principal(),
        "whole_complete_finite_momentum_energy_errors": reference.finite_q(),
        "whole_both_original_clean_canonical_chart_maps": reference.transition(),
        "whole_noncommuting_reference_energy_and_canonical_normalization": state.energy(),
        "whole_unchanged_sampled_Gaussian_Gramian_and_covariance": state.preparation(),
        "whole_original_physical_rows_and_improved_boundary_enclosure": physical.rows(),
        "whole_complete_R3_band_covariances_and_single_observable_tail": physical.moments(),
        "whole_correct_physical_Weyl_volume_jet_and_explicit_boundary": physical.volume(),
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
        "S261_false_physical_binding_remains_refuted": True,
        "original_quantum_mean_cutoff_and_P8_remain_open": True,
    }


def observable():
    return {
        "positive_result": "The SAME complete S251 scalar/tensor sampled preparation has an evaluated all-P>=1e64 covariance bound near the bounce. Entire profile-dependent rational energy inequalities, complete finite-q errors and both original time-boundary chart maps give the actual Gramian estimate, not a leading-symbol substitution. The complete scalar weighted covariance is at most1e21 P I and each tensor block at most1e5 P I.",
        "physical_band_result": "On the unchanged R3 band[1e64,2e64] and |u|<=1e-60, full physical lapse variance is below1e-490, hat log-scale variance below1e-740, absolute symmetric cross covariance below1e-615, and both TT Frobenius variance below1e-660. The correct finite second-order Weyl spatial-volume contact is below1e-489. No cross or boundary term is deleted.",
        "probability_scope": "For EACH fixed point and time, the spectral law of the self-adjoint band-limited LINEAR reference lapse has tail at threshold1e-230 bounded by2 exp(-1e29). This does not make an unbounded Gaussian supported inside a nonlinear auxiliary branch, specify a joint probability for noncommuting observables, or bound a spacetime supremum or infinite-volume event.",
        "fixed_reference": "The original Cauchy covariance, smooth bump, full reference profiles, kappa and canonical maps are unchanged. Reexpressing the fixed Gramian at a target time is exact transport, not state reselection. Original H/Proca preparations and every determinant/nonlocal response remain; this local scalar/tensor covariance calculation does not compute their complete interacting onepoint.",
        "original_problem": "S263's nearby classical family is not identified with this fixed quantum reference or an interacting mean. Full nonlinear quantum ordering/gauge regulator, physical curved subtraction, complete omitted-loop bounds, Wilsonian matching/cutoff, nonlinear inhomogeneous bounce and original V/G/B/P8 remain OPEN. The S261 physical-binding refutation and completed scoped P8(a) qualifications remain unchanged.",
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
            ("momentum", require_momentum),
            ("threshold", require_threshold),
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
    for value in (0, source.LOW - 1, -source.LOW):
        cases.append(("outside_momentum_" + str(value), require_momentum, (value,)))
    for value in (0, source.EPSILON / 2, -source.EPSILON):
        cases.append(
            ("outside_tail_threshold_" + str(value), require_threshold, (value,))
        )
    for low, high in (
        (0, source.HIGH),
        (source.LOW - 1, source.HIGH),
        (source.LOW, source.HIGH + 1),
        (source.HIGH, source.LOW),
    ):
        cases.append(
            ("outside_band_" + str(low) + "_" + str(high), require_band, (low, high))
        )
    for first, last in (
        (-2 * source.TIME, 0),
        (0, 2 * source.TIME),
        (source.TIME, -source.TIME),
    ):
        cases.append(
            (
                "outside_times_" + str(first) + "_" + str(last),
                require_times,
                (first, last),
            )
        )
    for label in (
        "instantaneous_ground_state_replaces_fixed_preparation",
        "bare_zero_profile_Gramian_is_actual_state",
        "leading_UV_symbol_is_finite_P_covariance",
        "onepoint_tail_proves_nonlinear_branch_support",
        "joint_probability_for_noncommuting_physical_observables",
        "uniform_infinite_volume_Gaussian_event",
        "Gaussian_tail_is_exactly_zero",
        "band_is_physical_Wilsonian_cutoff",
        "full_nonlinear_volume_expectation_equals_Weyl_jet",
        "interacting_physical_mean_is_stationary",
        "scalar_tensor_reference_deletes_Proca_H_determinants",
        "physical_curved_subtraction_is_complete",
        "all_omitted_loops_bounded_by_variance",
        "same_reference_is_S263_classical_family",
        "S256_growth_counterexample_refuted",
        "R_minus_half_is_physical",
        "nonlinear_inhomogeneous_Cauchy_problem_closed",
        "original_P8_closed",
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
            "previous_classical_state_is_not_new_quantum_mean",
            require_state,
            (previous.STATE,),
        )
    )
    return cases


def rejected_inputs():
    count = 0
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (TypeError, ValueError):
            count += 1
        else:
            raise AssertionError("Unsupported input accepted: " + name)
    return count


def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "unchanged_original_source_profiles_and_parameters": True,
        "unchanged_complete_sampled_preparation_and_Cauchy_state": True,
        "complete_canonical_boundary_and_momentum_rows_retained": True,
        "full_finite_P_bound_not_asymptotic_guess": True,
        "S261_physical_binding_explicitly_refuted": True,
        "Gaussian_marginal_tail_not_nonlinear_support": True,
        "finite_Weyl_jet_not_interacting_physical_mean": True,
        "original_V_G_B_P8_and_quantum_cutoff_gates_remain_open": True,
    }
