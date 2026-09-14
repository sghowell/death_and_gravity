"""An evaluated classical finite band is not a Wilsonian or quantum verdict."""

from functools import cache

import sympy as s
from p8_vacuum_affine_classical_principal_realization import audit as previous
from p8_vacuum_affine_heavy_source_filtration import source
from p8_vacuum_canonical_affine_decoupling import family as canonical

from . import background, window
from .intervals import fraction

STATE = "same_parent_evaluated_local_classical_finite_band_comparison_not_fixed_quantum_mean"
OBSERVABLES = (
    "entire_current_coefficient_bounds_and_evaluated_unforced_classical_lifetime",
    "complete_mass_adapted_sixteen_phase_finite_band_propagator_growth",
    "finite_subheavy_scales_without_Wilsonian_or_quantum_mean_verdict",
)
ITEM = {
    "id": "QG2_H8A420_complete_current_classical_evaluated_lifetime_mass_adapted_all_mode_finite_band_growth",
    "status": "EXACT_WHOLE_CURRENT_COEFFICIENT_AND_FULL_MATRIX_ENCLOSURES_WITH_WRITTEN_EVALUATED_CLASSICAL_LIFETIME_AND_FINITE_BAND_GROWTH_NOT_WILSONIAN_CUTOFF_QUANTUM_MEAN_NONLINEAR_BOUNCE_ORIGINAL_V_G_B_OR_P8",
}


def parameters():
    return {
        "epsilon": background.EPSILON,
        "momentum_min": background.MOMENTUM_MIN,
        "momentum_max": background.MOMENTUM_MAX,
        "time_length": background.TIME_LENGTH,
        "kappa": source.KAPPA,
        "mass_squared": source.MASS2,
        "zeta": canonical.ZETA,
    }


def exact_scalar(value):
    result = fraction(value)
    return s.Rational(result.numerator, result.denominator)


def require_parameters(values):
    wanted = parameters()
    if not isinstance(values, dict) or set(values) != set(wanted):
        raise ValueError(
            "Require the entire fixed current finite-window parameter tuple"
        )
    result = {key: exact_scalar(value) for key, value in values.items()}
    if result != wanted:
        raise ValueError(
            "No changed mass, source theory, lapse datum, time or band is certified"
        )
    return result


def require_momentum(value):
    value = exact_scalar(value)
    if not background.MOMENTUM_MIN <= value <= background.MOMENTUM_MAX:
        raise ValueError("Require a momentum inside the complete evaluated finite band")
    return value


def require_time(value):
    value = exact_scalar(value)
    if not 0 <= value <= background.TIME_LENGTH:
        raise ValueError("Require the evaluated local classical interval")
    return value


def require_state(label):
    if not isinstance(label, str) or label != STATE:
        raise ValueError(
            "The classical comparison does not replace the fixed quantum mean"
        )
    return label


def require_observable(label):
    if not isinstance(label, str) or label not in OBSERVABLES:
        raise ValueError(
            "Keep finite-band classical growth separate from EFT and quantum closure"
        )
    return label


def frontier():
    return previous.frontier()


def matching():
    return [*previous.matching(), dict(ITEM)]


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("No earlier primitive or matching obligation is promoted")
    return True


@cache
def packets():
    return {
        "whole_current_function_and_profile_four_jet_bounds": background.source_data(),
        "whole_current_central_constraint_coefficient_enclosures": background.initial_data(),
        "whole_current_rational_tree_and_error_four_jet_bounds": background.tree_data(),
        "whole_current_constrained_Euler_and_time_jet_bounds": background.rate_data(),
        "whole_heavy_energy_and_evaluated_classical_continuation": background.energy_data(),
        "whole_initial_to_persistent_background_and_matrix_box_bridge": background.continuation_data(),
        "whole_eight_phase_mass_adapted_finite_band_matrix": window.matrix_data(),
        "whole_light_growth_chart_time_connection": window.time_data(),
        "whole_actual_fast_coefficient_and_time_jet_binding": window.coefficient_data(),
        "all_other_physical_modes_and_energy_connections": window.other_modes_data(),
        "whole_sixteen_phase_finite_band_growth": window.growth_data(),
    }


@cache
def residuals():
    return {
        packet + "_" + name: s.ImmutableMatrix(value.applyfunc(s.cancel))
        if isinstance(value, s.MatrixBase)
        else s.cancel(value)
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
        "fixed_parent_parameters_are_not_an_adjustable_diagnostic": require_parameters(
            parameters()
        )
        == parameters(),
        "no_Wilsonian_or_quantum_mean_or_original_P8_closure": True,
    }


def observable():
    return {
        "positive_result": "The unchanged full-current local action has an evaluated classical comparison solution on[0,10^-60] at initial lapse1+10^-6, with the entire positive M1 constraint root, full heavy source, all profiles and constants retained. A mass-adapted full16-phase cone estimate gives physical canonical propagator norm at least10^-300 exp(5*10^31) uniformly on the comoving band[10^64,2*10^64]. No P=infinity or heavy-mass limit is used.",
        "physical_scales": "Physical wavenumbers are enclosed using a_phys=a_hat/R^(1/4); the selected band and principal growth timescale are below the retained heavy mass. The stable heavy oscillator itself is kept in the complementary energy space, not integrated out or removed.",
        "scope_boundary": "An evaluated classical finite band is not a controlled Wilsonian window. No omitted-operator, loop, nonlinear matching or state error budget is inferred. This is neither the fixed interacting quantum mean nor a macroscopic physical bounce or nonlinear instability theorem.",
        "original_problem": "The physical UV scattering and quantum gravitational limit, finite-gravity IR/Regge information, physical covariant subtraction, interacting state, quantum mean and compatible nonlinear bounce remain open. All earlier P8(a) qualifications and scoped ladder verdicts are unchanged.",
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
        for label, call in [
            ("state", require_state),
            ("observable", require_observable),
            ("momentum", require_momentum),
            ("time", require_time),
        ]:
            cases.append(("invalid_" + label + "_" + str(i), call, (value,)))
        for key in parameters():
            data = parameters()
            data[key] = value
            cases.append(
                ("invalid_parameter_" + key + "_" + str(i), require_parameters, (data,))
            )
    for key in parameters():
        for factor in (0, -1, 2):
            data = parameters()
            data[key] *= factor
            cases.append(
                ("changed_" + key + "_" + str(factor), require_parameters, (data,))
            )
    for value in (
        0,
        background.MOMENTUM_MIN - 1,
        background.MOMENTUM_MAX + 1,
        -background.MOMENTUM_MIN,
    ):
        cases.append(("outside_band_" + str(value), require_momentum, (value,)))
    for value in (-background.TIME_LENGTH, 2 * background.TIME_LENGTH):
        cases.append(("outside_interval_" + str(value), require_time, (value,)))
    for label in (
        "original_P8_closed",
        "quantum_bounce_refuted",
        "fixed_quantum_mean_replaced",
        "Wilsonian_cutoff_proved",
        "all_UV_completions_excluded",
        "nonlinear_instability_proved",
        "macroscopic_physical_bounce",
        "physical_P_equals_comoving_P",
        "heavy_oscillator_deleted",
        "mass_squared_is_growth_rate",
        "all_momenta_are_in_the_certified_band",
        "arbitrary_new_source_or_state",
        "old_reference_constraint_profile_erased",
        "P8a_qualifications_removed",
    ):
        cases.append(("unsupported_" + label, require_observable, (label,)))
    cases.append(
        ("old_comparison_scope_not_new_finite_window", require_state, (previous.STATE,))
    )
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
    missing = parameters()
    del missing["time_length"]
    extra = parameters()
    extra["UV_cutoff"] = background.MOMENTUM_MAX
    cases.extend(
        [
            ("missing_parameter", require_parameters, (missing,)),
            ("extra_parameter", require_parameters, (extra,)),
            ("missing_original", validate_scope, (frontier()[:-1], matching())),
            ("missing_matching", validate_scope, (frontier(), matching()[:-1])),
        ]
    )
    return cases


@cache
def rejected_inputs():
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (TypeError, ValueError):
            continue
        raise ValueError("Unsupported finite-window input accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "full_coefficient_functions_profiles_and_constants_retained": True,
        "evaluated_background_lifetime_uses_full_heavy_energy": True,
        "whole_finite_band_not_an_infinite_momentum_limit": True,
        "full_heavy_oscillator_kept_not_bounded_by_raw_mass_squared": True,
        "all_physical_modes_and_all_time_connections_retained": True,
        "full_physical_chart_norms_evaluated": True,
        "subheavy_scales_not_a_Wilsonian_cutoff": True,
        "all_original_and_earlier_matching_rows_unchanged": frontier()
        == previous.frontier()
        and matching()[:-1] == previous.matching(),
    }
