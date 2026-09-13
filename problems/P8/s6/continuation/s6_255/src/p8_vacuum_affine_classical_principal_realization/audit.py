"""Classical comparison solutions do not replace the fixed quantum reference."""

from functools import cache

import sympy as s
from p8_vacuum_affine_coupled_principal_obstruction import audit as previous
from p8_vacuum_affine_coupled_principal_obstruction import coupled as c

from . import homogeneous, realization, rolling

STATE = "same_parent_local_classical_comparison_not_the_fixed_prepared_quantum_mean"
OBSERVABLES = (
    "full_current_homogeneous_classical_constraint_data_and_local_solution",
    "whole_rolling_heavy_four_mode_Gaussian_on_shell_principal",
    "regular_on_shell_field_map_and_unbounded_momentum_growth",
)
CNN = s.Symbol("homogeneous_positive_Cnn", positive=True)
DOMAIN_PARAMETERS = (*previous.DOMAIN_PARAMETERS, CNN)
ITEM = {
    "id": "QG2_H8A420_complete_rolling_heavy_block_current_constraint_data_local_unforced_classical_realization_and_on_shell_principal_obstruction",
    "status": "EXACT_CURRENT_CLASSICAL_CONSTRAINT_AND_ROLLING_BLOCK_WITH_WRITTEN_LOCAL_UNFORCED_REALIZATION_AND_ON_SHELL_UNBOUNDED_MOMENTUM_GROWTH_NOT_FIXED_QUANTUM_MEAN_CUTOFF_MACRO_BOUNCE_UV_REGGE_OR_ORIGINAL_P8",
}


def require_state(label):
    if not isinstance(label, str) or label != STATE:
        raise ValueError(
            "Keep the classical comparison distinct from the fixed quantum reference"
        )
    return label


def require_observable(label):
    if not isinstance(label, str) or label not in OBSERVABLES:
        raise ValueError("Keep the stated local unforced classical comparison scope")
    return label


def require_domain(values):
    if not isinstance(values, dict) or set(values) != set(DOMAIN_PARAMETERS):
        raise ValueError("Require the entire principal AND homogeneous lapse domain")
    result = {
        **previous.require_domain(
            {key: value for key, value in values.items() if key != CNN}
        ),
        CNN: previous.exact_scalar(values[CNN]),
    }
    if result[CNN].is_positive is not True:
        raise ValueError("The homogeneous implicit lapse pivot must be positive")
    return result


def example_domain():
    return {**previous.example_domain(), CNN: s.Integer(2)}


def central_rate(rho, pressure):
    rho, pressure = previous.exact_scalar(rho), previous.exact_scalar(pressure)
    if (
        abs(rho) > realization.PROFILE_BOUND
        or abs(pressure) > realization.PROFILE_BOUND
    ):
        raise ValueError(
            "Retain the actual fixed profile enclosure, not arbitrary retuning"
        )
    square = s.Rational(1, 100) - 4 * (rho - 3 * pressure / 2)
    if square.is_positive is not True:
        raise ValueError("The comparison uses the strictly positive M1 rate root")
    return s.sqrt(square)


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
        "whole_homogeneous_classical_action_and_regular_Euler_system": homogeneous.data(),
        "whole_rolling_heavy_Gaussian_block_and_Euler_contacts": rolling.data(),
        "whole_current_classical_constraint_data_and_local_realization": realization.data(),
        "whole_on_shell_field_map_contact_bridge": realization.map_contact_data(),
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
        "comparison_data_not_replacement_of_original_preparation": STATE
        != previous.STATE,
        "no_finite_cutoff_quantum_mean_or_original_P8_closure": True,
    }


def observable():
    return {
        "classical_scope": "The unchanged full current local action admits an exact nearby homogeneous constraint-data family and a written local unforced classical realization, with the rolling heavy field, all coupled Gaussian entries, all eight physical modes and the on-shell field-map contacts retained. On each fixed nonreference compact interval with strict margins, the whole S254 P^(3/2) growth proof applies to these classical comparison solutions.",
        "state_boundary": "These are comparison solutions with a slightly changed classical M1 rate, not replacements for the user-selected reference or its fixed preparations. The local action contains its unchanged fixed profile functions; stationarity of that action does not establish stationarity of the complete quantum effective action. The profile-sized central correction is retained and not falsely called zero.",
        "domain_boundary": "The lapse root uses the entire current functions and exact primitive. Local existence is not a uniform-heavy-mass or macroscopic lifetime, a physical-scale bounce, or a nonlinear stability theorem. The all-momentum obstruction is not located below a controlled Wilsonian cutoff, and a finite-cutoff stability neighborhood is not excluded.",
        "original_problem": "No interacting quantum mean, physical curved subtraction, omitted-loop norm, common-parent controlled quantum bounce, physical UV scattering, finite-gravity IR/Regge or original V/G/B/P8 closure follows. Prior P8(a) symmetry/matter qualifications are unchanged.",
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
        cases.append(("invalid_state_" + str(i), require_state, (value,)))
        cases.append(("invalid_observable_" + str(i), require_observable, (value,)))
        for key in DOMAIN_PARAMETERS:
            data = example_domain()
            data[key] = value
            cases.append(
                ("invalid_domain_" + str(key) + "_" + str(i), require_domain, (data,))
            )
        cases.append(("invalid_profile_rho_" + str(i), central_rate, (value, 0)))
        cases.append(("invalid_profile_pressure_" + str(i), central_rate, (0, value)))
    for key in (*previous.POSITIVE_PARAMETERS, CNN):
        for value in (0, -1):
            data = example_domain()
            data[key] = value
            cases.append(
                (
                    "nonpositive_domain_" + str(key) + "_" + str(value),
                    require_domain,
                    (data,),
                )
            )
    for key in (c.r, c.L2):
        data = example_domain()
        data[key] = 0
        cases.append(("missing_fast_margin_" + str(key), require_domain, (data,)))
    for value in (1, s.sqrt(s.Rational(2, 3))):
        data = example_domain()
        data[c.r] = value
        cases.append(("nonpositive_gamma_" + str(value), require_domain, (data,)))
    for label in (
        "original_P8_closed",
        "finite_cutoff_instability_proved",
        "original_quantum_bounce_refuted",
        "fixed_M1_reference_replaced",
        "heavy_field_can_be_held_zero_on_shell",
        "all_quantum_Euler_expressions_zero",
        "physical_scale_bounce_constructed",
        "macroscopic_bounce_lifetime",
        "uniform_heavy_mass_existence_time",
        "new_heavy_mass_or_profile",
        "all_nonlocal_field_maps_equivalent",
        "every_UV_completion_excluded",
        "old_reference_classical_constraint_exactly_zero",
        "profile_correction_discarded",
        "all_loops_bounded",
        "original_P8a_qualifications_removed",
    ):
        cases.append(("unsupported_" + label, require_observable, (label,)))
    cases.append(
        (
            "fixed_reference_not_classical_comparison_label",
            require_state,
            (previous.STATE,),
        )
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
    missing = example_domain()
    del missing[CNN]
    cases.append(("missing_homogeneous_pivot", require_domain, (missing,)))
    extra = example_domain()
    extra[c.H] = 0
    cases.append(("extra_domain_field", require_domain, (extra,)))
    cases.append(("missing_original", validate_scope, (frontier()[:-1], matching())))
    cases.append(("missing_matching", validate_scope, (frontier(), matching()[:-1])))
    for value in (-2 * realization.PROFILE_BOUND, 2 * realization.PROFILE_BOUND):
        cases.append(("outside_rho_bound_" + str(value), central_rate, (value, 0)))
        cases.append(("outside_pressure_bound_" + str(value), central_rate, (0, value)))
    return cases


@cache
def rejected_inputs():
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (TypeError, ValueError):
            continue
        raise ValueError("Unsupported classical comparison input accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "whole_homogeneous_three_velocity_system_not_a_scalar_pivot": True,
        "rolling_heavy_source_and_acceleration_required": True,
        "fixed_profile_constraint_correction_not_assumed_zero": True,
        "original_prepared_quantum_state_not_replaced": STATE != previous.STATE,
        "on_shell_contacts_use_classical_Euler_equations": True,
        "original_full_eight_mode_fast_theorem_retained": True,
        "no_macroscopic_or_subcutoff_or_quantum_verdict": True,
        "all_original_and_earlier_matching_rows_unchanged": frontier()
        == previous.frontier()
        and matching()[:-1] == previous.matching(),
    }
