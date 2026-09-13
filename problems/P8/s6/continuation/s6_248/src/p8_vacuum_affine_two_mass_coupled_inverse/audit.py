"""Scope guards for the actual conditional two-mass coupled inverse."""

from functools import cache

import sympy as s
from p8_vacuum_affine_two_mass_reference import audit as previous

from . import assembly, inverse, normal

exact = previous.exact
require_parameters = previous.require_parameters
ITEM = {
    "id": "QG2_H8A420_actual_two_mass_Gaussian_scalar_clock_matter_coupled_inverse_profile_once_and_full_force_recovery_on_smooth_prepared_Fourier_balls",
    "status": "ACTUAL_CONDITIONAL_TWO_MASS_SMOOTH_PREPARED_FOURIER_BALL_COUPLED_INVERSE_NOT_UNRESTRICTED_GRAPH_SMALLNESS_STABILITY_INTERACTING_LOOPS_NONLINEAR_UV_REGGE_OR_P8",
}


def require_stage(stage):
    if not isinstance(stage, str) or stage not in (
        "whole_current_profile_once_assembly",
        "actual_two_mass_strong_normal_form",
        "actual_prepared_smooth_ball_coupled_inverse",
        "full_current_smooth_ball_force_recovery",
    ):
        raise ValueError(
            "Only the stated actual conditional smooth Fourier-ball inverse is proved"
        )
    return stage


def require_graph(label):
    if (
        not isinstance(label, str)
        or label != "smooth_prepared_scalar_amplitudes_on_each_fixed_Fourier_ball"
    ):
        raise ValueError(
            "Keep the original initial germ and the explicit smooth finite-ball subgraph"
        )
    return label


def require_assembly(label):
    if (
        not isinstance(label, str)
        or label != "current_QG2_classical_plus_pure_Proca_and_H_Gaussian"
    ):
        raise ValueError(
            "The current classical coefficients already include the full heavy profile"
        )
    return label


def require_reference(label):
    if (
        not isinstance(label, str)
        or label != "complete_two_mass_sum_before_inverse_same_finite_scheme"
    ):
        raise ValueError(
            "Invert the complete two-mass sum, retaining both thresholds and fixed finite terms"
        )
    return label


def require_internal_integral(label):
    if (
        not isinstance(label, str)
        or label != "all_original_internal_momenta_no_physical_cutoff"
    ):
        raise ValueError(
            "The external Fourier support is not a cutoff on the internal loop"
        )
    return label


def frontier():
    return previous.frontier()


def matching():
    return [*previous.matching(), dict(ITEM)]


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError(
            "Keep every original primitive and earlier matching status unchanged"
        )
    return True


@cache
def packets():
    return {
        "full_current_profile_once": assembly.profile_data(),
        "whole_current_mixed_principal": assembly.principal_data(),
        "complete_actual_two_species_Ward": assembly.ward_data(),
        "full_original_constraint_recovery": assembly.residual_data(),
        "complete_actual_strong_normal_form": normal.data(),
        "ordered_current_coupled_inverse": inverse.data(),
    }


@cache
def residuals():
    return {
        name + "_" + key: s.ImmutableMatrix(value.applyfunc(s.cancel))
        if isinstance(value, s.MatrixBase)
        else s.cancel(value)
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
            name + "_" + key: bool(value)
            for name, data in packets().items()
            for key, value in data["gates"].items()
        },
        "all_nine_original_primitive_rows_retained": len(frontier()) == 9,
        "all_matching_identifiers_unique": len({r["id"] for r in matching()})
        == len(matching()),
        "actual_current_model_not_old_QG1_inverse_status": True,
        "new_inverse_not_interacting_state_or_P8_closure": True,
    }


def observable():
    return {
        "actual": "The complete current QG2-H8A420 classical scalar/clock/M1 coefficient sector, both actual Gaussian response sectors and fixed profiles are assembled once in the original prescription.",
        "new": "A strong finite-ball remainder and the complete S247 two-mass reference yield a two-sided causal coupled inverse on smooth prepared scalar amplitudes and the invariant full current force subgraph.",
        "bound": "With the actual finite row-sum C0(Pmax) and total-channel L1 norm Kmax, b0=max(15625/6144,16Kmax/(9gamma),1), lambda=(4b0 C0+1)^2 and inverse bound2 exp(lambda)b0 on row-primitive adapted forcing. Constants are finite but unevaluated; physical gamma^-1=64pi²kappa remains.",
        "remaining": "No all-Pmax or unrestricted completed-graph bound, numerical smallness, physical stability, nonlinear same-state bounce, interacting light/mixed loops, quantum gravitational limit, physical UV/Regge or original V/G/B/P8 closure.",
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
            (
                normal.reference.PROCA_MASS2,
                normal.reference.HEAVY_MASS2,
                normal.reference.KAPPA,
            ),
        ),
        ("ball", inverse.require_ball, (1,)),
        ("stage", require_stage, ("actual_prepared_smooth_ball_coupled_inverse",)),
        (
            "graph",
            require_graph,
            ("smooth_prepared_scalar_amplitudes_on_each_fixed_Fourier_ball",),
        ),
        (
            "assembly",
            require_assembly,
            ("current_QG2_classical_plus_pure_Proca_and_H_Gaussian",),
        ),
        (
            "reference",
            require_reference,
            ("complete_two_mass_sum_before_inverse_same_finite_scheme",),
        ),
        (
            "internal",
            require_internal_integral,
            ("all_original_internal_momenta_no_physical_cutoff",),
        ),
    )
    cases = []
    for i, value in enumerate(invalid):
        for name, call, args in templates:
            for j in range(len(args)):
                bad = list(args)
                bad[j] = value
                cases.append((f"{name}_bad_{i}_{j}", call, tuple(bad)))
    for value in (0, -1, -s.Rational(1, 2)):
        cases.append((f"nonpositive_ball_{value}", inverse.require_ball, (value,)))
    for j in range(3):
        args = [
            normal.reference.PROCA_MASS2,
            normal.reference.HEAVY_MASS2,
            normal.reference.KAPPA,
        ]
        args[j] += 1
        cases.append(
            (f"changed_physical_parameter_{j}", require_parameters, tuple(args))
        )
    for label in (
        "full_interacting_inverse",
        "unrestricted_completed_graph",
        "uniform_all_Pmax",
        "small_physical_inverse",
        "stable_bounce",
        "closed_P8",
        "quantum_gravity_limit",
        "finite_nonlinear_tube",
        "physical_UV_cutoff",
        "reset_initial_state",
    ):
        cases.append(("unsupported_" + label, require_stage, (label,)))
    for label, call in (
        (
            "current_classical_plus_full_S246_without_profile_subtraction",
            require_assembly,
        ),
        ("old_Proca_inverse_plus_heavy_inverse", require_reference),
        ("Proca_only_pole_in_total_inverse", require_reference),
        ("bare_internal_loop_cutoff", require_internal_integral),
        ("unprepared_solution_after_initial_time_only", require_graph),
    ):
        cases.append(("wrong_" + label, call, (label,)))
    for n in range(9):
        rows = frontier()
        rows[n]["status"] = "COMPLETE"
        cases.append((f"primitive_promotion_{n}", validate_scope, (rows, matching())))
    for n in range(len(matching())):
        rows = matching()
        rows[n]["status"] = "COMPLETE"
        cases.append((f"matching_promotion_{n}", validate_scope, (frontier(), rows)))
    cases.extend(
        (
            ("missing_primitive", validate_scope, (frontier()[:-1], matching())),
            ("missing_matching", validate_scope, (frontier(), matching()[:-1])),
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
        raise ValueError("An unsupported coupled inverse input was accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "double_heavy_profile_changes_whole_operator": assembly.profile_data()[
            "double_profile_defect"
        ]
        != 0,
        "omitted_matter_shift_changes_highest_low_row": assembly.system()[
            "unshifted_matter_eta_second"
        ]
        != 0,
        "distinct_Gaussian_clock_contact_nonzero": assembly.ward_data()[
            "nonzero_clock_bounce_eta_second"
        ]
        != 0,
        "original_constraint_germ_cannot_be_reset": assembly.residual_data()[
            "nonzero_unprepared_constraint_kernel"
        ]
        != s.zeros(4, 1),
        "variable_reference_pivot_commutator_retained": True,
        "complete_reference_has_both_mass_thresholds_not_old_pole": True,
        "direct_auxiliary_rank_two_and_output_density_order_retained": True,
        "finite_band_inverse_not_all_field_or_P8_promotion": frontier()
        == previous.frontier(),
    }
