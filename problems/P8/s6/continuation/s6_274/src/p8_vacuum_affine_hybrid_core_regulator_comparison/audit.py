"""Finite hybrid comparisons, exact guards and unchanged original frontiers."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_selfconsistent_finite_feedback import audit as previous

from . import domain, dynamics, operators, quadratic, source

STATE, MEASURE, MODEL = previous.STATE, previous.MEASURE, previous.MODEL
OBSERVABLES = (
    "finite_hybrid_positive_coherent_core_probability",
    "coupled_finite_hybrid_cutoff_comparison",
    "coupled_finite_hybrid_ordering_comparison",
    "finite_hybrid_physical_volume_regulator_comparison",
)
ITEM = {
    "id": "QG2_H8A438_same_full_source_finite_hybrid_core_probability_and_coupled_regulator_comparisons",
    "status": "COMPLETE_SAME_FINITE_HYBRID_CORE_POVM_AND_TWO_CUTOFF_ORDERING_COMPARISONS_NOT_EXACT_SUPPORT_REGULATOR_REMOVAL_OR_ORIGINAL_V_G_B_P8",
}


def parameters():
    return previous.parameters()


def require_parameters(value):
    return previous.require_parameters(value)


def require_radius(value):
    return previous.require_radius(value)


def require_state(value):
    return previous.require_state(value)


def require_measure(value):
    return previous.require_measure(value)


def require_cutoff(value):
    return previous.require_cutoff(value)


def require_time(value):
    return previous.require_time(value)


def require_model(value):
    return previous.require_model(value)


def require_ordering(value):
    return previous.require_ordering(value)


def require_homogeneous_deviation(value):
    return previous.require_homogeneous_deviation(value)


def require_amplitude(value):
    return operators.require_amplitude(value)


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError("Require a proved same-model finite hybrid comparison")
    return value


def frontier():
    return previous.frontier()


def matching():
    return [*previous.matching(), dict(ITEM)]


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("Finite regulator comparisons do not close original P8")
    return True


@cache
def packets():
    return {
        "whole_original_source_and_refined_full_constraint": source.data(),
        "whole_literal_constrained_quadratic_energy": quadratic.data(),
        "whole_large_analysis_domain_and_complete_remainder": domain.data(),
        "whole_exact_operator_identity_and_all_phase_jets": operators.data(),
        "whole_coupled_core_probability_and_regulator_comparisons": dynamics.data(),
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
        "all_nine_original_primitive_rows_retained": len(frontier()) == 9,
        "all_matching_ids_distinct": len({row["id"] for row in matching()})
        == len(matching()),
        "same_original_parameters": require_parameters(parameters()) == parameters(),
        "same_finite_hybrid_model": require_model(MODEL) == MODEL,
        "same_original_physical_core": require_radius(source.RADIUS) == source.RADIUS,
        "no_frozen_ancestor_or_scoped_P8a_upgrade": True,
        "original_P8_remains_open": True,
    }


def observable():
    return {
        "established": "For the SAME four S273 finite hybrid solutions on |u|<=1e-180, the transported positive coherent outside-core probability is<1e-6. Between the original two cutoffs at fixed ordering, five-coordinate, phase-factored state and physical-volume differences are<1e-580,1e-12,1e-514. Between calibrated-coherent and Weyl at fixed cutoff they are<1e-630,1e-65,1e-566. All background feedback, scalar-center forces and complete operator remainders are retained.",
        "domain": "Unchanged L1,P1e64,kappa1e800,zeta1e-6,physical coreR1e20,original48-pair pure seed and fixed reference. The full nonlinear source and auxiliary/spatial fields are unchanged. Rstar1e150 is a complex analysis neighborhood only; the improved full centered Hamiltonian bound1e172 is not a replacement Hamiltonian.",
        "not_established": "No exact support, sharp joint phase projection, comparison of un-factored state vectors with different global phases, homogeneous quantization, unlocalized Hamiltonian, uniform mode/volume/cutoff limit, unique or strict volume minimum, physical matching, omitted-loop/Regge/UV control or nonlinear global completion. Original V/G/B/P8 remain OPEN; scoped P8(a) and historical refutations remain unchanged.",
    }


def bad_cases():
    cases = [
        ("inherited_" + name, call, args) for name, call, args in previous.bad_cases()
    ]
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
        for name, call in (
            ("amplitude", require_amplitude),
            ("observable", require_observable),
            ("analysis_radius", domain.bounds),
            ("radial_partition", operators.radial_partition),
            ("phase_order", operators.coefficient),
        ):
            cases.append(("invalid_comparison_" + name + "_" + str(i), call, (value,)))
    for i, value in enumerate((-1, 207, 208, s.Rational(1, 2))):
        for name, call in (
            ("coefficient", operators.coefficient),
            ("radial", operators.radial_partition),
        ):
            cases.append(("unproved_extended_" + name + "_" + str(i), call, (value,)))
        cases.append(
            (
                "unproved_extended_jet_" + str(i),
                operators.jet,
                (value, source.H_AMPLITUDE),
            )
        )
    for i, value in enumerate((-1, 3, True, 1.0, s.Rational(1, 2))):
        cases.append(
            (
                "unproved_homogeneous_order_" + str(i),
                operators.bounds,
                (source.H_AMPLITUDE, value),
            )
        )
    for i, value in enumerate((-1, 8, True, 1.0, s.Rational(1, 2))):
        cases.append(
            (
                "unproved_heat_order_" + str(i),
                operators.heat_term,
                (value, source.H_AMPLITUDE),
            )
        )
    for i, value in enumerate(
        (0, -1, 1, source.previous.previous.TIME, 10**280, 10**1000)
    ):
        cases.append(("wrong_new_amplitude_" + str(i), operators.bounds, (value,)))
    for i, value in enumerate(
        (0, -1, 1, source.RADIUS * 2, source.ANALYSIS_RADIUS * 2)
    ):
        cases.append(("unproved_analysis_radius_" + str(i), domain.bounds, (value,)))
    for i, value in enumerate(
        (
            "original_P8_closed",
            "exact_zero_core_leakage",
            "sharp_joint_phase_support",
            "homogeneous_quantization",
            "regulator_removed",
            "unlocalized_Hamiltonian",
            "identify_full_Hessian_with_free_reference",
            "delete_antiheat_remainder",
            "core_support_of_heat_remainder",
            "compare_unfactored_global_state_phases",
            "drop_scalar_center_force",
            "larger_physical_cutoff",
            "new_isotropic_seed",
            "externally_identical_background_histories",
            "unique_or_strict_minimum",
        )
    ):
        cases.append(("wrong_comparison_scope_" + str(i), require_observable, (value,)))
    cases.extend(
        [
            ("delete_original_frontier", validate_scope, ([], matching())),
            ("delete_new_matching_frontier", validate_scope, (frontier(), [])),
        ]
    )
    if len({name for name, _, _ in cases}) != len(cases):
        raise ValueError("Duplicate rejected finite comparison input")
    return cases


def rejected_inputs():
    count = 0
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (ValueError, TypeError):
            count += 1
        else:
            raise ValueError("Unsupported finite comparison input accepted: " + name)
    return count


def controls():
    return {
        "same_original_source_state_and_finite_hybrid_solutions": True,
        "full_trace_shear_and_complete_C1_not_free_Hessian": True,
        "all_full_tensor_curvature_and_primitive_contacts": True,
        "longitudinal_Gauss_and_transverse_curl_both_retained": True,
        "larger_analysis_domain_not_larger_physical_cutoff": True,
        "all206_phase_jets_without_more_real_time_regularity": True,
        "exact_positive_heat_identity_with_full_remainder": True,
        "positive_POVM_probability_not_exact_phase_support": True,
        "both_live_background_and_quantum_stability_rows": True,
        "own_scalar_phases_factored_for_vector_norm_only": True,
        "physical_volume_parameter_state_and_operator_changes": True,
        "historical_errata_and_all_original_frontiers_retained": True,
        "rejected_inputs": rejected_inputs(),
    }
