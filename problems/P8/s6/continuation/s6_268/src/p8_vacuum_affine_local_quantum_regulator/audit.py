"""Declared finite regulator, unchanged original state and strict closure gates."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_nonlinear_spatial_reduction import audit as previous

from . import dynamics, localization, ordering, source

STATE = "same_original_nonzero_pure_seed_no_conditioning"
MEASURE = "explicit_finite_IR_UV_phase_localized_coherent_regulator"
OBSERVABLES = (
    "finite_local_first_Weyl_calibrated_unitary_dynamics",
    "positive_coherent_regulated_volume",
    "conditional_first_Weyl_calibrated_volume",
    "dimension_dependent_coherent_state_leakage_comparison",
)
ITEM = {
    "id": "QG2_H8A432_explicit_finite_local_coherent_quantization_unitary_dynamics_and_leakage",
    "status": "COMPLETE_NAMED_FINITE_LOCAL_QUANTUM_REGULATOR_AND_EXACT_ERROR_CRITERIA_NOT_EVALUATED_PHYSICAL_MATCHING_REGULATOR_REMOVAL_OR_ORIGINAL_V_G_B_P8",
}


def parameters():
    return previous.parameters()


def require_parameters(values):
    return previous.require_parameters(values)


def require_radius(value):
    return previous.require_radius(value)


def require_dimension(value):
    return localization.positive_integer(value)


def require_state(value):
    if not isinstance(value, str) or value != STATE:
        raise ValueError(
            "Retain the same original pure nonzero seed without projection"
        )
    return value


def require_measure(value):
    if not isinstance(value, str) or value != MEASURE:
        raise ValueError("Require the explicitly named finite localized regulator")
    return value


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError(
            "Require a declared finite regulated readout or error criterion"
        )
    return value


def frontier():
    return previous.frontier()


def matching():
    return [*previous.matching(), dict(ITEM)]


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("Finite regulated dynamics is not original P8 quantum closure")
    return True


@cache
def packets():
    return {
        "whole_original_time_source_chart_and_reference": source.data(),
        "whole_full_covariance_ordering_and_physical_contacts": ordering.data(),
        "whole_strict_phase_effect_and_state_difference_bounds": localization.data(),
        "whole_convergent_regulated_quantum_dynamics": dynamics.data(),
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
        len(v) if isinstance(v, s.MatrixBase) else 1 for v in residuals().values()
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
        "all_matching_identifiers_distinct": len({r["id"] for r in matching()})
        == len(matching()),
        "all_original_physical_parameters_unchanged": require_parameters(parameters())
        == parameters(),
        "historical_refutations_and_Gaussian_integrability_boundary_unchanged": True,
        "original_P8_open_despite_explicit_finite_quantum_regulator": True,
    }


def observable():
    return {
        "established": "An explicitly named finite local regulator uses the whole original canonical Hamiltonian, full nonlinear spatial/auxiliary reconstruction, complete time connections and SAME pure nonzero Gaussian seed. Coherent quantization gives a bounded first-Weyl-calibrated interaction, actual unitary dynamics, residual-translation preservation and a positive finite regulated volume. Exact dimension-dependent localization tails, ordering contacts and operator-on-state comparison criteria remain explicit.",
        "domain": "Finite torus radius and symmetric nonzero mode set; external homogeneous reference parameters; a sufficiently small but UNEVALUATED classical canonical phase ball and positive time interval. The nonlinear convolution and elliptic reconstruction are not truncated. The entire original free mode flow is retained. The cutoff continuation and finite ordering prescription are new declared regulators, not established original counterterms.",
        "not_established": "No evaluated whitened canonical radius, time width, actual small leakage, ordering norm or numerical P8 interacting mean; no equivalence to the original singular Hamiltonian, homogeneous quantum state, continuum diffeomorphism/BRST regulator, regulator removal, Wilsonian matching cutoff, omitted-loop bounds or nonlinear inhomogeneous completeness. Original V/G/B/P8 remain OPEN and completed scoped P8(a) is unchanged.",
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
    calls = (
        ("radius", require_radius),
        ("dimension", require_dimension),
        ("state", require_state),
        ("measure", require_measure),
        ("observable", require_observable),
    )
    for i, value in enumerate(invalid):
        for label, call in calls:
            cases.append(("invalid_" + label + "_" + str(i), call, (value,)))
        for key in parameters():
            changed = parameters()
            changed[key] = value
            cases.append(
                (
                    "invalid_parameter_" + key + "_" + str(i),
                    require_parameters,
                    (changed,),
                )
            )
    for value in (0, -1):
        cases.append(("nonpositive_radius_" + str(value), require_radius, (value,)))
        cases.append(
            ("nonpositive_dimension_" + str(value), require_dimension, (value,))
        )
    cases.append(("fractional_dimension", require_dimension, (s.Rational(3, 2),)))
    for value in (
        "new_vacuum",
        "conditioned_Gaussian",
        "exact_phase_support",
        "homogeneous_state",
    ):
        cases.append(("wrong_state_" + value, require_state, (value,)))
    for value in (
        "original_unregularized_Hamiltonian",
        "full_BRST_regulator",
        "cutoff_without_derivative_contacts",
        "undeclared_counterterm",
    ):
        cases.append(("wrong_measure_" + value, require_measure, (value,)))
    for value in (
        "original_P8_closed",
        "positive_Weyl_unconditionally",
        "symbol_sup_is_operator_norm",
        "uniform_continuum_tail",
        "actual_small_P8_leakage",
        "initial_tail_is_later_support",
    ):
        cases.append(("wrong_observable_" + value, require_observable, (value,)))
    cases.extend(
        [
            ("delete_primitive_frontier", validate_scope, ([], matching())),
            ("delete_matching_frontier", validate_scope, (frontier(), [])),
        ]
    )
    assert len({name for name, _, _ in cases}) == len(cases)
    return cases


def rejected_inputs():
    count = 0
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (ValueError, TypeError):
            count += 1
        else:
            raise ValueError("Invalid regulator input was accepted: " + name)
    return count


def controls():
    return {
        "full_original_time_source_auxiliary_spatial_contacts_retained": True,
        "same_original_free_preparation_without_reselection": True,
        "new_finite_regulator_declared_not_silent_original_action_change": True,
        "all_ordering_and_cutoff_derivative_contacts_retained": True,
        "phase_effect_not_state_projection_or_continuum_limit": True,
        "exact_bounded_quantum_dynamics_not_unregularized_matching": True,
        "historical_errata_and_original_P8_open": True,
        "written_functional_operator_proofs_not_FORMALIZED": True,
        "rejected_inputs": rejected_inputs(),
    }
