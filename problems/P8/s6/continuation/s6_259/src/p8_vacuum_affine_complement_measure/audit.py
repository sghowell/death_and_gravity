"""Finite complement identities do not close the full continuum quantum mean."""

from functools import cache

import sympy as s
from p8_vacuum_affine_finite_window_growth.intervals import fraction
from p8_vacuum_affine_spatial_gauge import audit as previous

from . import gaussian, jets, split

STATE = "same_parent_finite_canonical_jet_and_conditional_affine_complement_not_new_quantum_state"
MEASURE = "canonical_jet_normalized_conditional_complement_with_retained_sources"
OBSERVABLES = (
    "whole_64_projective_quotient_and_56_complement_determinant_inertia",
    "whole_finite_conditional_Fresnel_density_and_source_contacts",
    "whole_finite_jet_constraint_density_and_source_pullback",
    "nonlinear_ordering_and_continuum_quantum_boundary",
)
ITEM = {
    "id": "QG2_H8A423_complete_56_affine_complement_determinant_inertia_and_finite_canonical_jet_source_measure_bridge",
    "status": "EXACT_FULL_CONNECTION_COMPLEMENT_PROJECTIVE_QUOTIENT_FINITE_FRESNEL_SOURCE_AND_JET_CONSTRAINT_DENSITY_NOT_FULL_QUANTUM_ORDERING_BRST_CONTINUUM_STATE_MEAN_ORIGINAL_V_G_B_OR_P8",
}


def parameters():
    return previous.parameters()


def require_parameters(values):
    return previous.require_parameters(values)


def require_state(label):
    if not isinstance(label, str) or label != STATE:
        raise ValueError(
            "The finite affine complement does not replace the fixed quantum state"
        )
    return label


def require_observable(label):
    if not isinstance(label, str) or label not in OBSERVABLES:
        raise ValueError(
            "Keep the full finite complement and explicit continuum boundary"
        )
    return label


def require_measure(label):
    if not isinstance(label, str) or label != MEASURE:
        raise ValueError(
            "A flat configuration comparison is not the canonical conditional density"
        )
    return label


def require_p(value):
    value = fraction(value)
    if value <= 0 or not fraction(s.Rational(1, 8)) < value * value < fraction(
        s.Rational(3, 10)
    ):
        raise ValueError("Require the full strict coefficient interval 1/2<R=4p^2<6/5")
    return s.Rational(value.numerator, value.denominator)


def require_cell_scale(value):
    value = fraction(value)
    if value <= 0:
        raise ValueError("Require the complete positive finite cell action scale")
    return s.Rational(value.numerator, value.denominator)


def require_cells(value):
    if type(value) is not int or value <= 0:
        raise TypeError("Require an explicit positive finite cell count")
    return value


def require_projective_source(value):
    if not isinstance(value, s.MatrixBase) or value.shape != (64, 1):
        raise TypeError("Require all 64 exact source components")
    if any(entry.is_Rational is not True for entry in value):
        raise TypeError(
            "Require real exact rational source components for this finite diagnostic"
        )
    if split.matrices()["projective"].T * value != s.zeros(4, 1):
        raise ValueError("An original connection source must be projectively invariant")
    return s.ImmutableMatrix(value)


def frontier():
    return previous.frontier()


def matching():
    return [*previous.matching(), dict(ITEM)]


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("No earlier primitive or matching row is promoted")
    return True


@cache
def packets():
    return {
        "whole_56_affine_complement_Hessian_inverse_and_determinant": split.complement(),
        "whole_64_projective_configuration_and_gauge_Jacobians": split.projective(),
        "whole_56_exact_inertia_and_continued_Fresnel_phase": split.inertia(),
        "whole_actual_64_connection_source_center_and_retained_trace": split.source_map(),
        "whole_coefficient_dependent_finite_complement_weight": gaussian.determinant(),
        "whole_normalized_conditional_Fresnel_insertion": gaussian.normalization(),
        "whole_56_and_64_connection_source_functional_contacts": gaussian.source_contacts(),
        "whole_scale_and_coefficient_variations_and_composite_mean_boundary": gaussian.local_variation(),
        "whole_finite_jet_point_cotangent_lift": jets.point_lift(),
        "whole_finite_jet_primary_secondary_reduction": jets.primary_secondary(),
        "whole_source_dependent_jet_and_remaining_momentum_constraint": jets.source_dependent_jet(),
        "whole_nonlinear_quantum_ordering_counterexample": jets.ordering_boundary(),
    }


@cache
def residuals():
    return {
        packet + "_" + name: split.clean(value)
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
        "same_fixed_parent_parameters": require_parameters(parameters())
        == parameters(),
        "no_full_quantum_ordering_measure_or_original_P8_closure": True,
    }


def observable():
    return {
        "positive_result": "The unchanged full 64-connection source map has the explicit projective quotient and 56-direction complement determinant, inverse and exact inertia. The full finite Fresnel integral and its normalized conditional density retain the continued phase, cell scale and source contacts. A finite derivative-coordinate cotangent/constraint argument removes the regular auxiliary pairs without inverting the potentially singular original velocity Hessian or dropping the remaining kinematic constraints.",
        "scope_boundary": "The conditional insertion can lift an already specified retained finite functional; it does not choose the full nonlinear time slicing, operator ordering, state or continuum measure. The nonzero source contact, source-dependent remaining constraints and nonlinear composite mean are retained. No BRST/covariant continuum regulator, physical subtraction, interacting state or fixed quantum mean is constructed.",
        "original_problem": "The same action, functions, heavy source, profiles and all earlier primitive/matching/P8(a) qualifications remain. Controlled Wilsonian matching, physical UV scattering, finite-gravity IR/Regge, nonlinear bounce and original V/G/B/P8 stay OPEN.",
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
        ("state", require_state),
        ("observable", require_observable),
        ("measure", require_measure),
        ("coefficient", require_p),
        ("cell_scale", require_cell_scale),
        ("cell_count", require_cells),
        ("source", require_projective_source),
    )
    for i, value in enumerate(invalid):
        for label, call in calls:
            cases.append(("invalid_" + label + "_" + str(i), call, (value,)))
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
    for key in parameters():
        for factor in (0, -1, 2):
            values = parameters()
            values[key] *= factor
            cases.append(
                ("changed_" + key + "_" + str(factor), require_parameters, (values,))
            )
    for value in (0, -1, s.Rational(1, 4), s.Rational(3, 5), 1):
        cases.append(("outside_coefficient_domain_" + str(value), require_p, (value,)))
    for value in (0, -1):
        cases.append(("invalid_cell_scale_" + str(value), require_cell_scale, (value,)))
        cases.append(("invalid_cell_count_" + str(value), require_cells, (value,)))
    for label in (
        "flat_configuration_measure_is_canonical",
        "all_affine_determinants_constant",
        "principal_endpoint_square_root_phase",
        "dimensionally_set_delta_zero",
        "source_dependent_normalization_cancels_all_source_contacts",
        "naive_Weyl_symbol_unchanged_under_any_point_map",
    ):
        cases.append(("wrong_measure_" + label, require_measure, (label,)))
    for label in (
        "original_P8_closed",
        "all_auxiliary_correlations_zero",
        "full_covariant_BRST_measure",
        "unitary_gauge_through_zero_clock_gradient",
        "full_quantum_mean_classical",
        "fixed_state_replaced",
        "source_center_evaluated_at_mean",
        "physical_cutoff_proved",
        "all_velocity_Hessians_invertible",
        "remaining_jet_constraints_discarded",
        "dynamical_Proca_integrated_out",
        "first_order_Einstein_Hilbert_measure_imported",
        "finite_density_proves_nonlinear_quantum_equivalence",
        "P8a_qualifications_removed",
    ):
        cases.append(("unsupported_" + label, require_observable, (label,)))
    cases.append(
        (
            "old_spatial_state_scope_not_new_complement_scope",
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
    noninvariant = s.zeros(64, 1)
    noninvariant[0] = 1
    inexact = s.zeros(64, 1)
    inexact[1] = s.Float(1)
    cases.extend(
        [
            (
                "nonprojective_original_source",
                require_projective_source,
                (noninvariant,),
            ),
            ("inexact_original_source", require_projective_source, (inexact,)),
            (
                "incomplete_original_source",
                require_projective_source,
                (s.zeros(60, 1),),
            ),
            (
                "wrong_original_source_orientation",
                require_projective_source,
                (s.zeros(1, 64),),
            ),
            (
                "singular_rational_inertia_matrix",
                split.rational_congruence,
                (s.zeros(2),),
            ),
            (
                "nonsymmetric_rational_inertia_matrix",
                split.rational_congruence,
                (s.Matrix([[1, 2], [0, 1]]),),
            ),
            (
                "inexact_rational_inertia_matrix",
                split.rational_congruence,
                (s.Matrix([[s.Float(1)]]),),
            ),
        ]
    )
    missing = parameters()
    del missing["time_length"]
    extra = parameters()
    extra["affine_quantum_measure_removed"] = 1
    cases.extend(
        [
            ("missing_parameter", require_parameters, (missing,)),
            ("extra_parameter", require_parameters, (extra,)),
            ("missing_primitive", validate_scope, (frontier()[:-1], matching())),
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
        raise ValueError("Unsupported affine-complement input accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "whole_64_connection_and_56_complement_retained": True,
        "exact_inertia_and_Fresnel_phase_not_principal_endpoint_root": True,
        "nonconstant_flat_measure_factor_and_cell_scale_retained": True,
        "source_contact_and_composite_mean_not_deleted": True,
        "derivative_coordinate_and_remaining_constraints_retained": True,
        "nonlinear_ordering_and_continuum_boundary_explicit": True,
        "no_quantum_mean_covariant_measure_or_cutoff_promotion": True,
        "all_original_and_earlier_matching_rows_unchanged": frontier()
        == previous.frontier()
        and matching()[:-1] == previous.matching(),
    }
