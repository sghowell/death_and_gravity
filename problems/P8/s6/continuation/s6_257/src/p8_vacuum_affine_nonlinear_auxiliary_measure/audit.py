"""A full classical auxiliary branch is not a completed quantum measure."""

from functools import cache

import sympy as s
from p8_vacuum_affine_finite_window_growth import audit as previous
from p8_vacuum_affine_finite_window_growth.intervals import fraction

from . import canonical, cotangent, measure

STATE = (
    "same_parent_local_unforced_classical_auxiliary_branch_not_fixed_quantum_reference"
)
OBSERVABLES = (
    "whole_retained_nonlinear_unitary_canonical_Hamiltonian_and_local_auxiliary_branch",
    "whole_extended_cotangent_and_primitive_boundary_canonical_measure",
    "single_regular_finite_regulator_second_class_density_with_gauge_factors_retained",
    "full_nonlinear_source_pullback_and_classical_quantum_centering_boundary",
)
ITEM = {
    "id": "QG2_H8A421_complete_retained_nonlinear_canonical_auxiliary_branch_cotangent_boundary_second_class_measure",
    "status": "EXACT_WHOLE_CURRENT_NONLINEAR_RETAINED_CANONICAL_AUXILIARY_AND_BOUNDARY_IDENTITIES_WITH_WRITTEN_LOCAL_BRANCH_AND_FINITE_REGULATOR_SECOND_CLASS_REDUCTION_NOT_FULL_GAUGE_AFFINE_COVARIANT_QUANTUM_MEASURE_FIXED_MEAN_ORIGINAL_V_G_B_OR_P8",
}


def parameters():
    return previous.parameters()


def require_parameters(values):
    return previous.require_parameters(values)


def require_finite_cells(value):
    if type(value) is not int or value < 1:
        raise TypeError("Require an explicit positive finite integer cell count")
    return value


def require_pivot(value):
    exact = fraction(value)
    if exact == 0:
        raise ValueError("A singular auxiliary Jacobian is not a regular root branch")
    return s.Rational(exact.numerator, exact.denominator)


def require_branch(label):
    if not isinstance(label, str) or label != "one_local_regular_auxiliary_root":
        raise ValueError("Do not identify separate roots or assert global uniqueness")
    return label


def require_state(label):
    if not isinstance(label, str) or label != STATE:
        raise ValueError("The classical branch is not a replacement quantum reference")
    return label


def require_observable(label):
    if not isinstance(label, str) or label not in OBSERVABLES:
        raise ValueError("Keep the retained classical and finite second-class scope")
    return label


def frontier():
    return previous.frontier()


def matching():
    return [*previous.matching(), dict(ITEM)]


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("No primitive or earlier matching obligation is promoted")
    return True


@cache
def packets():
    return {
        "whole_current_nonlinear_trace_temporal_and_lapse_Hamiltonian": canonical.full_trace(),
        "whole_homogeneous_full_potential_lapse_Schur_bridge": canonical.homogeneous_bridge(),
        "whole_actual_current_functions_profiles_and_regular_classical_datum": canonical.current_data(),
        "whole_extended_all_field_cotangent_map": cotangent.extended_map(),
        "whole_primitive_boundary_and_combined_canonical_map": cotangent.primitive_boundary_map(),
        "whole_covariant_and_unitary_metric_density_bridge": cotangent.metric_density_bridge(),
        "whole_spatial_metric_vector_matter_and_Gauss_Legendre_blocks": cotangent.spatial_kinetic_blocks(),
        "whole_finite_regulator_second_class_measure": measure.block_measure(),
        "whole_nonlinear_finite_constraint_diagnostic": measure.nonlinear_example(),
        "whole_source_pullback_and_fixed_reference_centering_boundary": measure.source_centering(),
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
        "unchanged_parent_parameters": require_parameters(parameters()) == parameters(),
        "no_original_P8_or_full_quantum_measure_closure": True,
    }


def observable():
    return {
        "positive_result": "The full retained spatially inhomogeneous unitary Hamiltonian admits a local algebraic lapse/normal-vector auxiliary branch about the actual full-current S256 classical datum. Its exact temporal pivot and full lapse Schur bridge are nonzero. All spatial kinetic, gradient, curvature, heavy-source, primitive and Gauss terms remain. The entire extended point cotangent lift and primitive boundary shift preserve canonical Liouville volume. At any fixed finite regulator on one regular root branch the second-class density cancels the full auxiliary delta Jacobian even with nonzero secondary brackets.",
        "source_centering": "Original generating-function sources pull back through the complete nonlinear map, including auxiliary Hessian contacts. The exact old classical reference constraint profile is retained. This is not the fixed interacting quantum mean and does not automatically reconstruct the specified off-shell S252 Gaussian.",
        "scope_boundary": "Remaining spatial gauge factors, affine-complement and projective-gauge quantization, covariant/BRST continuum regulation, ordering and local counterterms, physical quantum subtraction, interacting state and mean, cutoff and matching errors are not supplied. No arbitrary inhomogeneous nonlinear Cauchy theorem or global auxiliary uniqueness follows.",
        "original_problem": "All original V/G/B and P8(a) qualifications remain unchanged; original P8 is OPEN.",
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
        ("cells", require_finite_cells),
        ("pivot", require_pivot),
        ("branch", require_branch),
    )
    for i, value in enumerate(invalid):
        for label, call in calls:
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
    for value in (0, -1, s.Rational(1, 2), s.Integer(2)):
        cases.append(
            ("invalid_finite_cells_" + str(value), require_finite_cells, (value,))
        )
    cases.append(("singular_auxiliary_pivot", require_pivot, (0,)))
    for label in (
        "all_roots_identified",
        "global_unique_auxiliary_root",
        "singular_branch",
        "original_reference_is_unforced",
    ):
        cases.append(("invalid_branch_" + label, require_branch, (label,)))
    for label in (
        "original_P8_closed",
        "all_quantum_determinants_vanish",
        "full_covariant_BRST_measure",
        "spatial_gauge_determinant_deleted",
        "affine_complement_quantized",
        "fixed_quantum_mean_replaced",
        "old_profile_constraint_zero",
        "old_profile_constraint_nonzero",
        "all_sources_linear_after_map",
        "primitive_boundary_discarded",
        "instantaneous_vacuum_reset",
        "nonlinear_configuration_integral_Gaussian",
        "nonlinear_inhomogeneous_Cauchy_theorem",
        "all_auxiliary_roots_same",
        "cutoff_and_matching_closed",
        "P8a_qualifications_removed",
        "classical_comparison_equals_old_Gaussian",
    ):
        cases.append(("unsupported_" + label, require_observable, (label,)))
    cases.append(
        (
            "old_finite_band_scope_not_new_measure_scope",
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
    missing = parameters()
    del missing["time_length"]
    extra = parameters()
    extra["quantum_regulator_removed"] = 1
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
        raise ValueError("Unsupported nonlinear auxiliary input accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "whole_current_functions_profiles_constants_and_heavy_source_retained": True,
        "full_spatial_canonical_Hamiltonian_and_Gauss_density_retained": True,
        "complete_cotangent_and_primitive_boundary_maps_retained": True,
        "secondary_constraint_bracket_not_assumed_zero": True,
        "finite_single_regular_branch_not_continuum_quantum_gauge_measure": True,
        "full_source_contacts_and_fixed_reference_centering_boundary_retained": True,
        "no_classical_auxiliary_result_promoted_to_quantum_mean_or_cutoff": True,
        "all_original_and_earlier_matching_rows_unchanged": frontier()
        == previous.frontier()
        and matching()[:-1] == previous.matching(),
    }
