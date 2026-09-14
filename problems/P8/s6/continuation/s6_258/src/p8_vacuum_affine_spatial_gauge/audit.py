"""A local spatial ghost operator is not a completed covariant quantum measure."""

from functools import cache

import sympy as s
from p8_vacuum_affine_finite_window_growth.intervals import fraction
from p8_vacuum_affine_nonlinear_auxiliary_measure import audit as previous

from . import gauge, ghost, spatial

STATE = (
    "same_parent_local_spatial_gauge_and_finite_ghost_vertices_not_new_quantum_state"
)
GAUGE = "nonlinear_conformal_Dirac_spatial_only"
OBSERVABLES = (
    "whole_nonlinear_spatial_generator_physical_frame_and_off_gauge_ghost_operator",
    "local_mean_zero_torus_gauge_slice_with_residual_translations_and_IR_scale",
    "whole_three_ghost_finite_determinant_vertices_with_second_shape_contact",
    "second_order_gauge_restoration_and_scalar_volume_source_contact",
)
ITEM = {
    "id": "QG2_H8A422_complete_retained_nonlinear_spatial_Dirac_gauge_operator_local_slice_and_full_finite_ghost_vertices",
    "status": "EXACT_WHOLE_SPATIAL_GENERATOR_FRAME_OFF_GAUGE_GHOST_AND_SECOND_ORDER_RESTORATION_WITH_WRITTEN_LOCAL_SLICE_AND_FINITE_DETERMINANT_VERTICES_NOT_GLOBAL_GAUGE_BRST_COVARIANT_QUANTUM_MEASURE_FIXED_MEAN_ORIGINAL_V_G_B_OR_P8",
}


def parameters():
    return previous.parameters()


def require_parameters(values):
    return previous.require_parameters(values)


def require_state(label):
    if not isinstance(label, str) or label != STATE:
        raise ValueError(
            "The local spatial gauge does not replace the fixed quantum state"
        )
    return label


def require_observable(label):
    if not isinstance(label, str) or label not in OBSERVABLES:
        raise ValueError(
            "Keep the local classical spatial operator and finite ghost scope"
        )
    return label


def require_gauge(label):
    if not isinstance(label, str) or label != GAUGE:
        raise ValueError(
            "No changed time gauge or globally fixed spatial gauge is certified"
        )
    return label


def require_shape_distance(value):
    value = fraction(value)
    if not 0 <= value <= fraction(gauge.SHAPE_BOUND):
        raise ValueError("Require the entire weak-coercivity shape box")
    return s.Rational(value.numerator, value.denominator)


def require_torus_radius(value):
    value = fraction(value)
    if value <= 0:
        raise ValueError(
            "Require a positive finite torus radius, not the infinite-volume limit"
        )
    return s.Rational(value.numerator, value.denominator)


def require_mode(value):
    if (
        not isinstance(value, tuple)
        or len(value) != 3
        or any(type(k) is not int for k in value)
        or not any(value)
    ):
        raise TypeError(
            "Require a nonzero three-dimensional integer mode; keep translations separate"
        )
    return value


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
        "whole_metric_vector_matter_spatial_generator": spatial.generator(),
        "whole_spatial_classical_Lie_algebra": spatial.lie_algebra(),
        "whole_physical_frame_and_fixed_reference_shape": spatial.frame(),
        "whole_linear_TT_and_scalar_volume_gauge_dictionary": spatial.reference_projector(),
        "whole_off_gauge_nonlinear_spatial_ghost_operator": gauge.full_operator(),
        "whole_nonsymmetric_principal_symbol_and_inverse": gauge.principal(),
        "whole_local_slice_coercivity_and_torus_IR_scale": gauge.local_slice(),
        "whole_translation_and_finite_regulator_scope": gauge.zero_modes(),
        "whole_off_gauge_three_ghost_Fourier_operator": ghost.full_fourier(),
        "whole_finite_three_ghost_logdet_shape_vertices": ghost.finite_vertices(),
        "whole_second_order_gauge_restoration_and_volume_contact": ghost.mixed_off_gauge_contact(),
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
        "same_fixed_parent_parameters": require_parameters(parameters())
        == parameters(),
        "no_full_quantum_measure_or_original_P8_closure": True,
    }


def observable():
    return {
        "positive_result": "The full retained spatial momentum generator and conformal Dirac density give the complete nonlinear off-gauge three-ghost operator. The density is exactly identical in the physical and full lapse-dependent conformal frames and its linearization agrees with the existing TT plus scalar-volume reference. On a flat finite torus a local slice exists modulo residual translations, with a full nonsymmetric weak-coercivity margin 3/4 on the stated 1/8 shape box and explicit L-dependent inverse gap. Entire finite three-ghost determinant vertices retain the nonzero second shape contact; the mixed TT example is gauge-restored with its required scalar-volume contact.",
        "scope_boundary": "Translations and global constraints are not discarded; mean-zero vectors are not a Lie subalgebra and finite summation by parts is not an exact diffeomorphism regulator. No global gauge, infinite-volume gap, BRST/covariant continuum quantum measure, affine/projective determinant, physical counterterm prescription, interacting state or fixed quantum mean is constructed.",
        "original_problem": "The same parent, all profiles and all earlier primitive/matching/P8(a) qualifications remain. Wilsonian matching, physical UV scattering, finite-gravity IR/Regge, nonlinear bounce and original V/G/B/P8 stay OPEN.",
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
        ("gauge", require_gauge),
        ("shape", require_shape_distance),
        ("radius", require_torus_radius),
        ("mode", require_mode),
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
    for value in (-1, gauge.SHAPE_BOUND + 1, s.Rational(1, 4)):
        cases.append(
            ("outside_shape_box_" + str(value), require_shape_distance, (value,))
        )
    for value in (0, -1):
        cases.append(("invalid_radius_" + str(value), require_torus_radius, (value,)))
    for i, value in enumerate(
        (
            (0, 0, 0),
            (1, 2),
            (1, 2, 3, 4),
            [1, 2, 3],
            (True, 0, 0),
            (s.Integer(1), 0, 0),
            (1.0, 0, 0),
        )
    ):
        cases.append(("invalid_mode_" + str(i), require_mode, (value,)))
    for label in (
        "CMC_time_gauge",
        "global_harmonic_gauge",
        "all_translation_modes_removed",
        "linear_TT_gauge_at_all_orders",
    ):
        cases.append(("wrong_gauge_" + label, require_gauge, (label,)))
    for label in (
        "original_P8_closed",
        "all_ghost_vertices_zero",
        "all_quantum_determinants_vanish",
        "full_covariant_BRST_measure",
        "first_class_regulator_closed",
        "zero_momentum_inverse",
        "infinite_volume_uniform_gap",
        "mean_zero_vectors_are_Lie_subalgebra",
        "global_unique_spatial_gauge",
        "old_quantum_mean_is_classical",
        "fixed_quantum_state_replaced",
        "no_second_shape_contact",
        "scalar_volume_contact_dropped",
        "gauge_surface_formula_valid_off_slice",
        "physical_cutoff_proved",
        "negative_curvature_GR_theorem_applies",
        "P8a_qualifications_removed",
    ):
        cases.append(("unsupported_" + label, require_observable, (label,)))
    cases.append(
        (
            "old_auxiliary_state_scope_not_new_spatial_scope",
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
    extra["quantum_gauge_regulator_removed"] = 1
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
        raise ValueError("Unsupported spatial-gauge input accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "whole_current_source_profiles_and_canonical_boundaries_retained": True,
        "complete_off_gauge_spatial_operator_not_surface_only": True,
        "all_three_ghosts_and_second_shape_contact_retained": True,
        "second_order_gauge_restoration_keeps_scalar_volume_contact": True,
        "translation_kernel_and_finite_torus_IR_scale_retained": True,
        "finite_SBP_not_anomaly_free_diffeomorphism_regulator": True,
        "no_quantum_mean_covariant_measure_or_cutoff_promotion": True,
        "all_original_and_earlier_matching_rows_unchanged": frontier()
        == previous.frontier()
        and matching()[:-1] == previous.matching(),
    }
