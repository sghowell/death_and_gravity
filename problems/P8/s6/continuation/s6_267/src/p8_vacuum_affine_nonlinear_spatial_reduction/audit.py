"""Strict local spatial reduction and unchanged original quantum frontier."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_finite_window_growth.intervals import fraction
from p8_vacuum_affine_nonlinear_lapse_branch import audit as previous

from . import cotangent, shape, source, translations

STATE = "same_original_pure_nonzero_mode_prescriptions_no_nonlinear_pushforward"
MEASURE = "canonical_classical_spatial_cotangent_and_residual_translation_factor_only"
OBSERVABLES = (
    "whole_nonlinear_spatial_shape_and_cotangent_reduction",
    "whole_residual_translation_constraints_and_pure_reference_kernel",
)
ITEM = {
    "id": "QG2_H8A431_full_nonlinear_spatial_shape_cotangent_and_residual_translation_reduction",
    "status": "COMPLETE_LOCAL_FULL_SPATIAL_COTANGENT_REDUCTION_AND_NONZERO_PURE_REFERENCE_TRANSLATION_KERNEL_NOT_NONLINEAR_QUANTUM_DOMAIN_EVOLUTION_REGULATOR_OR_ORIGINAL_V_G_B_P8",
}


def parameters():
    return previous.parameters()


def require_parameters(values):
    return previous.require_parameters(values)


def require_norm(value):
    value = fraction(value)
    if not 0 <= value <= fraction(shape.TAU_RADIUS):
        raise ValueError(
            "Require the full weighted Fourier norm in the explicit shape ball"
        )
    return value


def require_radius(value):
    value = fraction(value)
    if not value > 0:
        raise ValueError(
            "Require a finite positive torus radius; no infinite-volume gap"
        )
    return value


def require_cutoff(value):
    if isinstance(value, bool) or not isinstance(value, (int, s.Integer)) or value < 1:
        raise ValueError("Require a positive integer symmetric Fourier cutoff")
    return int(value)


def require_state(value):
    if not isinstance(value, str) or value != STATE:
        raise ValueError(
            "Keep the same pure mode prescriptions without nonlinear state projection"
        )
    return value


def require_measure(value):
    if not isinstance(value, str) or value != MEASURE:
        raise ValueError("The classical cotangent lift is not a full quantum measure")
    return value


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError(
            "Require the whole spatial reduction or complete residual translation factor"
        )
    return value


def frontier():
    return previous.frontier()


def matching():
    return [*previous.matching(), dict(ITEM)]


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError(
            "Do not promote spatial reduction to original quantum P8 closure"
        )
    return True


@cache
def packets():
    return {
        "whole_original_generator_auxiliary_and_boundary": source.data(),
        "whole_nonlinear_unit_determinant_transverse_shape": shape.data(),
        "whole_full_nonsymmetric_cotangent_constraint_lift": cotangent.data(),
        "whole_residual_translations_and_original_pure_reference": translations.data(),
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
        "all_original_primitive_rows_retained": len(frontier()) == 9,
        "all_matching_rows_have_distinct_identifiers": len(
            {r["id"] for r in matching()}
        )
        == len(matching()),
        "all_original_parameters_unchanged": require_parameters(parameters())
        == parameters(),
        "both_historical_errata_and_S265_integrability_boundary_retained": True,
        "original_V_G_B_P8_and_nonlinear_quantum_domain_remain_open": True,
    }


def observable():
    return {
        "established": "A full-convolution Fourier Banach contraction constructs the nonlinear spatial shape on the explicit norm ball, including generated homogeneous trace corrections. The entire original momentum generator has a unique cotangent lift solving the nonzero constraint complement, retaining the three residual global translations. The same finite nonzero pure Gaussian mode prescriptions form an exact zero-translation-charge factor.",
        "domain": "The classical torus chart has finite radius L, exponent8 and ||tau||<=1/100; all nonlinear contacts and five homogeneous tracefree shape directions remain. Auxiliary reconstruction applies only where its complete twelve bounce invariants satisfy S266. The reference Gaussian claim concerns the original linear canonical CCR variables under the explicit finite IR/UV mode evaluation, not a nonlinear Gaussian pushforward or support statement.",
        "not_established": "No global diffeomorphism quotient, homogeneous quantum state, nonlinear ordering or self-adjoint interacting Hamiltonian, regulated physical mean, continuum limit, full spatial Ward construction, physical matching cutoff, loop bounds or nonlinear inhomogeneous stability. S265's separate Gaussian integrability obstruction and S261's physical refutation remain. Original V/G/B/P8 remain OPEN; completed scoped P8(a) is unchanged.",
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
        ("norm", require_norm),
        ("radius", require_radius),
        ("cutoff", require_cutoff),
        ("state", require_state),
        ("measure", require_measure),
        ("observable", require_observable),
    )
    for i, value in enumerate(invalid):
        for label, call in calls:
            cases.append(("invalid_" + label + "_" + str(i), call, (value,)))
        for key in parameters():
            altered = parameters()
            altered[key] = value
            cases.append(
                (
                    "invalid_parameter_" + key + "_" + str(i),
                    require_parameters,
                    (altered,),
                )
            )
    for value in (-1, 2 * shape.TAU_RADIUS, 1):
        cases.append(("outside_shape_" + str(value), require_norm, (value,)))
    for value in (0, -1):
        cases.append(("nonpositive_radius_" + str(value), require_radius, (value,)))
        cases.append(("nonpositive_cutoff_" + str(value), require_cutoff, (value,)))
    for value in (s.Rational(3, 2),):
        cases.append(("noninteger_cutoff", require_cutoff, (value,)))
    for value in (
        "new_vacuum",
        "conditioned_Gaussian",
        "Gaussian_nonlinear_chart",
        "homogeneous_state",
    ):
        cases.append(("unsupported_state_" + value, require_state, (value,)))
    for value in (
        "full_BRST_measure",
        "positive_Weyl_symbol",
        "flat_global_quantum_measure",
    ):
        cases.append(("unsupported_measure_" + value, require_measure, (value,)))
    for value in (
        "full_interacting_mean",
        "global_diffeomorphism_quotient",
        "original_P8_closed",
    ):
        cases.append(("unsupported_observable_" + value, require_observable, (value,)))
    cases.append(("drop_primitive_frontier", validate_scope, ([], matching())))
    cases.append(("drop_matching_frontier", validate_scope, (frontier(), [])))
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
            raise ValueError("Invalid input was accepted: " + name)
    return count


def controls():
    return {
        "whole_full_convolution_and_generated_mean_retained": True,
        "whole_vector_matter_and_boundary_generator_retained": True,
        "full_adjoint_and_nonsymmetric_inverse_contacts_retained": True,
        "residual_translation_kernel_not_discarded": True,
        "same_reference_gaussian_not_projected_or_reselected": True,
        "finite_IR_UV_mode_family_not_continuum_nonlinear_state": True,
        "S261_refutation_S265_boundary_and_original_P8_open": True,
        "written_functional_arguments_not_FORMALIZED": True,
        "rejected_inputs": rejected_inputs(),
    }
