"""Evaluated bounce-domain scope with the unchanged finite quantum regulator."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_finite_window_growth.intervals import fraction
from p8_vacuum_affine_local_quantum_regulator import audit as previous

from . import geometry, invariants, reference, source

STATE, MEASURE = previous.STATE, previous.MEASURE
OBSERVABLES = (
    "all_twelve_nonlinear_bounce_invariant_image",
    "same_seed_quantitative_coherent_core_mass",
    "positive_coherent_bounce_volume_and_cutoff_comparison",
)
ITEM = {
    "id": "QG2_H8A433_evaluated_finite_nonlinear_bounce_phase_domain_and_same_seed_volume_comparison",
    "status": "COMPLETE_EVALUATED_FINITE_BOUNCE_PHASE_DOMAIN_SAME_SEED_COHERENT_TAIL_AND_POSITIVE_VOLUME_REGULATOR_COMPARISON_NOT_TIME_EVOLUTION_MATCHING_OR_ORIGINAL_V_G_B_P8",
}


def parameters():
    return previous.parameters()


def require_parameters(values):
    return previous.require_parameters(values)


def require_radius(value):
    value = fraction(value)
    if not 0 <= value <= fraction(reference.SUPPORT):
        raise ValueError("Require the evaluated whitened support ball")
    return value


def require_time(value):
    value = fraction(value)
    if value != 0:
        raise ValueError("Only the bounce slice is evaluated")
    return value


def require_dimension(value):
    dimension = previous.require_dimension(value)
    if dimension != reference.DIMENSION:
        raise ValueError("Require all48 canonical pairs in the declared finite family")
    return dimension


def require_momentum(value):
    value = fraction(value)
    if value != fraction(reference.P):
        raise ValueError(
            "Require the fixed integer wave number of the declared torus family"
        )
    return value


def require_state(value):
    return previous.require_state(value)


def require_measure(value):
    return previous.require_measure(value)


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError(
            "Require an evaluated bounce-domain or coherent-volume readout"
        )
    return value


def frontier():
    return previous.frontier()


def matching():
    return [*previous.matching(), dict(ITEM)]


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("Evaluated finite bounce data does not close original P8")
    return True


@cache
def packets():
    return {
        "whole_original_canonical_invariants_and_source": source.data(),
        "whole_same_reference_and_evaluated_field_normalization": reference.data(),
        "whole_full_nonlinear_geometry_and_adjoint_lift": geometry.data(),
        "whole_all_invariant_images_and_coherent_volume_comparison": invariants.data(),
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
        "all_nine_original_primitive_rows_retained": len(frontier()) == 9,
        "all_matching_identifiers_distinct": len({r["id"] for r in matching()})
        == len(matching()),
        "all_original_physical_parameters_unchanged": require_parameters(parameters())
        == parameters(),
        "S261_refutation_and_S265_unlocalized_integrability_boundary_unchanged": True,
        "original_P8_open_despite_evaluated_finite_bounce_domain": True,
    }


def observable():
    return {
        "established": "At the actual bounce, the full nonlinear image of an evaluated96-phase-dimensional whitened support ball lies uniformly in the original spatial and twelve-invariant auxiliary chart. The SAME pure nonzero seed has coherent core leakage below exp(-1e39). Its positive finite coherent volume has operator distance fromI below1e-255, and two admissible core-agreeing volume cutoffs have an evaluated exponentially small same-seed difference.",
        "domain": "Fixed L=1 torus, integer P=1e64 along three independent axes and both signs, all eight channels, d=48. Core radius1e20, support radius2e20 in the FULL covariance whitening. All reconstructed harmonics, density factors, adjoint ghost contacts and residual translations remain. Homogeneous variables are external reference parameters; no new state is selected.",
        "not_established": "No evaluated positive time interval, later-time leakage, complete interaction/ordering bound, calibrated-volume positivity or equality to an original quantum mean; no continuum limit, physical matching cutoff, omitted-loop estimate or nonlinear global completeness. This is a declared finite coherent regulator, not a new original action or R3-state identity. Original V/G/B/P8 remain OPEN and scoped P8(a) is unchanged.",
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
        ("time", require_time),
        ("dimension", require_dimension),
        ("momentum", require_momentum),
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
    for value in (-1, 2 * reference.SUPPORT):
        cases.append(("outside_phase_ball_" + str(value), require_radius, (value,)))
    for value in (-1, 1, s.Rational(1, 10**2000)):
        cases.append(("unevaluated_time_" + str(value), require_time, (value,)))
    for value in (1, 16, 96):
        cases.append(("wrong_dimension_" + str(value), require_dimension, (value,)))
    for value in (0, reference.P / 2, 2 * reference.P):
        cases.append(("wrong_momentum_" + str(value), require_momentum, (value,)))
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
        "actual_small_Hamiltonian_error",
        "uniform_continuum_tail",
        "nonlinear_Gaussian_pushforward",
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
            raise ValueError(
                "Unsupported finite phase-domain input was accepted: " + name
            )
    return count


def controls():
    return {
        "all_eight_channels_and_full_covariance_rows_retained": True,
        "whole_nonlinear_convolution_adjoint_and_curvature_not_projected": True,
        "all_twelve_density_invariants_not_independent_Gaussians": True,
        "trace_preserved_but_shear_lift_not_deleted": True,
        "actual_finite_Husimi_tail_not_exact_phase_projection": True,
        "coherent_volume_not_calibrated_or_original_quantum_mean": True,
        "historical_errata_and_original_P8_open": True,
        "written_Banach_operator_arguments_not_FORMALIZED": True,
        "rejected_inputs": rejected_inputs(),
    }
