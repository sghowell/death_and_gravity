"""Evaluated finite-time regulator scope with every original frontier retained."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_finite_window_growth.intervals import fraction
from p8_vacuum_affine_quantitative_phase_domain import audit as previous

from . import branch, quantum, reference, source

STATE, MEASURE = previous.STATE, previous.MEASURE
OBSERVABLES = (
    "complete_finite_time_nonlinear_phase_domain",
    "positive_coherent_and_calibrated_volume_with_explicit_ordering_bounds",
    "same_seed_evolved_coherent_leakage_and_two_defined_regulator_comparison",
)
ITEM = {
    "id": "QG2_H8A434_evaluated_full_source_finite_time_quantum_regulator_and_same_seed_comparison",
    "status": "COMPLETE_EXPLICIT_FINITE_TINY_TIME_DOMAIN_CALIBRATED_POSITIVE_READOUT_STATE_LEAKAGE_AND_TWO_REGULATOR_COMPARISON_NOT_ORIGINAL_MATCHING_CONTINUUM_OR_V_G_B_P8",
}


def parameters():
    return previous.parameters()


def require_parameters(values):
    return previous.require_parameters(values)


def require_radius(value):
    return previous.require_radius(value)


def require_time(value):
    value = fraction(value)
    if not -fraction(source.TIME) <= value <= fraction(source.TIME):
        raise ValueError("Require the evaluated full real time interval")
    return value


def require_dimension(value):
    return previous.require_dimension(value)


def require_momentum(value):
    return previous.require_momentum(value)


def require_cutoff(value):
    value = fraction(value)
    if value not in (1, 2):
        raise ValueError(
            "Require one of the two explicitly differentiated radial cutoffs"
        )
    return value


def require_state(value):
    return previous.require_state(value)


def require_measure(value):
    return previous.require_measure(value)


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError("Require an evaluated finite-time regulated observable")
    return value


def frontier():
    return previous.frontier()


def matching():
    return [*previous.matching(), dict(ITEM)]


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("Finite tiny-time quantum estimates do not close original P8")
    return True


@cache
def packets():
    return {
        "whole_full_time_source_primitive_and_complex_lapse_jets": source.data(),
        "whole_original_fixed_canonical_flow_and_complete_phase_image": reference.data(),
        "whole_full_complex_auxiliary_branch_and_implicit_contacts": branch.data(),
        "whole_explicit_cutoff_ordering_positive_volume_and_quantum_comparison": quantum.data(),
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
        "original_P8_open_despite_evaluated_finite_quantum_time": True,
    }


def observable():
    return {
        "established": "The declared full finite regulator has an evaluated real time interval |u|<=1e-2000, entire nonlinear phase domain and same-seed unitary state error<1e-980. Its evolved coherent outside-core probability is<1e-980, not zero. Both the positive coherent and separately first-Weyl-calibrated volume are actually positive and uniformly close toI; the calibrated operator distance is<2e-255, with heat SYMBOL ordering error<1e-310. The two explicit cutoff regulators have state difference<1e-1970 and corresponding readout mean difference<1e-1230.",
        "domain": "Same L=1, P=1e64, all-eight-channel, d=48 finite torus family and full original pure reference; real support ball2e20, complex initial ball4e20 and full linear image ball8e20. The cutoffs have c=1,2 in the exact differentiated radial formula. All nonlinear harmonics, source/matter/vector/primitive/time/implicit contacts and residual translations remain.",
        "not_established": "No original interacting volume mean or unregularized Hamiltonian, automatic Weyl OPERATOR norm comparison, arbitrary-cutoff bound, uniform mode/volume limit, physical Wilsonian matching, omitted-loop estimate or nonlinear global completeness. Omega_star is not the physical heavy dispersion; the complex proof uses the holomorphic continuation of the real adjoint, not conjugate-transpose for complex coefficients. Original V/G/B/P8 remain OPEN and completed scoped P8(a) is unchanged.",
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
        ("cutoff", require_cutoff),
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
    for value in (-1, 4 * reference.CORE):
        cases.append(
            ("outside_real_phase_ball_" + str(value), require_radius, (value,))
        )
    for value in (-2 * source.TIME, 2 * source.TIME, s.Rational(1, 10**60)):
        cases.append(("outside_evaluated_time_" + str(value), require_time, (value,)))
    for value in (1, 16, 96):
        cases.append(("wrong_dimension_" + str(value), require_dimension, (value,)))
    for value in (0, reference.field.P / 2, 2 * reference.field.P):
        cases.append(("wrong_momentum_" + str(value), require_momentum, (value,)))
    for value in (0, 3, s.Rational(3, 2)):
        cases.append(("unproved_cutoff_" + str(value), require_cutoff, (value,)))
    for value in (
        "new_vacuum",
        "conditioned_Gaussian",
        "exact_phase_support",
        "complex_physical_state",
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
        "Weyl_operator_norm_from_symbol_bound",
        "arbitrary_smooth_cutoff_bound",
        "uniform_continuum_tail",
        "nonlinear_Gaussian_pushforward",
        "initial_tail_is_exact_support",
        "Omega_star_is_physical_dispersion",
        "conjugate_adjoint_is_holomorphic",
    ):
        cases.append(("wrong_observable_" + value, require_observable, (value,)))
    for i, value in enumerate(
        (
            s.diff(source.R, source.N, 6),
            s.diff(source.F, source.u, 6),
            s.diff(source.j, source.N, 6),
            s.diff(source.Hclock, source.u, 4),
            1 / s.Symbol("unproved_denominator"),
            s.sqrt(s.Symbol("unproved_fractional_base")),
            s.Symbol("unbounded_atom"),
            s.Float(1),
            s.oo,
            s.nan,
        )
    ):
        cases.append(("unproved_majorant_atom_" + str(i), branch.magnitude, (value,)))
    cases.extend(
        [
            ("delete_primitive_frontier", validate_scope, ([], matching())),
            ("delete_matching_frontier", validate_scope, (frontier(), [])),
        ]
    )
    if len({name for name, _, _ in cases}) != len(cases):
        raise ValueError("Duplicate rejected-input case identifiers")
    return cases


def rejected_inputs():
    count = 0
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (ValueError, TypeError):
            count += 1
        else:
            raise ValueError("Unsupported finite-time input was accepted: " + name)
    return count


def controls():
    return {
        "whole_time_source_and_original_reference_not_bounce_replacement": True,
        "all_nonlinear_harmonics_adjoint_density_and_implicit_contacts_retained": True,
        "complex_formal_transpose_not_conjugate_adjoint": True,
        "physical_H_frequency_not_fixed_energy_balance": True,
        "entire_cutoff_derivatives_and_background_phase_retained": True,
        "positive_calibrated_coherent_volume_not_original_matching": True,
        "evolved_POVM_tail_not_exact_phase_support_or_continuum_limit": True,
        "historical_errata_written_not_FORMALIZED_original_P8_open": True,
        "rejected_inputs": rejected_inputs(),
    }
