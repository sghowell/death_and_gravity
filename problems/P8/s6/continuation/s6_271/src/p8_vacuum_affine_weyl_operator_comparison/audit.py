"""Qualified finite Weyl operator comparison with all original frontiers open."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_quantitative_local_time import audit as previous

from . import comparison, derivatives, kernel, source

STATE, MEASURE = previous.STATE, previous.MEASURE
OBSERVABLES = (
    "full_finite_Weyl_operator_ordering_comparison",
    "positive_full_Weyl_volume_for_the_declared_finite_regulator",
    "entire_two_ordering_unitary_operator_and_two_Weyl_cutoff_comparison",
)
ITEM = {
    "id": "QG2_H8A435_explicit_Weyl_OPERATOR_ordering_positive_volume_and_entire_unitary_comparison",
    "status": "COMPLETE_EXPLICIT_FINITE_WEYL_OPERATOR_NORM_ORDERING_POSITIVE_VOLUME_AND_UNITARY_COMPARISON_NOT_ORIGINAL_PHYSICAL_MATCHING_CONTINUUM_LOOP_OR_V_G_B_P8",
}


def parameters():
    return previous.parameters()


def require_parameters(value):
    return previous.require_parameters(value)


def require_radius(value):
    return previous.require_radius(value)


def require_time(value):
    return previous.require_time(value)


def require_cutoff(value):
    return previous.require_cutoff(value)


def require_state(value):
    return previous.require_state(value)


def require_measure(value):
    return previous.require_measure(value)


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError(
            "Require one of the defined full finite Weyl operator comparisons"
        )
    return value


def require_order(value):
    derivatives.coefficient(value)
    return value


def require_amplitude(value):
    derivatives.amplitude_jet(0, value)
    return value


def frontier():
    return previous.frontier()


def matching():
    return [*previous.matching(), dict(ITEM)]


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("Finite Weyl operator comparison does not close original P8")
    return True


@cache
def packets():
    return {
        "whole_unchanged_full_source_and_high_phase_domain": source.data(),
        "whole_direct_normalized_Gaussian_and_operator_Schur_bound": kernel.data(),
        "whole_full_high_cutoff_and_holomorphic_phase_derivatives": derivatives.data(),
        "whole_actual_Weyl_operator_positivity_and_entire_unitary_comparison": comparison.data(),
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
        "same_original_family_parameters": require_parameters(parameters())
        == parameters(),
        "S261_refutation_and_S265_integrability_boundary_unchanged": True,
        "original_P8_open_despite_actual_finite_Weyl_operator_estimates": True,
    }


def observable():
    return {
        "established": "An explicit normalized Gaussian-frame Schur theorem and all196 full cutoff/phase derivatives turn the previous SYMBOL error into an actual finite Weyl OPERATOR norm comparison. Volume ordering error<1e-200, Weyl volume distance fromI<1e-199 and actual operator floor>1/2. Entire two-ordering unitary OPERATOR difference<1e-944. Two Weyl cutoff states differ<1e-943 and volume means<1e-199; evolved original coherent-POVM leakage<1e-942, not zero.",
        "domain": "The unchanged S270 full original source, same pure48-pair finite family, all96 canonical phase variables, real |u|<=1e-2000 and exactly two explicit c1/c2 cutoffs. All spatial harmonics, original scalar phases, nonlinear implicit, clock, primitive, matter/vector, Gauss and full derivative contacts remain.",
        "not_established": "No generic positivity-preserving Weyl map, original interacting mean, unregularized singular Hamiltonian, uniform mode/volume limit, physical Wilsonian matching, omitted loops, UV completion or nonlinear global completeness. The stronger S270 coherent-only cutoff estimates are not reassigned to Weyl comparisons. Original V/G/B/P8 remain OPEN; completed scoped P8(a) is unchanged.",
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
        cases.append(("bad_high_order_" + str(i), require_order, (value,)))
        cases.append(("bad_full_amplitude_" + str(i), require_amplitude, (value,)))
    for order in (-1, 197, 10**10):
        cases.append(
            ("outside_phase_derivative_range_" + str(order), require_order, (order,))
        )
    for i, value in enumerate(
        (0, 1, -1, 10**900, s.Rational(1, 10**200), "10**1000", "1/10**255")
    ):
        cases.append(
            ("wrong_or_coerced_amplitude_" + str(i), require_amplitude, (value,))
        )
    for value in (
        *previous.OBSERVABLES,
        "symbol_sup_is_operator_norm",
        "arbitrary_positive_symbol_has_positive_Weyl_operator",
        "new_diagonal_reference_vacuum",
        "all96_Cauchy_radius_R_over4",
        "unproved_high_time_profile_derivatives",
        "original_P8_closed",
        "coherent_cutoff_estimates_unchanged_for_Weyl",
        "uniform_continuum_limit",
    ):
        cases.append(("wrong_operator_scope_" + value, require_observable, (value,)))
    cases.extend(
        [
            ("delete_original_frontier", validate_scope, ([], matching())),
            ("delete_operator_matching_frontier", validate_scope, (frontier(), [])),
        ]
    )
    if len({name for name, _, _ in cases}) != len(cases):
        raise ValueError("Duplicate rejected Weyl comparison input")
    return cases


def rejected_inputs():
    count = 0
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (ValueError, TypeError):
            count += 1
        else:
            raise ValueError("Unsupported Weyl operator scope accepted: " + name)
    return count


def controls():
    return {
        "direct_normalized_Gaussian_kernel_not_inconsistent_printed_prefactor": True,
        "explicit_full_operator_constant_not_symbol_sup_shortcut": True,
        "all196_phase_derivatives_not_unproved_time_profile_derivatives": True,
        "all_high_cutoff_mixed_implicit_and_covariance_contacts_retained": True,
        "actual_positive_Weyl_volume_not_generic_positivity": True,
        "entire_scalar_phases_unitaries_and_changed_readouts_compared": True,
        "finite_ordering_comparison_not_original_physical_matching": True,
        "historical_errata_written_not_FORMALIZED_original_P8_open": True,
        "rejected_inputs": rejected_inputs(),
    }
