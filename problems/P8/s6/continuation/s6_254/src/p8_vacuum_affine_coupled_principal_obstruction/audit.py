"""The current off-reference obstruction does not promote original P8 rows."""

from functools import cache

import sympy as s
from p8_vacuum_affine_physical_background_vertices import audit as previous

from . import background, coupled, principal

STATE = previous.STATE
OBSERVABLES = (
    "complete_homogeneous_off_reference_four_mode_Gaussian_symbol",
    "all_eight_physical_modes_fast_characteristic",
    "fixed_nonreference_time_dependent_unbounded_momentum_growth",
)
ITEM = {
    "id": "QG2_H8A420_complete_homogeneous_off_reference_Gaussian_fast_principal_and_written_time_dependent_growth_obstruction",
    "status": "EXACT_COMPLETE_HOMOGENEOUS_GAUSSIAN_FAST_SYMBOL_AND_WRITTEN_FIXED_NONREFERENCE_GROWTH_OBSTRUCTION_NOT_FINITE_CUTOFF_ON_SHELL_NONLINEAR_QUANTUM_UV_REGGE_OR_ORIGINAL_P8",
}
POSITIVE_PARAMETERS = (
    coupled.D,
    coupled.Z,
    coupled.J,
    coupled.K,
    coupled.a,
    coupled.zeta,
    coupled.Y,
    coupled.Yv,
)
DOMAIN_PARAMETERS = (*POSITIVE_PARAMETERS, coupled.r, coupled.L2)


def require_state(label):
    if not isinstance(label, str) or label != STATE:
        raise ValueError(
            "No new state, profile, finite matching prescription or parent"
        )
    return label


def require_observable(label):
    if not isinstance(label, str) or label not in OBSERVABLES:
        raise ValueError("Keep the fixed nonreference homogeneous Gaussian scope")
    return label


def exact_scalar(value):
    if isinstance(value, (bool, float, str)) or not isinstance(value, (int, s.Expr)):
        raise TypeError("An exact finite real scalar is required")
    value = s.sympify(value)
    if (
        value.free_symbols
        or value.has(s.Float)
        or value.is_real is not True
        or value.is_finite is not True
    ):
        raise ValueError("An exact finite real scalar is required")
    return value


def require_domain(values):
    if not isinstance(values, dict) or set(values) != set(DOMAIN_PARAMETERS):
        raise ValueError("The complete exact nonreference domain is required")
    result = {key: exact_scalar(value) for key, value in values.items()}
    if not all(result[key].is_positive is True for key in POSITIVE_PARAMETERS):
        raise ValueError("Every positive domain coefficient must be strictly positive")
    if (
        result[coupled.r].is_zero is not False
        or result[coupled.L2].is_zero is not False
    ):
        raise ValueError("The nonreference fast proof requires r and L2 nonzero")
    gamma = s.factor(coupled.GAMMA.subs(result, simultaneous=True))
    if gamma.is_positive is not True:
        raise ValueError(
            "The joint velocity and constraint chart requires gamma positive"
        )
    return result


def example_domain():
    return {
        **dict.fromkeys(POSITIVE_PARAMETERS, s.S.One),
        coupled.r: s.Rational(1, 10),
        coupled.L2: -s.S.One,
    }


def frontier():
    return previous.frontier()


def matching():
    return [*previous.matching(), dict(ITEM)]


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("No original or prior matching obligation is promoted")
    return True


@cache
def packets():
    return {
        "whole_current_background_and_heavy_source": background.data(),
        "whole_four_mode_joint_constraint_reduction": coupled.reduction(),
        "whole_time_dependent_physical_central_chart": coupled.chart(),
        "whole_eight_mode_fast_symbol_and_exact_remainder": principal.fast_data(),
        "independent_whole_rank_one_determinant_diagnostic": principal.determinant_data(),
        "whole_time_dependent_growth_normal_form": principal.growth_data(),
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
        "no_finite_cutoff_or_on_shell_exclusion_inferred": True,
        "original_P8_remains_open": True,
        "fixed_affine_clock_chart_not_every_covariant_offshell_Hessian": True,
    }


def observable():
    return {
        "fixed_chart": "This is the second variation in S253's fixed affine-clock fluctuation chart. A nonlinear field redefinition can add its second derivative times a nonzero background onepoint. No exclusion of every covariant or connection-corrected off-shell Hessian is claimed; the exact time-dependent canonical transformations within the specified quadratic system retain its full operator.",
        "scope": "The unchanged QG2-H8A420 parent on homogeneous aligned W=S backgrounds with Hbar0, including the ENTIRE off-clock heavy source and both fixed profiles. All four coupled physical Gaussian modes and four remaining tensor/transverse-vector modes are retained. The whole exact time-dependent canonical generator has a P^(3/2) fast real branch on small fixed nonreference backgrounds; its written smooth-chart cone proof excludes an all-momentum finite-derivative Sobolev bound on those prescribed backgrounds.",
        "order_of_limits": "Fix a nonzero smooth background displacement and compact interval with strict positive margins first; then let momentum tend to infinity. Constants and the required momentum threshold need not be uniform as the displacement tends to zero or the heavy mass grows. No exchange with the reference limit or with a Wilsonian cutoff is made.",
        "boundary": "Not an unforced on-shell nonlinear bounce, not the nonlinear or quantum corrected principal operator, not a sub-cutoff instability claim, not an exclusion of the original witness or every UV completion. The earlier finite probe jets, reference Gaussian and finite Fourier-ball inverses retain their scopes. No physical curved subtraction, loop norm, interacting state, UV scattering, finite-gravity IR/Regge or original V/G/B/P8 closure is established.",
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
    for i, value in enumerate(invalid):
        cases.append(("invalid_state_" + str(i), require_state, (value,)))
        cases.append(("invalid_observable_" + str(i), require_observable, (value,)))
        for key in DOMAIN_PARAMETERS:
            data = example_domain()
            data[key] = value
            cases.append(
                ("invalid_domain_" + str(key) + "_" + str(i), require_domain, (data,))
            )
    for key in POSITIVE_PARAMETERS:
        for value in (0, -1):
            data = example_domain()
            data[key] = value
            cases.append(
                (
                    "nonpositive_domain_" + str(key) + "_" + str(value),
                    require_domain,
                    (data,),
                )
            )
    for key in (coupled.r, coupled.L2):
        data = example_domain()
        data[key] = 0
        cases.append(("missing_fast_margin_" + str(key), require_domain, (data,)))
    for value in (1, s.sqrt(s.Rational(2, 3))):
        data = example_domain()
        data[coupled.r] = value
        cases.append(("nonpositive_gamma_" + str(value), require_domain, (data,)))
    for label in (
        "original_P8_closed",
        "finite_cutoff_instability_proved",
        "unforced_bounce_refuted",
        "new_parent_allowed",
        "heavy_source_is_zero_off_clock",
        "frozen_eigenvalues_alone_prove_growth",
        "all_momentum_reference_bound_continues_to_neighborhood",
        "quantum_loops_controlled",
        "interacting_state_constructed",
        "uniform_as_r_tends_to_zero",
        "indefinite_Hpp_alone_is_a_ghost",
        "current_parent_has_no_UV_completion",
        "every_covariant_offshell_Hessian_excluded",
    ):
        cases.append(("unsupported_" + label, require_observable, (label,)))
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
    data = example_domain()
    del data[coupled.K]
    cases.append(("missing_domain_parameter", require_domain, (data,)))
    data = example_domain()
    data[coupled.H] = 0
    cases.append(("extra_domain_parameter", require_domain, (data,)))
    cases.append(("missing_primitive", validate_scope, (frontier()[:-1], matching())))
    cases.append(("missing_matching", validate_scope, (frontier(), matching()[:-1])))
    return cases


@cache
def rejected_inputs():
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (ValueError, TypeError):
            continue
        raise ValueError("Unsupported principal or scope input accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "complete_heavy_source_not_deleted_off_clock": background.data()["gates"][
            "entire_heavy_source_is_not_identically_zero"
        ],
        "all_four_coupled_modes_required": coupled.PHASE.shape == (8, 1),
        "all_sixteen_phases_in_fast_symbol": principal.fast_data()[
            "whole_sixteen_phase_fast_symbol"
        ].shape
        == (16, 16),
        "whole_time_connections_required": True,
        "whole_regular_remainder_required": True,
        "slow_nilpotent_block_not_assumed_diagonalizable": True,
        "fixed_background_limit_not_finite_cutoff_verdict": True,
        "all_original_and_prior_matching_rows_unchanged": frontier()
        == previous.frontier()
        and matching()[:-1] == previous.matching(),
    }
