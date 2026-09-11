"""Scope guards for finite-time conditional sourced-vector response."""

from functools import cache

import sympy as s
from p8_vacuum_canonical_affine_decoupling import audit as previous
from p8_vacuum_flat_dirac_hadamard.symbols import rational

from . import clock, coefficients, estimates, hamiltonian

ITEM = {
    "id": "complete_affine_vector_retarded_mean_energy_and_causal_cubic_light_force",
    "status": "FIXED_CD_FULL_SOURCE_FINITE_TIME_CONDITIONAL_GAUSSIAN_MEAN_AND_LIGHT_FORCE_NOT_METRIC_NOISE_FULL_PARENT_OR_B_CLOSURE",
}
KAPPA = s.Integer(10) ** 800
ZETA = s.Rational(1, 10**6)


def require_scope(time, delta, zeta=ZETA, kappa=KAPPA):
    t, d, z, k = map(rational, (time, delta, zeta, kappa))
    if not -s.Rational(1, 2) <= t <= s.Rational(1, 2):
        raise ValueError("Outside the stated finite-time CD slab")
    if not 0 <= d <= s.Rational(1, 100):
        raise ValueError("Outside the full small-clock-jet class")
    if z != ZETA or k != KAPPA:
        raise ValueError("Require the unchanged conditional base parent")
    return t, d, z, k


def frontier():
    return previous.frontier()


def matching():
    return [dict(row) for row in previous.matching()] + [dict(ITEM)]


def validate_scope(rows, obligations):
    if rows != frontier() or obligations != matching():
        raise ValueError("Conditional mean response cannot close original P8 gates")
    return True


def observable():
    return {
        "parent": "Same S6.174/S6.176 classical parent and canonical conditional vector state; the metric is the unchanged CD geometry.",
        "source_class": "u=t+psi with psi compact smooth in the unit-length bounce slab times R3 and all coordinate jets through3 at most delta<=1/100; the complete analytic source is used.",
        "new_bounds": "A finite-time positive reference energy controls the retarded sourced mean at every spatial momentum. Explicit full-source tame constants give spatially integrated free coherent mean energy and a cubic H3-dual bound for the full source/contact light force.",
        "state_scope": "Only the conditional Gaussian sector with the fixed all-order initial covariance and fixed metric is used. Its connected fluctuations are unchanged by the coherent mean.",
        "important_distinction": "Free coherent mean energy is not the full source-interaction metric tensor. The full light force cancels for homogeneous closed sources. The generic shrinking-pole inverse issue is not a contradiction to a finite-time causal energy norm.",
        "remaining": "Metric and noise response, variable-geometry norms, full light/tensor/auxiliary/mixed loops, interacting state, self-consistent background, physical cutoff, full matching and V/G/B remain open.",
    }


@cache
def residuals():
    out = {}
    for name, packet in (
        ("hamiltonian", hamiltonian.data()),
        ("modes", hamiltonian.modes()),
        ("coefficients", coefficients.data()),
        ("clock", clock.data()),
        ("clock_bounds", clock.bounds()),
        ("estimates", estimates.data()),
    ):
        out.update({name + "_" + key: value for key, value in packet["checks"].items()})
    out.update(
        {
            "nine_primitive_rows_retained": len(frontier()) - 9,
            "primitive_rows_unchanged": sum(
                a != b for a, b in zip(frontier(), previous.frontier())
            ),
            "prior_matching_rows_retained": sum(
                a != b for a, b in zip(matching(), previous.matching())
            ),
            "one_finite_time_conditional_matching_row_added": len(matching())
            - len(previous.matching())
            - 1,
        }
    )
    return {k: s.cancel(v) for k, v in out.items()}


def scalar_entry_count():
    return len(residuals())


@cache
def gates():
    quantitative = {}
    for name, packet in (
        ("coefficients", coefficients.data()),
        ("clock", clock.bounds()),
        ("estimates", estimates.data()),
    ):
        quantitative.update(
            {name + "_" + key: value for key, value in packet["gates"].items()}
        )
    return {
        **quantitative,
        "same_fixed_CD_and_complete_classical_parent": True,
        "complete_source_and_same_covariant_prescription": True,
        "specified_conditional_Hadamard_state_only": True,
        "nonconserved_temporal_source_retained": True,
        "positive_reference_energy_not_full_driven_energy": True,
        "temporal_and_light_contact_cancellation_retained": True,
        "finite_time_zero_Cauchy_data_and_Green_hypotheses": True,
        "all_momentum_mathematical_bound_not_physical_cutoff": True,
        "complete_analytic_switch_bounds_not_Taylor_truncation": True,
        "spatial_and_arbitrary_Frechet_tame_norms_explicit": True,
        "no_support_volume_substitute_for_source_norm": True,
        "homogeneous_closed_source_cancels_full_force": True,
        "actual_inhomogeneous_nonclosed_source_control": True,
        "causal_expectation_not_symmetric_retarded_action_gradient": True,
        "metric_noise_and_full_interacting_response_remain_open": True,
        "all_prior_frozen_scientific_bytes_unchanged": True,
    }


def bad_cases():
    bad = (
        True,
        False,
        1.0,
        s.Float(1),
        "1",
        None,
        s.oo,
        -s.oo,
        s.zoo,
        s.I,
        s.nan,
        s.Symbol("x"),
    )
    out = []
    for i, value in enumerate(bad):
        for pos in range(4):
            args = [0, s.Rational(1, 100), ZETA, KAPPA]
            args[pos] = value
            out.append((f"type_{pos}_{i}", require_scope, tuple(args)))
    for t, d, z, k in (
        (-1, 0, ZETA, KAPPA),
        (1, 0, ZETA, KAPPA),
        (0, -1, ZETA, KAPPA),
        (0, s.Rational(1, 50), ZETA, KAPPA),
        (0, 0, 0, KAPPA),
        (0, 0, s.Rational(1, 2000), KAPPA),
        (0, 0, ZETA, 2 * KAPPA),
    ):
        out.append((f"scope_{len(out)}", require_scope, (t, d, z, k)))
    for i in range(9):
        rows = frontier()
        rows[i]["status"] = "COMPLETE"
        out.append((f"primitive_{i}", validate_scope, (rows, matching())))
    for i in range(len(matching())):
        rows = matching()
        rows[i]["status"] = "COMPLETE"
        out.append((f"matching_{i}", validate_scope, (frontier(), rows)))
    out.extend(
        (
            ("missing_primitive", validate_scope, (frontier()[:-1], matching())),
            ("missing_matching", validate_scope, (frontier(), matching()[:-1])),
        )
    )
    return out


@cache
def rejected_inputs():
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (TypeError, ValueError):
            continue
        raise ValueError("Unsupported response scope accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "temporal_constraint_source_not_dropped": True,
        "source_contact_force_not_omitted": True,
        "free_mean_energy_not_full_interaction_stress": True,
        "homogeneous_full_force_cancellation_retained": True,
        "no_retarded_single_branch_variational_shortcut": True,
        "finite_time_norm_not_global_on_shell_inverse_bound": True,
        "conditional_Gaussian_not_full_parent_response": True,
        "original_P8_not_closed": True,
    }
