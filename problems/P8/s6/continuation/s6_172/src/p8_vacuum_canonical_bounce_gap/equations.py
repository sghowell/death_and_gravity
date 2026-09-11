"""Exact physical-frame null Einstein residual of a named minimal diagnostic."""

from functools import cache

import sympy as s
from p8_exceptional_vacuum import analytic
from p8_vacuum_flat_dirac_hadamard.symbols import rational

from . import nullstress


def margins(additional_absolute_null_bound):
    rest = rational(additional_absolute_null_bound)
    if rest < 0:
        raise ValueError("An absolute additional null-stress bound must be nonnegative")
    free = nullstress.data()
    k = analytic.KAPPA
    return {
        "central_interval_null_equation_residual_lower": s.Rational(121, 25) * k
        - free["full_uniform_absolute_free_null_stress_upper"]
        - rest,
        "actual_bounce_null_equation_residual_lower": 9 * k
        - free["actual_bounce_absolute_free_null_stress_upper"]
        - rest,
    }


@cache
def data():
    t = s.Symbol("t", real=True)
    k = s.Symbol("kappa", positive=True)
    b, N = s.symbols("log_a lapse", real=True)
    v, vd, Nd, phid, V = s.symbols("v vdot Ndot Phidot V", real=True)
    # Einstein and canonical scalar densities after one boundary integration.
    Lg = -3 * k * s.exp(3 * b) * v * v / N
    Lm = s.exp(3 * b) * (phid * phid / (2 * N) - N * V)
    dt = lambda f: s.diff(f, b) * v + s.diff(f, v) * vd + s.diff(f, N) * Nd
    rhoG = s.simplify(-s.diff(Lg, N).subs(N, 1) / s.exp(3 * b))
    PG = s.simplify(
        (s.diff(Lg, b) - dt(s.diff(Lg, v))).subs({N: 1, Nd: 0}) / (3 * s.exp(3 * b))
    )
    rhoM = s.simplify(-s.diff(Lm, N).subs(N, 1) / s.exp(3 * b))
    PM = s.simplify(s.diff(Lm, b).subs(N, 1) / (3 * s.exp(3 * b)))
    H = 4 * t / (1 + t * t)
    Hd = s.diff(H, t)
    residual = 2 * k * Hd + k
    free = nullstress.data()
    actual_k = analytic.KAPPA
    zero = margins(0)
    half = margins(actual_k / 2)
    repair = -9 * actual_k + free["actual_bounce_absolute_free_null_stress_upper"]
    return {
        "minimal_diagnostic_name": "CD-MIN-EQ: physical metric Einstein term -kappa R/2, tree-canonical clock Phi=sqrt(kappa)t, any positive-kinetic homogeneous classical scalars/healthy isotropic Yang-Mills, arbitrary derivative-free potential, and exactly the S6.171 free quadratic/reference fermion tensor.",
        "additional_terms_boundary": "All other quantum/source/finite-reference terms, nonminimal curvature, leading derivative interactions, higher operators and any change of gravity dictionary are collected in a signed additional null contribution. They are not asserted absent in the full GY14 candidate or the original DHOST target.",
        "classical_null_equation_residual": residual,
        "zero_additional_null_bound_margins": zero,
        "conditional_half_kappa_additional_bound_margins": half,
        "necessary_upper_on_additional_bounce_null_stress": repair,
        "necessary_upper_on_additional_bounce_null_stress_over_kappa": repair
        / actual_k,
        "meaning": "With the named free sector, the exact CD trajectory cannot solve this minimally Einstein-coupled canonical diagnostic. Any repair at the bounce requires an additional null contribution <=-9kappa+the displayed free remainder allowance, more negative if extra classical matter carries positive null stress. In particular a remainder bounded in absolute value by kappa/2 cannot repair it.",
        "not_a_full_candidate_exclusion": "No physical bound on all omitted interacting/gravity terms is supplied. The original CD action has leading nonminimal/derivative structure outside this diagnostic, and is not excluded. The result is a necessary signed budget for common-parent matching, not a new assumption imposed on every admissible P8 completion.",
        "checks": {
            "Einstein_lapse_energy_sign": s.factor(rhoG + 3 * k * v * v),
            "Einstein_scale_pressure_sign": s.factor(PG - k * (2 * vd + 3 * v * v)),
            "canonical_scalar_energy": s.factor(rhoM - phid * phid / 2 - V),
            "canonical_scalar_pressure": s.factor(PM - phid * phid / 2 + V),
            "null_potential_cancels": s.factor(rhoM + PM - phid * phid),
            "complete_classical_null_equation": s.factor(
                rhoG + PG + rhoM + PM - 2 * k * vd - phid * phid
            ),
            "actual_CD_Hdot": s.factor(Hd - 4 * (1 - t * t) / (1 + t * t) ** 2),
            "central_interval_Hdot_lower_factor": s.factor(
                Hd
                - s.Rational(48, 25)
                - 4 * (1 - 4 * t * t) * (13 + 3 * t * t) / (25 * (1 + t * t) ** 2)
            ),
            "actual_bounce_required_null_stress": (-2 * k * Hd).subs(t, 0) + 8 * k,
            "actual_clock_bounce_null_residual": residual.subs(t, 0) - 9 * k,
            "central_classical_lower_constant": 2 * s.Rational(48, 25)
            + 1
            - s.Rational(121, 25),
            "conditional_half_budget_slack": s.Rational(121, 25)
            - s.Rational(1, 2)
            - 4
            - s.Rational(17, 50),
        },
        "gates": {
            "free_only_central_residual_above_four_kappa": bool(
                zero["central_interval_null_equation_residual_lower"] > 4 * actual_k
            ),
            "free_only_bounce_residual_above_eight_kappa": bool(
                zero["actual_bounce_null_equation_residual_lower"] > 8 * actual_k
            ),
            "half_budget_central_residual_above_four_kappa": bool(
                half["central_interval_null_equation_residual_lower"] > 4 * actual_k
            ),
            "half_budget_bounce_residual_above_eight_kappa": bool(
                half["actual_bounce_null_equation_residual_lower"] > 8 * actual_k
            ),
            "necessary_extra_bounce_null_more_negative_than_eight_kappa": bool(
                repair < -8 * actual_k
            ),
            "conditional_additional_budget_not_claimed_derived": True,
        },
    }
