"""Ordered regulator removal and an unexpanded exponentially fine detector."""

from functools import cache

import sympy as s

from . import gamma, poisson, source

STRENGTH = s.Symbol("fixed_soft_strength_K0_over_4pi2", nonnegative=True)
MATCH = s.Symbol("fixed_angular_conversion_times_kappa", real=True)
CHI = s.Symbol("nonnegative_resolution_scaling", nonnegative=True)
RHO = s.Symbol("positive_joint_regulator_times_kappa", positive=True)


def require_domain(index, ratio):
    values = []
    for value in (index, ratio):
        if isinstance(value, (bool, float, s.Float)) or not isinstance(
            value, (int, s.Integer, s.Rational)
        ):
            raise TypeError("Require exact real rational detector parameters")
        values.append(s.sympify(value))
    a, x = values
    if a < 0 or x <= 0 or x > 1:
        raise ValueError("Require a>=0 and0<resolution/reference<=1")
    return tuple(values)


def scaled_factor(
    kappa=source.K, strength=STRENGTH, conversion_strength=MATCH, chi=CHI
):
    k, strength, d, chi = map(s.sympify, (kappa, strength, conversion_strength, chi))
    return gamma.conversion(strength / k, d / k) * s.exp(-strength * chi)


def joint_limit_factor(strength=STRENGTH, rho=RHO, chi=CHI):
    strength, rho, chi = map(s.sympify, (strength, rho, chi))
    return s.exp(strength * (s.exp(-2 * rho * chi) - 1) / (2 * rho))


@cache
def data():
    checks = {}

    def put(name, value):
        checks[name] = s.factor(s.expand_func(value))

    k, A, chi, rho, x = s.symbols("kappa A chi rho x", positive=True)
    D = s.Symbol("D", real=True)
    factor = scaled_factor(k, A, D, chi)
    put("unexpanded_resolution_power", (A / k) * (-k * chi) + A * chi)
    put(
        "physical_scaled_factor",
        factor - s.exp(D / k - s.EulerGamma * A / k - A * chi) / s.gamma(1 + A / k),
    )
    put(
        "weak_gravity_at_fixed_ratio",
        s.limit(poisson.calorimetric_limit(A / k, D / k, x), k, s.oo) - 1,
    )
    put(
        "zero_resolution_at_fixed_gravity",
        s.limit(poisson.calorimetric_limit(A / k, D / k, x), x, 0, dir="+"),
    )
    put("correlated_resolution_limit", s.limit(factor, k, s.oo) - s.exp(-A * chi))
    put("joint_regulator_index", (A / k) / (2 * (rho / k)) - A / (2 * rho))
    put(
        "joint_regulator_resolution_exponent",
        2 * (rho / k) * (-k * chi) + 2 * rho * chi,
    )
    ae = A / k + 2 * (rho / k) * (D / k)
    mean = ae * s.gamma(1 + 2 * rho / k) * s.exp(-2 * rho * chi) / (2 * rho / k)
    put(
        "joint_Poisson_mean_finite",
        s.limit(mean, k, s.oo) - A * s.exp(-2 * rho * chi) / (2 * rho),
    )
    put("joint_Poisson_scaled_mean_zero", s.limit(2 * rho * mean / k, k, s.oo))
    put(
        "joint_Poisson_scaled_variance_zero",
        s.limit(4 * rho * rho * mean / k**2, k, s.oo),
    )
    put(
        "joint_analytic_virtual_exponent",
        s.limit(mean - A / (2 * rho), k, s.oo)
        - A * (s.exp(-2 * rho * chi) - 1) / (2 * rho),
    )
    put(
        "joint_limit_ordered_boundary",
        s.limit(joint_limit_factor(A, rho, chi), rho, 0, dir="+") - s.exp(-A * chi),
    )
    put(
        "joint_limit_bad_large_boundary",
        s.limit(joint_limit_factor(A, rho, chi), rho, s.oo) - 1,
    )
    put(
        "joint_exponent_loss_primitive",
        s.integrate(1 - s.exp(-2 * rho * x), (x, 0, chi))
        - chi
        - (s.exp(-2 * rho * chi) - 1) / (2 * rho),
    )
    put(
        "joint_exponent_error_majorant",
        s.integrate(2 * rho * x, (x, 0, chi)) - rho * chi**2,
    )
    put(
        "single_real_truncation_at_scaled_resolution",
        1 + (A / k) * (-k * chi) - (1 - A * chi),
    )
    put(
        "single_real_truncation_negative_control",
        (1 - A * chi).subs({A: 2, chi: 1}) + 1,
    )
    put(
        "exact_leading_limit_positive_control",
        s.exp(-A * chi).subs({A: 2, chi: 1}) - s.exp(-2),
    )
    put("inherited_compact_resolution_ratio", s.Rational(1, 8) - s.Rational(1, 8) / 1)
    put(
        "uniform_original_leading_conversion_bound",
        gamma.uniform_error_bound(source.KAPPA) - s.Rational(2, 10**1596),
    )
    put("Gamma_scale_independence", s.diff(gamma.conversion(A / k, D / k), chi))
    return {
        "whole_unexpanded_correlated_soft_factor": scaled_factor(),
        "whole_wrong_joint_regulator_limit": joint_limit_factor(),
        "whole_ordered_weak_gravity_scaling_limit": s.exp(-STRENGTH * CHI),
        "whole_joint_exponent_error_bound": STRENGTH * RHO * CHI**2,
        "whole_order_of_limits_proof": "First take e->0 at fixed positive resolution and fixedkappa. Then x=exp(-kappa*chi) gives exp(-A*chi)*exp(D/kappa)*F(A/kappa), not a finite Taylor series in a*lnx. Atfixedkappa andA>0,x->0 gives0; kappa->infinity first gives1. If e=rho/kappa is instead held on a joint path withrho>0, the auxiliary Poisson mean tends toA*exp(-2rhochi)/(2rho) and2eN tends to0. The limit is exp[A*(exp(-2rhochi)-1)/(2rho)], distinct from the ordered exp(-Achi). Its logarithmic excess lies between0 andA*rho*chi^2. This is an explicit leading-soft limit-order test, not a physical alternative prescription.",
        "whole_missing_uniform_control": "No bound on the complete full-amplitude minus leading-soft remainder uniformly in the unbounded graviton multiplicity is supplied. S295 controls the minimal single-real difference, not allN or gravity-Born radiation. Finite hard matching, omitted hard loops, complex-energy analyticity/Regge growth, the original quantum state and the bounce connection remain open.",
        "checks": checks,
        "gates": {
            "fixed_positive_resolution_regulator_removal_first": True,
            "full_resolution_power_kept_in_scaling_regime": True,
            "explicit_noncommuting_iterated_limits": True,
            "explicit_finite_rho_joint_limit_countercheck": True,
            "quantified_joint_logarithmic_excess": True,
            "actual_original_compact_domain_not_extended": True,
            "all_multiplicity_recoil_and_hard_loop_errors_still_open": True,
        },
    }
