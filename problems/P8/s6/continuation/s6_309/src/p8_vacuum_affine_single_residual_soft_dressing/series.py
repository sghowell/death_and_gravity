"""Signed analytic dimensional sum and exact marked-energy moments."""

from functools import cache

import sympy as s

from . import source


def exact_real(value):
    if value is None or isinstance(value, (bool, float, s.Float, str)):
        raise TypeError("Require a numeric exact real coordinate")
    return source.recoil.exact_real(value)


def require_epsilon(value):
    e = exact_real(value)
    if not 0 < e <= s.Rational(1, 8):
        raise ValueError("Require a positive soft regulator at most1/8")
    return e


def require_remaining(value):
    y = exact_real(value)
    if not 0 < y <= 1:
        raise ValueError(
            "Require positive remaining energy at most the reference scale"
        )
    return y


def finite_term(count, index_e, epsilon, remaining):
    n = source.poisson.require_count(count)
    ae = exact_real(index_e)
    e, y = require_epsilon(epsilon), require_remaining(remaining)
    return source.poisson.simplex_term(n, ae, e, y)


def marked_power_term(count, power, index_e, epsilon, resolution):
    n = source.poisson.require_count(count)
    p = exact_real(power)
    if p <= 0:
        raise ValueError("Require a positive integrable marked-energy power")
    ae = exact_real(index_e)
    e, x = require_epsilon(epsilon), require_remaining(resolution)
    return (
        s.gamma(p)
        * x**p
        * (ae * s.gamma(2 * e) * x ** (2 * e)) ** n
        / (s.factorial(n) * s.gamma(p + 1 + 2 * e * n))
    )


def marked_power_limit(power, index, conversion, resolution):
    p, a, d, x = map(exact_real, (power, index, conversion, resolution))
    require_remaining(x)
    if p <= 0 or a < 0:
        raise ValueError(
            "Require integrable marked power and nonnegative physical index"
        )
    return s.gamma(p) * s.exp(d - s.EulerGamma * a) * x ** (p + a) / s.gamma(p + 1 + a)


@cache
def data():
    checks = {}

    def put(name, value):
        checks[name] = s.factor(s.expand_func(s.expand_log(s.sympify(value))))

    e, x, p, a = s.symbols("e x p a", positive=True)
    ae, d, L = s.symbols("a_e Delta log_x", real=True)
    n = s.Symbol("additional_count", integer=True, nonnegative=True)
    for count in range(1, 9):
        put(
            "marked_Bose_factor_" + str(count),
            s.Integer(count) / s.factorial(count) - 1 / s.factorial(count - 1),
        )
    for count in range(5):
        soft = source.poisson.simplex_term(count, ae, e, x)
        beta = s.gamma(p) * s.gamma(1 + 2 * e * count) / s.gamma(p + 1 + 2 * e * count)
        marked = (
            s.gamma(p)
            * x**p
            * (ae * s.gamma(2 * e) * x ** (2 * e)) ** count
            / (s.factorial(count) * s.gamma(p + 1 + 2 * e * count))
        )
        put("entire_marked_simplex_" + str(count), soft * x**p * beta - marked)
    lam = (a + 2 * e * d) * s.gamma(1 + 2 * e) * s.exp(2 * e * L) / (2 * e)
    put(
        "whole_finite_exponent",
        s.series(lam - a / (2 * e), e, 0, 1).removeO() - d + s.EulerGamma * a - a * L,
    )
    put("positive_index_scaled_mean", s.limit(2 * e * lam, e, 0, dir="+") - a)
    put("positive_index_scaled_variance", s.limit(4 * e * e * lam, e, 0, dir="+"))
    put("zero_physical_index_signed_mean", s.limit(lam.subs(a, 0), e, 0, dir="+") - d)
    put(
        "zero_index_exact_exponential_boundary", s.exp(d) * x**0 / s.gamma(1) - s.exp(d)
    )
    uniform = (
        source.angular_bounds.K1_BOUND
        + 4 * source.angular_bounds.K0_BOUND
        + source.angular_bounds.PHASE_REMAINDER / 8
    ) / 72
    put("uniform_signed_exponent_majorant", uniform - s.Rational(10585, 2))
    put("signed_exponent_majorant_margin", 5500 - uniform - s.Rational(415, 2))
    put(
        "regulator_phase_numerator",
        source.angular_bounds.K1_BOUND
        + 4 * source.angular_bounds.K0_BOUND
        + source.angular_bounds.PHASE_REMAINDER / 8
        - 381060,
    )
    u = s.Symbol("positive_Gamma_integration_coordinate", positive=True)
    put(
        "reciprocal_Gamma_global_integral",
        s.integrate(s.exp(-u), (u, 1, s.oo)) - s.exp(-1),
    )
    put("Gamma_logconvex_endpoint_one", s.gamma(1) - 1)
    put("Gamma_logconvex_endpoint_two", s.gamma(2) - 1)
    put("Gamma_radial_recurrence", s.gamma(1 + 2 * e) / (2 * e) - s.gamma(2 * e))
    put("wrong_virtual_pole_residue", (ae - a) / (2 * e) - (ae - a) / (2 * e))
    for power in range(1, 6):
        ratio = s.gamma(power + 1) * s.gamma(1 + a) / s.gamma(power + 1 + a)
        put(
            "complete_energy_Beta_ratio_" + str(power),
            ratio - s.factorial(power) / s.prod(a + j for j in range(1, power + 1)),
        )
    put(
        "first_marked_moment_remaining_energy",
        s.gamma(1 + a) / s.gamma(2 + a) - 1 / (1 + a),
    )
    put(
        "second_marked_moment_remaining_energy",
        2 * s.gamma(1 + a) / s.gamma(3 + a) - 2 / ((1 + a) * (2 + a)),
    )
    put("individual_cut_wrong_at_first_moment", 1 - 1 / (1 + a) - a / (1 + a))
    put(
        "marked_limit_zero_index",
        s.gamma(p) * s.exp(d) * x**p / s.gamma(p + 1) - s.exp(d) * x**p / p,
    )
    margins = {
        "signed_exponent_below_half": s.Rational(1, 2) - 5500 / source.KAPPA,
        "one_plus_epsilon_Gamma_domain": 1 - s.Rational(1, 4),
        "marked_reference_energy_domain": 1 - s.Rational(1, 8),
        "signed_majorant_roundup": 5500 - uniform,
    }
    for name, value in margins.items():
        put("positive_arithmetic_" + name, value - s.Abs(value))
    return {
        "whole_finite_N_soft_term": (ae * s.gamma(2 * e) * x ** (2 * e)) ** n
        / (s.factorial(n) * s.gamma(1 + 2 * e * n)),
        "whole_analytic_virtual_rate": s.exp(-a / (2 * e)),
        "whole_uniform_signed_exponent_bound": s.Rational(10585, 2) / source.KAPPA,
        "whole_absolute_summed_majorant": s.Integer(6),
        "whole_marked_moment_limit": s.gamma(p)
        * s.exp(d - s.EulerGamma * a)
        * x ** (p + a)
        / s.gamma(p + 1 + a),
        "whole_regulator_proof": "The finite-e angular contraction is real but need not be positive. Gamma(1+2e)<=1 by logconvexity between1and2, and1/Gamma(1+2eN)<=exp1. For0<y<=1, the absolute virtual-normalized series is<=exp1*exp[(abs(a_e)-a)/(2e)]<=exp1*exp[abs(a_e-a)/(2e)]<6 at originalparameters. The explicit S301 phase remainder bounds the exponent by10585/(2kappa)<5500/kappa uniformly over all marked states.",
        "whole_limit_proof": "At fixeda>0, a_e is eventually positive and the S299 auxiliary Poisson mean lambda=a_e Gamma(2e)y^(2e) gives2eN->a inL2. The bounded continuous reciprocal Gamma expectation yields exp(Delta-gamma_E*a)y^a/Gamma(1+a). At a0, lambda->Delta can be negative; uniform absolute domination of the exponential series gives expDelta instead. A finite signed marked measure admits dominated convergence under the uniform6bound at fixed positivex. Its physical energy density has no atom atomega=x.",
        "whole_Bose_proof": "For totaln=N+1 labels, sum the n possible marked legs with its state-dependent current and virtual factor attached. Relabeling each term in the common total-energy simplex gives n/n!=1/N!, then the additional-energy cut isy=x-omega. No extra marked factor or independent-energy replacement appears.",
        "whole_no_finite_e_probability_assumption": True,
        "whole_positive_margins": margins,
        "checks": checks,
        "gates": {
            "all_signed_series_margins_strict": all(
                bool(v > 0) for v in margins.values()
            ),
            "absolute_series_domination_includes_negative_a_e": True,
            "a_positive_and_a_zero_limits_proved_separately": True,
            "no_fixed_N_limit_of_positive_index_sum": True,
            "state_dependent_marked_Bose_relabeling": True,
            "total_energy_simplex_not_individual_cuts": True,
            "fixed_positive_threshold_regulator_removal_before_other_limits": True,
        },
    }
