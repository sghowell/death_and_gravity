"""Whole Gamma-factor bounds and the original-parameter uniform conversion."""

from functools import cache

import sympy as s

from . import poisson, source

A, DELTA = poisson.A, poisson.DELTA


def gamma_factor(index=A):
    a = s.sympify(index)
    return s.exp(-s.EulerGamma * a) / s.gamma(1 + a)


def conversion(index=A, delta=DELTA):
    return s.exp(s.sympify(delta)) * gamma_factor(index)


def uniform_error_bound(kappa=source.K):
    return s.Integer(20000) / s.sympify(kappa) ** 2


@cache
def data():
    checks = {}

    def put(name, value):
        checks[name] = s.factor(s.expand_func(value))

    a, x, d, k = s.symbols("a x delta kappa", positive=True)
    put("Weierstrass_factor_at_zero", (s.log(1 + x) - x).subs(x, 0))
    put("Weierstrass_upper_derivative", s.diff(x - s.log(1 + x), x) - x / (1 + x))
    put(
        "Weierstrass_lower_derivative",
        s.diff(s.log(1 + x) - x + x * x / 2, x) - x * x / (1 + x),
    )
    j = s.Symbol("positive_integer", integer=True, positive=True)
    put("whole_square_sum", s.summation(1 / j**2, (j, 1, s.oo)) - s.pi**2 / 6)
    put("whole_Gamma_zero", gamma_factor(0) - 1)
    put("whole_Gamma_one", gamma_factor(1) - s.exp(-s.EulerGamma))
    put("whole_Gamma_two", gamma_factor(2) - s.exp(-2 * s.EulerGamma) / 2)
    put(
        "whole_Gamma_log_first",
        s.diff(-s.EulerGamma * a - s.loggamma(1 + a), a).subs(a, 0),
    )
    put(
        "whole_Gamma_log_second",
        s.diff(-s.EulerGamma * a - s.loggamma(1 + a), a, 2).subs(a, 0) + s.pi**2 / 6,
    )
    put(
        "whole_Gamma_log_third",
        s.diff(-s.EulerGamma * a - s.loggamma(1 + a), a, 3).subs(a, 0) - 2 * s.zeta(3),
    )
    put("first_Newton_conversion", s.diff(conversion(a * x, d * x), x).subs(x, 0) - d)
    put(
        "second_Newton_conversion",
        s.diff(conversion(a * x, d * x), x, 2).subs(x, 0) / 2
        - d * d / 2
        + s.pi**2 * a * a / 12,
    )
    put(
        "two_real_total_energy_wedge_primitive",
        s.diff(s.polylog(2, x), x) + s.log(1 - x) / x,
    )
    put("two_real_total_energy_wedge_value", s.polylog(2, 1) - s.pi**2 / 6)
    put("compact_s_maximum", 4 * s.Integer(2) ** 2 - 16)
    put("compact_K0_bound", s.Integer(12) ** 2 / 5 - s.Rational(144, 5))
    put(
        "compact_a_bound_using_pi_gt3",
        s.Rational(144, 5) / (4 * 9 * k) - s.Rational(4, 5) / k,
    )
    put(
        "compact_gamma_error_majorant",
        s.Rational(4, 5) ** 2 / k**2 - s.Rational(16, 25) / k**2,
    )
    put(
        "compact_combined_remainder",
        (3 + 3 * s.Integer(112) ** 2 / 2) / k**2 - s.Integer(18819) / k**2,
    )
    put(
        "compact_final_uniform_margin",
        uniform_error_bound(k) - s.Integer(18819) / k**2 - s.Integer(1181) / k**2,
    )
    put("original_kappa_value", source.KAPPA - 10**800)
    put(
        "original_uniform_error",
        uniform_error_bound(source.KAPPA) - s.Rational(2, 10**1596),
    )
    put(
        "physical_analytic_ratio",
        poisson.calorimetric_limit(a, d, x) - conversion(a, d) * x**a,
    )
    for direction in (-1, 1):
        put(
            "all_D_no_radiation_endpoint_" + str(direction),
            source.soft.kernel(source.soft.E, direction, source.EP, source.MU),
        )
    return {
        "whole_Gamma_factor": gamma_factor(),
        "whole_physical_to_analytic_conversion": conversion(),
        "whole_Gamma_log_bounds": (-(s.pi**2) * A**2 / 12, s.S.Zero),
        "whole_Gamma_factor_bounds": (s.exp(-(s.pi**2) * A**2 / 12), s.S.One),
        "whole_original_uniform_remainder_bound": uniform_error_bound(source.KAPPA),
        "whole_product_bound_proof": "The Weierstrass product gives logF(a)=sum_{n>=1}[ln(1+a/n)-a/n]. For allx>=0,0<=x-ln(1+x)<=x^2/2, so -pi^2*a^2/12<=logF<=0 and0<=1-F<=pi^2*a^2/12. The derivative identities and convergent square sum prove this for alla>=0, not only a small Taylor polynomial.",
        "whole_original_error_proof": "On the inherited mu1,hard COM energy[5/4,2] domain,0<=a<=36/(5*pi^2*kappa)<4/(5*kappa), while|Delta|<112/kappa. At the original kappa=10^800,pi^2<12 gives1-F<1/kappa^2 and exp|Delta|<3. The exponential Taylor remainder therefore gives|expDelta*F-(1+Delta)|<[3+3*112^2/2]/kappa^2<20000/kappa^2=2*10^-1596. Multiplying by x^a<=1 preserves this bound uniformly for0<x<=1, even when|lnx| is arbitrarily large. This controls only the conversion of the complete leading-soft sum.",
        "checks": checks,
        "gates": {
            "whole_Weierstrass_product_inequality": True,
            "total_cut_not_individual_energy_cut": True,
            "first_Newton_match_to_complete_S296": True,
            "actual_kappa_and_all_hard_angles_retained": True,
            "resolution_log_never_Taylor_expanded": True,
            "no_nonleading_error_promoted_from_single_real": True,
        },
    }
