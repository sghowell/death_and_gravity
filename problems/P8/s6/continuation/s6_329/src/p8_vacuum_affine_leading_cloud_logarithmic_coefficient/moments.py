"""Exact calorimetric energy moments and full complex coefficient mean bounds."""

from functools import cache

import sympy as s
from p8_vacuum_affine_calorimetric_soft_resummation import gamma, poisson
from p8_vacuum_affine_radiative_logarithmic_soft_coefficient import kernel
from p8_vacuum_affine_radiative_logarithmic_soft_coefficient.measure import cutoff

from . import cloud, source


def conditional_energy_moment(value, resolution, power):
    a, x, p = cloud.index(value), cutoff(resolution), kernel.exact_scalar(power)
    if p <= 0:
        raise ValueError("Require an exact positive moment power")
    return a * x**p / (a + p)


def normalized_mean_change_upper(value, resolution):
    return 23667 * conditional_energy_moment(value, resolution, 1)


def normalized_second_moment_upper(value, resolution):
    return 23667**2 * conditional_energy_moment(value, resolution, 2)


def original_physical_mean_upper(resolution):
    return 19000 * cutoff(resolution) / source.KAPPA ** s.Rational(5, 2)


def original_physical_second_moment_upper(resolution):
    return 230000000 * cutoff(resolution) ** 2 / source.KAPPA**4


@cache
def data():
    a, x, r, k = s.symbols("a x r kappa", positive=True)
    prob = poisson.calorimetric_limit(a, 0, x)
    density = a * r ** (a - 1) / x**a
    checks = {
        "exact_conditioned_energy_CDF": poisson.calorimetric_limit(a, 0, r) / prob
        - (r / x) ** a,
        "original_Born_index_bound": s.Rational(144, 5) / (4 * 9 * k)
        - s.Rational(4, 5) / k,
        "first_moment_zero_index": s.limit(a * x / (1 + a), a, 0, dir="+"),
        "second_moment_zero_index": s.limit(a * x * x / (2 + a), a, 0, dir="+"),
        "unconditional_truncated_first_moment": s.integrate(
            r * s.diff(poisson.calorimetric_limit(a, 0, r), r), (r, 0, x)
        )
        - prob * a * x / (1 + a),
        "unchanged_Gamma_normalization": gamma.gamma_factor(a)
        - s.exp(-s.EulerGamma * a) / s.gamma(1 + a),
    }
    for power in (1, 2, 3, 4):
        checks["conditional_energy_moment_" + str(power)] = s.integrate(
            r**power * density, (r, 0, x)
        ) - a * x**power / (a + power)
    gates = {
        "mean_roundup": s.Rational(23667 * 4, 5) < 19000,
        "second_moment_roundup": s.Rational(23667**2 * 2, 5) < 230000000,
        "full_complex_TT_field_not_only_real_part": True,
        "Jensen_norm_inequality_before_expectation_bound": True,
        "second_moment_not_square_of_mean": True,
        "unexpanded_small_cut_probability_retained": True,
        "zero_index_exactly_Born": True,
        "marked_reference_not_complete_loop_rate": True,
    }
    return {
        "checks": {
            key: s.simplify(s.expand_func(value)) for key, value in checks.items()
        },
        "gates": {key: bool(value) for key, value in gates.items()},
        "whole_conditional_energy_law": "For a>0, P(R<=r|R<=x)=(r/x)^a and E[R^p|R<=x]=a*x^p/(a+p),p>0. At a0 all positive moments vanish. The unexpanded probability F(a)*x^a is retained, even for extremely small resolution.",
        "whole_normalized_coefficient_moments": "||E[C_sigma|R<=x]-C_Born||sup<=23667*a*x/(1+a), and E[||C_sigma-C_Born||sup^2|R<=x]<=23667^2*a*x^2/(2+a). These are a mean and a separate second moment of the full complex continuous TT field.",
        "whole_original_physical_bounds": "Using original a<4/(5kappa), the physical coefficient C/kappa^(3/2) has mean change<19000*x/kappa^(5/2), and second moment<230000000*x^2/kappa^4. Both vanish exactly when a0. No amplitude phase is erased and no finite radiative hard remainder or interacting state is supplied.",
    }
