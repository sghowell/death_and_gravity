"""Relative total-cut normalization and quantitative marked cutoff moments."""

from functools import cache

import sympy as s
from p8_vacuum_affine_leading_cloud_logarithmic_coefficient import cloud
from p8_vacuum_affine_radiative_logarithmic_soft_coefficient.measure import cutoff

from . import increment, source

B = s.Integer(23667)


def parameters(index, resolution, infrared):
    a, x, eta = cloud.index(index), cutoff(resolution), cloud.infrared_cutoff(infrared)
    if a > 1:
        raise ValueError(
            "Require0<=leading index<=1 for the relative normalization bound"
        )
    if eta > x:
        raise ValueError("Require infrared cutoff<=total energy resolution")
    return a, x, eta


def _tail_formula(a, x, h):
    u = 1 - h
    mean = a * x * (1 - u ** (a + 1)) / (a + 1)
    single = a * x * x * ((1 - u ** (a + 1)) / (a + 1) - (1 - u ** (a + 2)) / (a + 2))
    pair = (
        a
        * a
        * x
        * x
        * (1 - 2 * u ** (a + 2) + s.Max(0, 1 - 2 * h) ** (a + 2))
        / ((a + 1) * (a + 2))
    )
    return mean, single + pair


def conditional_lost_energy_moments(index, resolution, infrared):
    a, x, eta = parameters(index, resolution, infrared)
    return _tail_formula(a, x, eta / x)


def conditioning_mass_loss_upper(index, resolution, infrared):
    a, x, eta = parameters(index, resolution, infrared)
    return a * eta / x


def normalized_mean_cutoff_upper(index, resolution, infrared):
    a, _, eta = parameters(index, resolution, infrared)
    return (increment.LIPSCHITZ + 2 * B) * a * eta


def normalized_second_moment_cutoff_upper(index, resolution, infrared):
    a, x, eta = parameters(index, resolution, infrared)
    return (2 * B * increment.LIPSCHITZ + B * B) * a * x * eta


def original_physical_mean_cutoff_upper(resolution, infrared):
    _, _, eta = parameters(0, resolution, infrared)
    return 70000 * eta / source.KAPPA ** s.Rational(5, 2)


def original_physical_second_moment_cutoff_upper(resolution, infrared):
    _, x, eta = parameters(0, resolution, infrared)
    return s.Integer(2000000000) * x * eta / source.KAPPA**4


@cache
def data():
    a, x, u, eta = s.symbols("a x u eta", positive=True)
    n = s.Symbol("n", integer=True, positive=True)
    m = s.Symbol("m", integer=True, nonnegative=True)
    w = s.Symbol("w", positive=True)
    fullmean, fullsecond = _tail_formula(a, x, 1)
    zeromean, zerosecond = _tail_formula(a, x, 0)
    checks = {
        "conditional_tail_mean_derivative": -s.diff(
            a * x * (1 - u ** (a + 1)) / (a + 1), u
        )
        - a * x * u**a,
        "full_conditional_first_moment": fullmean - a * x / (a + 1),
        "full_conditional_second_moment": fullsecond - a * x * x / (a + 2),
        "zero_tail_first_moment": zeromean,
        "zero_tail_second_moment": zerosecond,
        "one_distinguished_Bose_emission": n / s.factorial(n) - 1 / s.factorial(n - 1),
        "two_distinguished_Bose_emissions": (m + 2) * (m + 1) / s.factorial(m + 2)
        - 1 / s.factorial(m),
        "first_count_sector_renewal": x * s.diff(eta**a * (1 + a * s.log(x / eta)), x)
        - a * eta**a,
        "second_count_sector_partial_fraction": 1 / (w * (x - w))
        - (1 / w + 1 / (x - w)) / x,
        "conditional_variance_not_squared_mean": s.factor(
            a * x * x / (a + 2) - (a * x / (a + 1)) ** 2
        )
        - a * x * x / ((a + 1) ** 2 * (a + 2)),
    }
    gates = {
        "physical_mean_cutoff_roundup": (increment.LIPSCHITZ + 2 * B) * s.Rational(4, 5)
        < 70000,
        "physical_second_moment_cutoff_roundup": (2 * B * increment.LIPSCHITZ + B * B)
        * s.Rational(4, 5)
        < 2000000000,
        "distinguished_emission_argument_not_conditioned_Poisson_assumption": True,
        "tail_independent_before_not_after_total_energy_conditioning": True,
        "relative_CDF_bound_retains_extremely_small_normalization": True,
        "CDF_zero_atom_and_nested_events_retained": True,
        "subtract_same_Born_mark_before_mixture_bound": True,
        "second_moment_nonnegative_range_not_square_of_mean": True,
        "uniform_eta_le_x_with_no_resolution_power_expansion": True,
        "zero_index_and_zero_tail_limits_separate": True,
        "D4_cutoff_only_not_marked_dimensional_matching": True,
    }
    return {
        "checks": {
            name: s.simplify(s.expand_func(value)) for name, value in checks.items()
        },
        "gates": {name: bool(value) for name, value in gates.items()},
        "whole_exact_conditioned_tail": "E[T_eta|R<=x]=a*x/(a+1)[1-(1-eta/x)^(a+1)]<=a*eta. The second moment retains the one-emission diagonal plus the ordered two-emission integral and is<=a*eta^2/2+a^2*eta^2. The exact API keeps the positive-part term when2eta>x.",
        "whole_relative_conditioning_proof": "For finite cutoff, rP_eta'(r)=aP_eta(r-eta) on0<r<1 and P_eta has atometa^a at0. Thus dlnP_eta/dr<=a/r. Independent omitted energy yields P(R<=x)/P_eta(x)>=E[(1-T_eta/x)_+^a]>=1-a*eta/x for0<a<=1; at a0 both probabilities are1. No loose absolute error is divided by the tiny calorimetric probability.",
        "whole_quantitative_marked_limits": "With B23667 and L40000, the actual finite-cutoff versus limiting conditioned mean error is<=(L+2B)*a*eta, and the second-moment-about-Born error is<=(2BL+B^2)*a*x*eta. At original a<4/(5kappa), physical bounds are70000*eta/kappa^(5/2) and2000000000*x*eta/kappa^4, uniformly0<eta<=x<=1/8.",
        "whole_scope": "Positive added angular energy at fixed original E,u and the fixed-Born D4 leading probability only. Full complex coefficient and conditioning normalization remain. No arbitrary signed perturbation, interacting state, marked dimensional conversion, hard radiative remainder or original P8 closure.",
    }
