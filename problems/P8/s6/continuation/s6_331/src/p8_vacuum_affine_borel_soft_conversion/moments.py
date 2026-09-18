"""Conditioned total-energy logarithms and the same-state soft insertion."""

from functools import cache

import sympy as s
from p8_vacuum_affine_leading_cloud_logarithmic_coefficient import cloud
from p8_vacuum_affine_radiative_logarithmic_soft_coefficient.measure import cutoff

from . import source


def parameters(index, resolution):
    a, x = cloud.index(index), cutoff(resolution)
    if a > 1:
        raise ValueError("Require0<=fixed-Born leading index<=1")
    return a, x


def weighted_log_moments(index, resolution):
    a, x = parameters(index, resolution)
    if a == 0:
        return s.S.Zero, s.S.Zero, s.S.Zero
    mean = a * x / (a + 1)
    return (
        mean,
        mean * (1 - s.log(x) + 1 / (a + 1)),
        mean * (-s.log(x) + s.polygamma(0, a + 2) + s.EulerGamma),
    )


def conditional_index_change_upper(index, resolution):
    mean, _, _ = weighted_log_moments(index, resolution)
    return 1400 * mean / source.KAPPA


def conditional_conversion_change_upper(index, resolution):
    _, log, _ = weighted_log_moments(index, resolution)
    return 10000 * log / source.KAPPA


def conditional_soft_transfer_tv_upper(index, resolution):
    _, log, remaining = weighted_log_moments(index, resolution)
    return (10000 * log + 1400 * remaining) / source.KAPPA


def conditional_regulator_dominator_upper(index, resolution):
    _, log, remaining = weighted_log_moments(index, resolution)
    return (15000 * log + 1400 * remaining) / source.KAPPA


def original_soft_transfer_tv_upper(resolution):
    x = cutoff(resolution)
    return 19000 * x * (1 - s.log(x)) / source.KAPPA**2


@cache
def data():
    a, x, p, y, e = s.symbols("a x p y e", positive=True)
    da, dD, dq = s.symbols("deltaa deltaDelta deltaq", real=True)
    moment = a * x**p / (a + p)
    mean = moment.subs(p, 1)
    log = mean * (1 - s.log(x) + 1 / (a + 1))
    b = s.Symbol("b", positive=True)
    beta_derivative = (s.polygamma(0, b) - s.polygamma(0, a + 1 + b)) * s.beta(a + 1, b)
    checks = {
        "weighted_log_from_moment_derivative": mean
        - s.diff(moment, p).subs(p, 1)
        - log,
        "remaining_energy_beta_derivative": s.expand_func(
            beta_derivative.subs(b, 1).replace(s.beta(a + 1, 1), 1 / (a + 1))
            + (s.polygamma(0, a + 2) + s.EulerGamma) / (a + 1)
        ),
        "same_state_soft_kernel_limit": s.limit(
            ((da + 2 * e * dD + e * e * dq) * y ** (2 * e) - da) / (2 * e), e, 0
        )
        - dD
        - da * s.log(y),
        "retained_remaining_energy_log": s.expand_log(
            s.log(x * y) - s.log(x) - s.log(y), force=True
        ),
        "zero_index_weighted_log": log.subs(a, 0),
        "harmonic_upper_endpoint": s.polygamma(0, 3) + s.EulerGamma - s.Rational(3, 2),
    }
    gates = {
        "conditional_transfer_dominator": 2 * 10000 + s.Rational(3, 2) * 1400 < 23000,
        "conditional_regulated_dominator": 2 * 15000 + s.Rational(3, 2) * 1400 < 33000,
        "physical_transfer_TV_roundup": s.Rational(4, 5) * 23000 < 19000,
        "same_radiative_state_real_and_virtual_pairing": True,
        "same_remaining_energy_for_Born_subtraction": True,
        "total_variation_regulator_limit_on_fixed_conditioned_reference": True,
        "boundary_R_equals_x_is_null_for_positive_index": True,
        "zero_index_handled_without_endpoint_singularity": True,
        "no_expansion_of_tiny_calorimetric_probability": True,
        "not_dimensionally_continued_outer_hard_amplitudes": True,
    }
    return {
        "checks": {name: s.simplify(value) for name, value in checks.items()},
        "gates": {name: bool(value) for name, value in gates.items()},
        "whole_exact_conditioned_log_moments": (
            mean,
            log,
            mean * (-s.log(x) + s.polygamma(0, a + 2) + s.EulerGamma),
        ),
        "whole_marked_regulator_limit": "For the fixed-Born D4 leading probability conditioned onR<=x, use Z_e(sigma;x-R)=[deltaa_e(sigma)*(x-R)^(2e)-deltaa(sigma)]/(2e). The measures converge in TV to[deltaDelta+deltaa*ln(x-R)]dP(.|cut); an integrable regulator-independent mean dominator is<33000a*x*(1+ln1/x)/kappa. The limiting TV is<23000a*x*(1+ln1/x)/kappa, physically<19000x*(1+ln1/x)/kappa^2.",
        "whole_same_prescription_boundary": "Only the additional soft factor is dimensionally continued. The full amplitude Coulomb phase is not erased. This named cloud insertion is not full hard real-virtual matching or a positive full detector rate.",
    }
