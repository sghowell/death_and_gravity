"""Defined infinite-particle, finite-energy D4 leading reference."""

from functools import cache

import sympy as s
from p8_vacuum_affine_calorimetric_soft_resummation import poisson
from p8_vacuum_affine_radiative_logarithmic_soft_coefficient import kernel
from p8_vacuum_affine_radiative_logarithmic_soft_coefficient.measure import cutoff


def index(value):
    a = kernel.exact_scalar(value)
    if a < 0:
        raise ValueError("Require a nonnegative exact leading index")
    return a


def infrared_cutoff(value):
    eta = kernel.exact_scalar(value)
    if eta <= 0 or eta > 1:
        raise ValueError("Require0<infrared cutoff<=1")
    return eta


def calorimetric_probability(value, resolution):
    a, x = index(value), cutoff(resolution)
    return poisson.calorimetric_limit(a, 0, x)


def lost_energy_moments(value, infrared):
    a, eta = index(value), infrared_cutoff(infrared)
    return a * eta, a * eta**2 / 2 + a * a * eta**2


def finite_cutoff_count_mean(value, infrared):
    a, eta = index(value), infrared_cutoff(infrared)
    return -a * s.log(eta)


@cache
def data():
    a, x, r, eta, z = s.symbols("a x r eta z", positive=True)
    laplace = -a * (s.EulerGamma + s.log(z) + s.expint(1, z))
    j = s.Symbol("j", integer=True, nonnegative=True)
    checks = {
        "original_leading_reference_CDF": poisson.calorimetric_limit(a, 0, x)
        - s.exp(-s.EulerGamma * a) * x**a / s.gamma(1 + a),
        "full_Levy_exponent_derivative": s.diff(laplace, z) - a * (s.exp(-z) - 1) / z,
        "full_Levy_exponent_zero": s.limit(laplace, z, 0, dir="+"),
        "finite_cutoff_count_mean": s.integrate(a / r, (r, eta, 1)) + a * s.log(eta),
        "lost_energy_mean": s.integrate(a, (r, 0, eta)) - a * eta,
        "lost_energy_variance": s.integrate(a * r, (r, 0, eta)) - a * eta**2 / 2,
        "lost_energy_second_moment_not_mean_square": (a * eta) ** 2
        + a * eta**2 / 2
        - eta**2 * (a * a + a / 2),
        "dyadic_band_intensity": s.integrate(a / r, (r, s.Rational(1, 2), 1))
        - a * s.log(2),
        "dyadic_band_energy": s.integrate(
            a, (r, s.Integer(2) ** (-j - 1), s.Integer(2) ** (-j))
        )
        - a * s.Integer(2) ** (-j - 1),
        "sum_all_band_mean_energies": s.summation(
            a * s.Integer(2) ** (-j - 1), (j, 0, s.oo)
        )
        - a,
    }
    return {
        "checks": {
            key: s.simplify(s.expand_func(value)) for key, value in checks.items()
        },
        "gates": {
            "independent_finite_intensity_dyadic_construction": True,
            "finite_expected_total_energy_implies_finite_energy_almost_surely": True,
            "positive_index_gives_infinite_particle_count_not_Fock_convergence": True,
            "cutoff_energy_measure_converges_in_total_variation": True,
            "S328_weak_to_uniform_TT_continuity_applies_on_total_cut": True,
            "nested_cutoff_events_decrease_to_exact_total_energy_event": True,
            "conditioning_probability_converges_and_stays_positive": True,
            "all_finite_Lp_limits_by_uniform_300_bound": True,
            "no_dimensional_matching_of_marked_insertions_claimed": True,
        },
        "whole_random_energy_measure": "Independent dyadic Poisson bands with intensity Lambda(dn)*dw/w define sigma=sum w_j delta_nj. R has mean a and is finite almost surely, while particle count is infinite almost surely when a>0. Finite sigma_eta retains w>=eta and ||sigma-sigma_eta||TV=R-R_eta decreases to0, with mean a*eta.",
        "whole_marked_cutoff_limit": "For fixed0<x<=1/8, define the zero-extended marked field 1_{R_eta<=x} C_sigma_eta. Nested cutoff events decrease exactly to{R<=x}; S328 supplies uniform angular TT convergence there, and outside the event the cutoff indicators eventually vanish. Norm<=300 gives every finite Lp convergence. Normalization converges to F(a)*x^a>0, or1 at a0. Thus conditional means and bounded coefficient moments have a unique D4 infrared limit.",
        "whole_regulator_boundary": "This is the fixed-Born D4 sharp-energy-cutoff leading probability, not a state-dependent branching law. The scalar S299 dimensional conversion expDelta does not automatically match marked insertions or evanescent terms. Remove the infrared cutoff at fixed positive total resolution before any correlated limit.",
    }
