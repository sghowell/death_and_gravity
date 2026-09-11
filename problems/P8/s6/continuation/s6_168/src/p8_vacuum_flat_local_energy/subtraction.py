"""Dimensional finite parts, with the full spatial angular measure retained."""

from functools import cache

import sympy as s


@cache
def data():
    eps, ell = s.symbols("epsilon log_mass_ratio")
    Q = 16 * s.pi**2
    common = s.exp((s.EulerGamma - ell) * eps) * s.gamma(1 + eps)
    z0 = 2 * common / ((eps - 1) * (eps - 2))
    z2 = (1 - 2 * eps / 3) * common
    M, rate = s.symbols("positive_mass real_mass_rate", positive=True)
    m, y, x, invariant = s.symbols("mean_mass Yukawa x invariant", positive=True)
    N = s.Symbol("N", positive=True)
    potential = -(M**4) * (s.log(M * M / (m * m)) - s.Rational(3, 2)) / Q
    kinetic_energy = -(rate**2) * (s.log(M * M / (m * m)) + s.Rational(2, 3)) / Q
    A = x * (1 - x)
    kernel_integrand = (4 * m * m - invariant) * (-s.log(1 - A * invariant / (m * m)))
    slope_zero = (
        2
        * N
        * y
        * y
        / Q
        * s.integrate(s.diff(kernel_integrand, invariant).subs(invariant, 0), (x, 0, 1))
    )
    return {
        "Dirac_trace_and_regularization": "Four-component Dirac trace, two occupied helicities, d=3-2epsilon spatial dimensions, full d-dimensional angular measure and MSbar subtraction at mu=mF. These are the same trace and scale conventions as the frozen fermion two-point and potential anchors.",
        "dimensionally_regularized_vacuum_energy": "-2 mu^(2epsilon) integral d^d p/(2pi)^d sqrt(p^2+M^2) = M^4/Q [Ibar+3/2-log(M^2/mu^2)]+O(epsilon)",
        "dimensionally_regularized_derivative_energy": "Mdot^2/4 mu^(2epsilon) integral d^d p/(2pi)^d p^2/(p^2+M^2)^(5/2) = Mdot^2/Q [Ibar-2/3-log(M^2/mu^2)]+O(epsilon)",
        "normalized_vacuum_gamma_factor": z0,
        "normalized_derivative_gamma_factor": z2,
        "MSbar_potential_per_color_flavor": potential,
        "MSbar_derivative_energy_per_color_flavor": kinetic_energy,
        "covariant_counterterm_dictionary": "The relevant pole action per Dirac copy is Ibar/Q integral sqrt(-g)[M^4-(partial M)^2], together with the required curvature counterterms. It cancels both displayed energy poles. M is a prescribed scalar under metric variation. The flat homogeneous 00 variation of a local curvature term R F(M) is zero; this is not true for all pressure components.",
        "finite_part_warning": "Keeping only the three-dimensional angular prefactor would lose the -2/3 evanescent contribution. No cutoff finite part or instantaneous-vacuum normal-ordering replaces the fixed MSbar prescription.",
        "relation_to_external_adiabatic_formula": "In del Rio et al. arXiv:1703.00908v2, section IV equations44-48, set a=1. The derivative term is p^2 sdot^2/(4(p^2+m^2)^(5/2)); the remaining terms are the fourth-order expansion in s of -2sqrt(p^2+(m+s)^2). Here the full M is kept in both terms, and the dimensionally regulated finite parts are restored. No derivative expansion approximates the exact state remainder.",
        "checks": {
            "vacuum_UV_pole_normalization": s.simplify(z0.subs(eps, 0) - 1),
            "derivative_UV_pole_normalization": s.simplify(z2.subs(eps, 0) - 1),
            "vacuum_MS_finite_part": s.simplify(
                s.diff(z0, eps).subs(eps, 0) - (s.Rational(3, 2) - ell)
            ),
            "derivative_MS_finite_part": s.simplify(
                s.diff(z2, eps).subs(eps, 0) - (-s.Rational(2, 3) - ell)
            ),
            "spatial_dimension_factor_not_omitted": s.diff(1 - 2 * eps / 3, eps)
            + s.Rational(2, 3),
            "same_frozen_zero_momentum_slope": s.simplify(
                slope_zero - 4 * N * y * y / (3 * Q)
            ),
            "same_MS_kinetic_energy_and_two_point_slope": s.simplify(
                N * kinetic_energy.subs({M: m, rate: y}) + slope_zero / 2
            ),
            "same_vacuum_potential_constant": s.simplify(
                potential.subs(M, m) - s.Rational(3, 2) * m**4 / Q
            ),
            "same_pair_quadratic_potential_coefficient": s.simplify(
                s.diff(potential, M, 2).subs(M, m) * N * y * y
                - 4 * N * y * y * m * m / Q
            ),
        },
    }
