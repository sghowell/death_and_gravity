"""Pressure finite parts with the d-dimensional isotropic stress factor."""

from functools import cache

import sympy as s


@cache
def data():
    e, ell = s.symbols("epsilon log_mass_ratio")
    dim = 3 - 2 * e
    common = s.exp((s.EulerGamma - ell) * e) * s.gamma(1 + e)
    J0 = (1 - 2 * e / 3) * common
    J1 = (1 - 2 * e / 3) * (1 - 2 * e / 5) * common
    MM = -2 * J0 / dim
    vv = J1 / dim
    M, v, a, Q = s.symbols("M Mdot Mddot Q", positive=True)
    pressure2 = (2 * M * a * ell / 3 - v * v * ell / 3 + 2 * v * v / 3) / Q
    rho2 = -v * v * (ell + s.Rational(2, 3)) / Q
    x = s.Symbol("x", positive=True)
    radial0 = s.integrate(x**3 / (1 + x * x) ** 3, (x, 0, s.oo))
    radial1 = s.integrate(x**5 / (1 + x * x) ** 5, (x, 0, s.oo))
    D = lambda F: s.diff(F, M) * v + s.diff(F, v) * a
    mu = s.Symbol("mu", positive=True)
    improvement = M * M * (s.log(M * M / (mu * mu)) - 1)
    diff = (pressure2 - rho2).subs(ell, s.log(M * M / (mu * mu)))
    # At a=1, equation57's two derivative terms agree with the full-M formula.
    paper_rate = (
        -5 * M**4 / (12 * (x * x + M * M) ** s.Rational(7, 2))
        + M * M / (3 * (x * x + M * M) ** s.Rational(5, 2))
        + 1 / (12 * (x * x + M * M) ** s.Rational(3, 2))
    )
    our_rate = x * x * (x * x + 6 * M * M) / (12 * (x * x + M * M) ** s.Rational(7, 2))
    paper_acc = M**3 / (6 * (x * x + M * M) ** s.Rational(5, 2)) - M / (
        6 * (x * x + M * M) ** s.Rational(3, 2)
    )
    return {
        "dimensional_spatial_pressure_factor": 1 / dim,
        "normalized_J0_pole_factor": J0,
        "normalized_J1_pole_factor": J1,
        "J0_integral": "integral d^dp/(2pi)^d p^2/omega^5=4/Q[Ibar-log(M^2/mu^2)-2/3]+O(epsilon)",
        "J1_integral": "integral d^dp/(2pi)^d p^4/omega^7=4/Q[Ibar-log(M^2/mu^2)-16/15]+O(epsilon)",
        "J2_integral_at_four_dimensions": "integral d^3p/(2pi)^3 p^2/omega^7=8/(5Q M^2)",
        "full_dimensional_second_pressure": "-M Mddot/(2d) J0+Mdot^2/(4d) J1+3M^2 Mdot^2/(2d) J2, d=3-2epsilon. Retain the epsilon dependence of 1/d before the finite part.",
        "restored_MSbar_second_pressure": pressure2,
        "restored_MSbar_zeroth_pressure": "-V_MS(M), since dimensional integration by parts gives integral p^2/omega=-d integral omega.",
        "curvature_improvement_identity": "p2_MS-rho2_MS=partial_t^2{M^2[log(M^2/mu^2)-1]}/(3Q). This is a nonzero local improvement in homogeneous flat pressure even though its 00 component vanishes. It must not be omitted just because the background curvature itself is zero.",
        "renormalization_prescription": "Same four-component Dirac dimensional trace, MSbar scale and zero additional finite curvature coupling as the minimal extension of S6.168. Curvature pole counterterms are retained. The displayed flat tensor difference can equivalently be written -(eta_mu_nu Box-partial_mu partial_nu)F/(3Q), F=M^2(log(M^2/mu^2)-1); this fixes the tensor convention without guessing a Riemann-sign convention for an action coefficient.",
        "source_crosscheck": "The Mddot coefficient agrees with del Rio et al. arXiv:1703.00908v2 equation56 at a=1; the Mdot^2 coefficient agrees with equation57. Those are checks of the full-M UV terms expanded to their source order, not an approximation to the exact state.",
        "checks": {
            "J0_pole_normalization": s.simplify(J0.subs(e, 0) - 1),
            "J1_pole_normalization": s.simplify(J1.subs(e, 0) - 1),
            "J0_finite_part": s.simplify(
                s.diff(J0, e).subs(e, 0) + ell + s.Rational(2, 3)
            ),
            "J1_finite_part": s.simplify(
                s.diff(J1, e).subs(e, 0) + ell + s.Rational(16, 15)
            ),
            "pressure_mass_acceleration_pole": s.simplify(
                MM.subs(e, 0) + s.Rational(2, 3)
            ),
            "pressure_mass_acceleration_finite": s.simplify(
                s.diff(MM, e).subs(e, 0) - 2 * ell / 3
            ),
            "pressure_rate_squared_pole": s.simplify(vv.subs(e, 0) - s.Rational(1, 3)),
            "pressure_rate_squared_finite_including_J2": s.simplify(
                s.diff(vv, e).subs(e, 0)
                + s.Rational(4, 5)
                - (s.Rational(2, 3) - ell / 3)
            ),
            "epsilon_isotropic_factor_not_dropped": s.diff(1 / dim, e).subs(e, 0)
            - s.Rational(2, 9),
            "exact_local_improvement_difference": s.simplify(
                diff - D(D(improvement)) / (3 * Q)
            ),
            "complete_linear_x_radial": radial0 - s.Rational(1, 4),
            "complete_cubic_x_radial": radial1 - s.Rational(1, 24),
            "full_linear_x_factor": s.Rational(1, 2)
            * 2
            * s.Rational(2, 3)
            * 2
            * radial0
            - s.Rational(1, 3),
            "full_cubic_x_factor": s.Rational(1, 2)
            * 2
            * s.Rational(2, 3)
            * 3
            * 8
            * radial1
            - s.Rational(2, 3),
            "external_pressure_acceleration_coefficient": s.simplify(
                paper_acc + x * x * M / (6 * (x * x + M * M) ** s.Rational(5, 2))
            ),
            "external_pressure_rate_coefficient": s.simplify(paper_rate - our_rate),
        },
    }
