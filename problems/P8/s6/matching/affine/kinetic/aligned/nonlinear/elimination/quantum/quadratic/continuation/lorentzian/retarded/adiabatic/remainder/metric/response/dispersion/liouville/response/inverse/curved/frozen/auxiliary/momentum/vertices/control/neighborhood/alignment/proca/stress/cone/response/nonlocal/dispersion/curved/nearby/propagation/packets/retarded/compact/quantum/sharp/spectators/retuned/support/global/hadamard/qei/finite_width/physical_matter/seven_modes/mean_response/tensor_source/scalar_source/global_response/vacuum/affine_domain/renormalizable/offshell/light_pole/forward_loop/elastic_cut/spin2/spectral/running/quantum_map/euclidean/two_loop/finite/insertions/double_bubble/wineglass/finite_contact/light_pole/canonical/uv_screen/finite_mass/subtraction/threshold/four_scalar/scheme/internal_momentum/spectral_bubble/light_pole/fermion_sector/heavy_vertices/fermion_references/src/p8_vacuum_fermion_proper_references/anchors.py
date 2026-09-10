"""Finite scalar/gauge local anchors with the dimensional terms retained."""

from functools import cache

import sympy as sp


@cache
def data():
    x, r = sp.symbols("x r", positive=True)
    m, mu, b = sp.symbols("m mu boson_mass_squared", positive=True)
    eps, Jg0, Jg1 = sp.symbols("epsilon Jg0 Jg1", real=True)
    Delta = x * m * m + (1 - x) * b
    log = sp.log(Delta / mu**2)
    normalized = r + (1 - r) * x
    J0 = -1 - r * sp.log(r) / (1 - r)
    J1 = (-sp.Rational(3, 4) - (r - r * r / 2) * sp.log(r) + r - r * r / 4) / (
        1 - r
    ) ** 2
    R = 1 / (1 - r) + r * sp.log(r) / (1 - r) ** 2
    t = sp.Symbol("t", positive=True)
    primitive0 = t * sp.log(t) - t
    primitive1 = (t - t * t / 2) * sp.log(t) - t + t * t / 4
    primitiveR = x / (1 - r) - r * sp.log(normalized) / (1 - r) ** 2
    gauge_mass = sp.expand((4 - 2 * eps) * (1 / eps - Jg0) - 4 / eps).coeff(eps, 0)
    gauge_kinetic = sp.expand((2 - 2 * eps) * (1 / (2 * eps) - Jg1) - 1 / eps).coeff(
        eps, 0
    )
    gauge_log = sp.log(m * m / mu**2) - 1
    finite_mass = m * (-4 * gauge_log - 2)
    fixed_mu_vertex = sp.diff(finite_mass, m).subs(mu, m)
    checks = {
        "fixed_mu_mass_derivative_integrand": sp.simplify(
            sp.diff(m * log, m) - log - 2 * x * m * m / Delta
        ),
        "zeroth_log_primitive": sp.diff(primitive0, t) - sp.log(t),
        "first_weighted_log_primitive": sp.simplify(
            sp.diff(primitive1, t) - (1 - t) * sp.log(t)
        ),
        "scalar_J0_closed_integral": sp.simplify(
            (primitive0.subs(t, 1) - primitive0.subs(t, r)) / (1 - r) - J0
        ),
        "scalar_J1_closed_integral": sp.simplify(
            (primitive1.subs(t, 1) - primitive1.subs(t, r)) / (1 - r) ** 2 - J1
        ),
        "mass_derivative_ratio_primitive": sp.simplify(
            sp.diff(primitiveR, x) - x / normalized
        ),
        "mass_derivative_ratio_closed_integral": sp.simplify(
            primitiveR.subs(x, 1) - primitiveR.subs(x, 0) - R
        ),
        "normalized_denominator_above_x": normalized - x - r * (1 - x),
        "normalized_denominator_below_one": 1 - normalized - (1 - r) * (1 - x),
        "ratio_above_x": sp.factor(
            x / normalized - x - x * (1 - normalized) / normalized
        ),
        "ratio_below_one": sp.factor(1 - x / normalized - r * (1 - x) / normalized),
        "gauge_J0": sp.integrate(sp.log(x), (x, 0, 1)) + 1,
        "gauge_J1": sp.integrate((1 - x) * sp.log(x), (x, 0, 1)) + sp.Rational(3, 4),
        "gauge_R": sp.integrate(1, (x, 0, 1)) - 1,
        "scalar_kinetic_UV_half_weight": sp.integrate(1 - x, (x, 0, 1))
        - sp.Rational(1, 2),
        "gauge_mass_finite_epsilon_product": gauge_mass - (-4 * Jg0 - 2),
        "gauge_kinetic_finite_epsilon_product": gauge_kinetic - (-2 * Jg1 - 1),
        "gauge_mass_anchor_at_named_scale": gauge_mass.subs(Jg0, -1) - 2,
        "gauge_kinetic_anchor_at_named_scale": gauge_kinetic.subs(
            Jg1, -sp.Rational(3, 4)
        )
        - sp.Rational(1, 2),
        "gauge_Yukawa_vertex_fixed_mu_derivative": fixed_mu_vertex + 6,
        "resetting_mu_before_derivative_loses_eight": sp.diff(
            finite_mass.subs(mu, m), m
        )
        - fixed_mu_vertex
        - 8,
    }
    return {
        "parameter_Delta_at_zero_fermion_momentum": Delta,
        "scalar_ratio_r": "1/mF^2; gauge ratio zero",
        "J0_closed": J0,
        "J1_closed": J1,
        "R_closed": R,
        "J0_interval": [-1, 0],
        "J1_interval": [-sp.Rational(3, 4), 0],
        "R_interval": [sp.Rational(1, 2), 1],
        "scalar_local_Yukawa_integrand_at_fixed_mu": log + 2 * x * m * m / Delta,
        "gauge_dimensional_mass_finite_coefficient": gauge_mass,
        "gauge_dimensional_kinetic_finite_coefficient": gauge_kinetic,
        "gauge_finite_mass_ratio_over_aCf_div_Q": sp.Integer(2),
        "gauge_finite_kinetic_over_aCf_div_Q": sp.Rational(1, 2),
        "gauge_finite_Yukawa_ratio_over_aCf_div_Q": sp.Integer(-6),
        "scope": "Off-shell local reference at zero fermion/scalar momenta in Feynman gauge. Keep mu fixed when differentiating the background-dependent mass, then evaluate at mu=mF. Not a pole-mass or physical-spectrum assertion.",
        "checks": checks,
    }
