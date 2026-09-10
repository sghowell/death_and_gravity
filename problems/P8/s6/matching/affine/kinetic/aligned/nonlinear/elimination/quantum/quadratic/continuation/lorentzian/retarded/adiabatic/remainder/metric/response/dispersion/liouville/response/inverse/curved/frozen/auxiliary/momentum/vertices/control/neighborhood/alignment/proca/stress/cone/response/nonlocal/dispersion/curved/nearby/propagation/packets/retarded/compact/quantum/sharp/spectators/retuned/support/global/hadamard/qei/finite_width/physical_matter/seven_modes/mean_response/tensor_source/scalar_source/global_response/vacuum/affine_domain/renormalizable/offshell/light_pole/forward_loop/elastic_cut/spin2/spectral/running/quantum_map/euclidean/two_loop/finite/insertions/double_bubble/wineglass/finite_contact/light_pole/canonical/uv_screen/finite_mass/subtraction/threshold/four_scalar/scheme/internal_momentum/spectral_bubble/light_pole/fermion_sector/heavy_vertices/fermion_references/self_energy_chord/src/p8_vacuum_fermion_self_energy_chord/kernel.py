"""Whole one-loop MS proper fermion kernels, not only their local anchors."""

from functools import cache

import sympy as sp
from p8_vacuum_fermion_proper_references import conversion


@cache
def data():
    prior = conversion.data()
    Y, a, Cf, Q = sp.symbols("Y a Cf Q", positive=True)
    m, x = sp.symbols("m x", positive=True)
    p2 = sp.Symbol("Euclidean_p_squared")
    eps = sp.Symbol("epsilon")
    As0, As1, Ag0, Ag1 = sp.symbols(
        "scalar_log0 scalar_log1 gauge_log0 gauge_log1", real=True
    )
    Delta_s = x * m * m + (1 - x) + x * (1 - x) * p2
    Delta_g = x * m * m + x * (1 - x) * p2
    Js = sp.log(Delta_s / (m * m))
    Jg = sp.log(Delta_g / (m * m))
    scalar_mass = Y * As0 / Q
    scalar_kin = -Y * As1 / Q
    gauge_mass = a * Cf * (-4 * Ag0 - 2) / Q
    gauge_kin = a * Cf * (-2 * Ag1 - 1) / Q
    reg_mass = (-Y * (1 / eps - As0) + a * Cf * (4 - 2 * eps) * (1 / eps - Ag0)) / Q
    reg_kin = (
        Y * (1 / (2 * eps) - As1) + a * Cf * (2 - 2 * eps) * (1 / (2 * eps) - Ag1)
    ) / Q
    paired_mass = reg_mass + prior["fermion_total_mass_UV_counterterm_over_m"] / eps
    paired_kin = reg_kin + prior["fermion_total_kinetic_UV_counterterm"] / eps
    J0, J1 = sp.symbols("J0 J1", real=True)
    at_zero = {As0: J0, As1: J1, Ag0: -1, Ag1: -sp.Rational(3, 4)}
    return {
        "scalar_parameter_Delta": Delta_s,
        "gauge_parameter_Delta": Delta_g,
        "MS_scalar_mass_inverse_correction": m * Y * sp.Integral(Js, (x, 0, 1)) / Q,
        "MS_scalar_kinetic_inverse_coefficient": -Y
        * sp.Integral((1 - x) * Js, (x, 0, 1))
        / Q,
        "MS_gauge_mass_inverse_correction": a
        * Cf
        * m
        * (-2 - 4 * sp.Integral(Jg, (x, 0, 1)))
        / Q,
        "MS_gauge_kinetic_inverse_coefficient": a
        * Cf
        * (-1 - 2 * sp.Integral((1 - x) * Jg, (x, 0, 1)))
        / Q,
        "mass_coefficient_convention": "Coefficient of identity in the Euclidean inverse; kinetic coefficient multiplies i slash(p)",
        "subgraph_ownership": "Complete one-loop proper MS self-energy, including its scalar/gauge mass and kinetic counterterms exactly once. The outer marked propagator is -S Sigma_MS S.",
        "scope": "All internal momenta are retained. These are formal one-loop proper kernels, not exact charged-particle propagators or a MOM-only Taylor replacement.",
        "checks": {
            "same_scalar_zero_momentum_parameter_mass": sp.expand(
                Delta_s.subs(p2, 0) - x * m * m - (1 - x)
            ),
            "same_gauge_zero_momentum_parameter_mass": sp.expand(
                Delta_g.subs(p2, 0) - x * m * m
            ),
            "proper_mass_pole_cancels": sp.expand(paired_mass).coeff(eps, -1),
            "proper_kinetic_pole_cancels": sp.expand(paired_kin).coeff(eps, -1),
            "full_MS_mass_finite_kernel": sp.expand(paired_mass).coeff(eps, 0)
            - scalar_mass
            - gauge_mass,
            "full_MS_kinetic_finite_kernel": sp.expand(paired_kin).coeff(eps, 0)
            - scalar_kin
            - gauge_kin,
            "same_frozen_finite_mass_anchor": sp.factor(
                (scalar_mass + gauge_mass).subs(at_zero)
                - prior["finite_mass_relative_eta"]
            ),
            "same_frozen_finite_kinetic_anchor": sp.factor(
                (scalar_kin + gauge_kin).subs(at_zero)
                - prior["finite_fermion_kinetic_z"]
            ),
            "gauge_mass_epsilon_constant_not_omitted": gauge_mass
            - (-4 * a * Cf * Ag0 / Q)
            + 2 * a * Cf / Q,
            "gauge_kinetic_epsilon_constant_not_omitted": gauge_kin
            - (-2 * a * Cf * Ag1 / Q)
            + a * Cf / Q,
        },
    }
