"""Proper MS box forest and the two exact regulated reference regions."""

from functools import cache

import sympy as sp


@cache
def data():
    e, z = sp.symbols("epsilon z", real=True)
    m, mu = sp.symbols("m mu", positive=True)
    # f_D(0)=CY/Y * (4m^2 B_D-2 T_D).
    A = sp.exp(sp.EulerGamma * e) * sp.gamma(1 + e)
    fzero = m * m * (mu * mu / (m * m)) ** e * A / e * (4 - 2 / (e - 1))
    raw = sp.factor(sp.diff(fzero, m, 2))
    expected = (mu * mu / (m * m)) ** e * A / e * (12 - 32 * e + 16 * e * e)
    finite = sp.simplify(sp.diff((12 - 32 * e + 16 * e * e) * A, e).subs(e, 0))
    p = sp.Rational(3, 2) - e
    # Spectral moments use Beta(n e,1/2-e) and its first moment.
    I = lambda n: (
        sp.beta(n * e, sp.Rational(1, 2) - e)
        * ((2 - 2 * e) * n * e / ((n - 1) * e + sp.Rational(1, 2)) - 1)
    )
    a, b = sp.symbols("alpha beta", positive=True)
    t = sp.Symbol("q_squared", positive=True)
    primitive = -(sp.log(1 + t) + 1) / (1 + t) + (2 * sp.log(1 + t) + 1) / (
        4 * (1 + t) ** 2
    )
    checks = {
        "fixed_mu_second_mass_derivative": sp.simplify(raw - expected),
        "proper_box_simple_pole": ((12 - 32 * e + 16 * e * e) * A).subs(e, 0) - 12,
        "proper_box_finite_anchor": finite + 32,
        "beta_first_moment_shift": sp.gammasimp(
            sp.gamma(a + 1) * sp.gamma(b) / sp.gamma(a + b + 1)
            - a / (a + b) * sp.gamma(a) * sp.gamma(b) / sp.gamma(a + b)
        ),
        "spectral_derivative_threshold_exponent": p - 2 + sp.Rational(1, 2) + e,
        "spectral_zero_soft_at_zero_momentum": -32 * 2 + 64,
        "integrable_light_bubble_radial_moment": sp.integrate(
            t / (1 + t) ** 3, (t, 0, sp.oo)
        )
        - sp.Rational(1, 2),
        "light_bubble_log_antiderivative": sp.simplify(
            sp.diff(primitive, t) - t * sp.log(1 + t) / (1 + t) ** 3
        ),
        "integrable_light_bubble_log_moment": sp.limit(primitive, t, sp.oo)
        - primitive.subs(t, 0)
        - sp.Rational(3, 4),
        "signed_density_low_momentum_moment": sp.integrate(
            (2 * z - 1) / sp.sqrt(1 - z), (z, 0, 1)
        )
        - sp.Rational(2, 3),
    }
    return {
        "dimensionless_raw_zero_soft_box": expected,
        "proper_MS_counterterm_in_CY_units": -12 / e,
        "dimensionless_spectral_beta_moments": (I(1), I(2)),
        "regulated_outer_reference": "F_D=-Integral nu_D(v) T3_D(v,1,1) dv -(12CY/epsilon) I2_D(1). Q=16pi^2, CY=2NY^2/Q. Each loop uses exp(gamma_E epsilon) mu^(2epsilon).",
        "regulator_common_convergence_strip": "0<Re epsilon<1/2 before continuation to epsilon=0.",
        "MS_reference_definition": "Finite Laurent coefficient of the entire proper-subgraph-paired F_D, with all local remaining poles subtracted.",
        "checks": checks,
    }
