"""Parameter-cut density and convergent moments of an unnormalizable measure."""

from functools import cache

import sympy as sp


@cache
def data():
    beta, x, u, T, C = sp.symbols("beta x u T C", positive=True)
    low, high = (1 - beta) / 2, (1 + beta) / 2
    A = x * (1 - x)
    weight = C * u * beta**3 / (u - 1) ** 2
    checks = {
        "lower_threshold_root": sp.expand((A - (1 - beta**2) / 4).subs(x, low)),
        "upper_threshold_root": sp.expand((A - (1 - beta**2) / 4).subs(x, high)),
        "allowed_parameter_interval_length": sp.expand(high - low - beta),
        "cut_parameter_interval_factorization": sp.expand(
            A - (1 - beta**2) / 4 - (x - low) * (high - x)
        ),
        "complete_density_beta_cubed": sp.factor((u - T) * beta - u * beta**3).subs(
            T, u * (1 - beta**2)
        ),
        "weight_positive_numerator": sp.factor(weight * (u - 1) ** 2 - C * u * beta**3),
        "weight_upper_factor_margin": sp.factor(
            4 * (u - 1) ** 2 - u**2 - (u - 2) * (3 * u - 2)
        ),
        "large_u_weight_lower_beta_test": sp.Rational(1, 2) ** 3
        - sp.Rational(1, 4) ** 2
        - sp.Rational(1, 16),
        "high_mass_gap_two_factor": sp.expand((u - 4) - u / 2 - (u - 8) / 2),
        "integrable_second_moment_tail": sp.integrate(u ** (-3), (u, T, sp.oo))
        - 1 / (2 * T * T),
        "integrable_first_moment_tail": sp.integrate(u ** (-2), (u, T, sp.oo)) - 1 / T,
    }
    return {
        "spectral_threshold": T,
        "parameter_roots": [low, high],
        "positive_spectral_weight": weight,
        "beta_dictionary": "beta=sqrt(1-T/u), u>=T=4mF^2.",
        "measure_upper": "0<=w(u)<=4C/u for u>=max(T,2), C=2NY/Q.",
        "total_measure_warning": "For u>=2T and T>=1, beta^3>=1/(2sqrt(2))>1/4 and w(u)>=C/(4u). Its total mass is infinite when C>0.",
        "convergence_scope": "The inserted propagator integral w(u)/(u-s) and the subtracted outer-bubble integrals converge. The total continuum weight without its propagator does not.",
        "checks": checks,
    }
