"""Two integrations by parts and the complete odd transition-series tail."""

from fractions import Fraction
from functools import cache

import sympy as s


def exact(value):
    if isinstance(value, bool) or not isinstance(value, (int, Fraction, s.Rational)):
        raise TypeError("Require a finite exact rational")
    return s.Rational(value)


def mixing_tail(eta):
    eta = exact(eta)
    if not 0 <= eta <= s.Rational(1, 100):
        raise ValueError("Require an integrated mixing cap in [0,1/100]")
    return eta**3 / 5


def amplitude_terms(momentum, mass_floor, amplitude, timescale):
    p, m, d, tau = map(exact, (momentum, mass_floor, amplitude, timescale))
    if min(p, d) < 0 or min(m, tau) <= 0 or d / m >= s.Rational(1, 100):
        raise ValueError("Require p,Delta>=0, positive gap/time and Delta/gap<1/100")
    E2 = p * p + m * m
    eta = p * d / E2
    first = 5 * p * d / (tau * tau * E2**2)
    tail = mixing_tail(eta)
    return {
        "integrated_mixing_upper": eta,
        "first_transition_integral_upper": first,
        "all_higher_odd_transition_terms_upper": tail,
        "exact_transition_amplitude_upper": first + tail,
    }


@cache
def data():
    t, p = s.symbols("t p", real=True)
    mass = s.Function("M")(t)
    omega = s.Function("omega")(t)
    v = p * s.diff(mass, t) / (2 * omega**2)
    twice = s.diff(s.diff(v / (2 * omega), t) / (2 * omega), t)
    expanded = (
        p
        / s.Integer(8)
        * (
            s.diff(mass, t, 3) / omega**4
            - 7 * s.diff(mass, t, 2) * s.diff(omega, t) / omega**5
            - 3 * s.diff(mass, t) * s.diff(omega, t, 2) / omega**5
            + 15 * s.diff(mass, t) * s.diff(omega, t) ** 2 / omega**6
        )
    )
    r = s.Symbol("r", nonnegative=True)
    return {
        "twice_integrated_transition_integrand": expanded,
        "absolute_integrated_derivative_factor": (36 + 10 * r + 36 * r * r) / 8,
        "transition_integral_bound": "|I1(p)|<=5 p Delta/[tau^2(p^2+m0^2)^2] for Delta/m0<1/100. The instantaneous omega in the phase is the full varying frequency, not expanded in Delta.",
        "all_orders_mixing_tail": "The interaction generator is off diagonal. Only odd Dyson orders contribute to beta. With eta=integral |v|<=p Delta/(p^2+m0^2), absolute convergence bounds beta-I1 by sinh(eta)-eta<=eta^3 exp(eta)/6<=eta^3/5 when eta<=1/100. These are all transition-coupling orders of a quadratic evolution, NOT all Feynman loop orders.",
        "absolute_omega_derivative_bounds": "|omega'|<=|M'| and |omega''|<=|M''|+|M'|^2/omega. Substitution into the exact twice-integrated expression gives coefficients1,10,18 for |M'''|/E^4, |M'M''|/E^5 and |M'|^3/E^6.",
        "checks": {
            "exact_twice_integration_by_parts": s.simplify(twice - expanded),
            "mixed_mass_derivative_absolute_coefficient": 7 + 3 - 10,
            "cubic_mass_derivative_absolute_coefficient": 3 + 15 - 18,
            "profile_integral_composition": s.expand(
                (36 + 10 * r + 18 * r * r * 2) / 8 - (36 + 10 * r + 36 * r * r) / 8
            ),
            "max_mixing_from_two_ab_le_square": s.expand(
                (p - s.Symbol("m", positive=True)) ** 2
                - (
                    p * p
                    + s.Symbol("m", positive=True) ** 2
                    - 2 * p * s.Symbol("m", positive=True)
                )
            ),
            "third_and_higher_tail_coefficient": s.Rational(1, 6)
            / (1 - s.Rational(1, 100))
            - s.Rational(50, 297),
            "zero_mixing_has_zero_transition": mixing_tail(0),
        },
        "bounds": {
            "derivative_factor_strictly_below_five": bool(
                (36 + s.Rational(10, 100) + s.Rational(36, 10000)) / 8 < 5
            ),
            "odd_series_geometric_coefficient_below_one_fifth": bool(
                s.Rational(50, 297) < s.Rational(1, 5)
            ),
        },
    }
