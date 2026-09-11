"""Full-resolvent common-class action comparison of two polynomial maps."""

from fractions import Fraction
from functools import cache

import sympy as s

from . import gate


def exact(v):
    if isinstance(v, bool) or not isinstance(v, (int, Fraction, s.Rational)):
        raise TypeError("Require a finite exact rational")
    return s.Rational(v)


def gate_difference_cap(Xcap):
    x = exact(Xcap)
    if x < 0:
        raise ValueError("Require a nonnegative gradient-invariant cap")
    return sum(
        abs(c) * x**j for j, c in gate.data()["step_nonzero_coefficients"].items()
    )


def action_difference(C, delta, L, g, M):
    C, delta, L, g, M = map(exact, (C, delta, L, g, M))
    if min(C, delta, L, g) < 0 or M <= 66**2:
        raise ValueError("Require nonnegative interaction/norm caps and M>66^2")
    A = 1 + C
    Jdelta = delta * (2 * A + delta)
    free = 100 * A * delta + 545 * delta**2
    local = L * delta * (2 * A + delta) * (A * A + (A + delta) ** 2) / 24
    heavy = g * (2 * A * A * Jdelta + Jdelta**2) / (8 * (M - 66**2))
    return {
        "old_field_norm_factor": A,
        "new_field_norm_factor": A + delta,
        "quadratic_source_difference_norm_factor": Jdelta,
        "free_action_difference_coefficient": free,
        "local_quartic_difference_coefficient": local,
        "full_heavy_resolvent_difference_coefficient": heavy,
        "complete_action_difference_coefficient": free + local + heavy,
    }


@cache
def data():
    A, d = s.symbols("A delta", positive=True)
    X = s.Symbol("X", nonnegative=True)
    cap = sum(
        abs(c) * X**j for j, c in gate.data()["step_nonzero_coefficients"].items()
    )
    return {
        "gate_difference_coefficient_majorant": cap,
        "common_class": "Same real Schwartz Psi with Fourier support in the Euclidean four-momentum unit ball and Fourier L1 norm <=1. Let U=||Psi||_2. All jets have sup norm <=1 and L2 norm <=U. The canonical target gradient invariant obeys |X|<=4/kappa.",
        "map_norm_argument": "R_old has norms <=C_R and C_R U. Every monomial of (q8-1)R is a finite product of Psi jets; place one factor in L2 and the rest in L-infinity. Thus its sup/L2 norms are bounded by delta=qcap*C_R and delta U, with Fourier radius at most33.",
        "full_resolvent_argument": "Both squared sources have Fourier radius at most66, so |p_Minkowski^2|<=4356. For M>4356 the complete heavy inverse has operator norm <=1/(M-4356) on their joint support. No heavy kernel is expanded. Old source radius6 is retained in its earlier estimate; only the map difference uses the enlarged support.",
        "checks": {
            "free_cross_derivative_and_mass_factor": 3 * 33 + 1 - 100,
            "free_difference_square_factor": s.Rational(33**2 + 1, 2) - 545,
            "quadratic_source_difference_factor": s.expand(
                (A + d) ** 2 - A * A - d * (2 * A + d)
            ),
            "quartic_L2_Linfinity_factor": s.expand(
                (A + d) ** 4 - A**4 - d * (2 * A + d) * (A * A + (A + d) ** 2)
            ),
            "heavy_bilinear_difference_factor": s.expand(
                (A * A + d * (2 * A + d)) ** 2
                - A**4
                - 2 * A * A * d * (2 * A + d)
                - (d * (2 * A + d)) ** 2
            ),
            "first_gate_cap_degree": min(
                j for j in gate.data()["step_nonzero_coefficients"]
            )
            - 8,
            "free_linearized_corner_factor": 100 * s.Rational(5, 4)
            + 545 * s.Rational(1, 4)
            - s.Rational(1045, 4),
            "local_linearized_corner_factor": s.Rational(11, 4)
            * (s.Rational(25, 16) + s.Rational(9, 4))
            / 24
            - s.Rational(671, 1536),
            "heavy_linearized_corner_factor": (
                2 * s.Rational(25, 16) * s.Rational(11, 4)
                + s.Rational(1, 4) * s.Rational(11, 4) ** 2
            )
            / 8
            - s.Rational(671, 512),
        },
    }


def linear_enclosure(C, delta, L, g, M):
    C, delta, L, g, M = map(exact, (C, delta, L, g, M))
    if min(C, delta, L, g) < 0 or max(C, delta) > s.Rational(1, 4) or M <= 66**2:
        raise ValueError(
            "Require C and delta in [0,1/4], nonnegative couplings and M>66^2"
        )
    free = 300 * delta
    local = L * delta
    heavy = 2 * g * delta / (M - 66**2)
    return {
        "free_action_difference_coefficient": free,
        "local_quartic_difference_coefficient": local,
        "full_heavy_resolvent_difference_coefficient": heavy,
        "complete_action_difference_coefficient": free + local + heavy,
    }
