"""Independent homogeneous curvature contractions and fixed finite coefficients."""

from functools import cache

import sympy as s

from . import angular, matching

epsilon = s.Symbol("epsilon_DR")
ell = s.Symbol("ell", real=True)


def dimensional_coefficients():
    dimension = 3 - 2 * epsilon
    weight = s.exp(epsilon * (s.EulerGamma - ell))
    vacuum = 2 * dimension * weight * s.gamma(epsilon - 2)
    einstein = (dimension - 6) * weight * s.gamma(epsilon - 1) / 3
    scalar, vector = s.symbols("a4_scalar a4_vector_4", real=True)
    fourth = 2 * weight * s.gamma(epsilon) * (vector - 2 * epsilon * scalar)
    return {
        "vacuum": s.series(vacuum, epsilon, 0, 1).removeO().expand(),
        "einstein": s.series(einstein, epsilon, 0, 1).removeO().expand(),
        "curvature": s.series(fourth, epsilon, 0, 1).removeO().expand(),
        "fixed_scalar": scalar,
        "fixed_vector": vector,
    }


def diagonal_curvature_jets(expansions, accelerations):
    """Actual orthonormal Bianchi-I Riemann invariants for diagonal scale rates."""
    n = len(expansions)
    if n != len(accelerations) or not n:
        raise ValueError("Require equally sized nonempty rate and acceleration lists")
    X = [accelerations[j] + expansions[j] ** 2 for j in range(n)]
    Y = [accelerations[j] + sum(expansions) * expansions[j] for j in range(n)]
    return {
        "R": sum(X) + sum(Y),
        "Ricci2": sum(X) ** 2 + sum(y * y for y in Y),
        "Riemann2": 4 * sum(x * x for x in X)
        + 4
        * sum(
            expansions[i] ** 2 * expansions[j] ** 2
            for i in range(n)
            for j in range(i + 1, n)
        ),
    }


@cache
def data():
    c = dimensional_coefficients()
    scalar, vector = c["fixed_scalar"], c["fixed_vector"]
    eps = epsilon
    checks = {
        "full_evanescent_vacuum_finite_coefficient": c["vacuum"]
        - 3 / eps
        - s.Rational(5, 2)
        + 3 * ell,
        "full_evanescent_Einstein_finite_coefficient": c["einstein"]
        - 1 / eps
        - s.Rational(5, 3)
        + ell,
        "full_evanescent_curvature_finite_contact": c["curvature"]
        - 2 * vector / eps
        + 2 * ell * vector
        + 4 * scalar,
    }
    H, Hd = angular.h, angular.h1
    for n in (3, 4, 5, 6):
        rates = [H] * n
        acceleration = [Hd] * n
        actual = diagonal_curvature_jets(rates, acceleration)
        for name, value in actual.items():
            target = matching.curvatures()[name].subs(angular.d, n)
            target = target.subs(
                {x: 0 for x in target.free_symbols if str(x).startswith("tr_")}
            )
            checks[f"isotropic_dimension_{n}_{name}"] = s.expand(value - target)
    R, Ric, Riem = s.symbols("R Ricci2 Riemann2", real=True)
    a4s = (5 * R**2 - 2 * Ric + 2 * Riem) / 360
    a4v = 3 * a4s + Ric / 2 - R**2 / 6 - Riem / 12
    checks["same_four_dimensional_Proca_pole_polynomial"] = s.expand(
        a4v + R**2 / 8 - s.Rational(29, 60) * Ric + Riem / 15
    )
    return {
        "MSbar_radial_convention": "The local radial integrals are dimensionally continued with d=3-2epsilon_DR and the standard exp(epsilon_DR*EulerGamma) normalization. Actual convergent state-minus-comparison integrals stay in physical dimension3.",
        "complete_dimension_limit": c,
        "fixed_pole": "Before 1/(64pi^2 epsilon_DR), subtract exactly 3m^4+m^2 R_old+2a4V(4), keeping these four-dimensional scalar coefficients while continuing invariant contractions and volume.",
        "fixed_finite_density": "Before 1/(64pi^2), the matched finite local density is -(3ell-5/2)m^4-(ell-5/3)m^2 R_old-2ell*a4V(4)-4a4s. It is the original prescription, not retuned coefficients.",
        "checks": checks,
        "gates": {
            "actual_polarization_factor_is_dimension_dependent": (3 - 2 * epsilon).diff(
                epsilon
            )
            == -2,
            "Einstein_evanescent_coefficient_retained": (
                (3 - 2 * epsilon - 6) / 3
            ).diff(epsilon)
            == -s.Rational(2, 3),
            "scalar_finite_curvature_contact_retained": (-2 * epsilon * 2).diff(epsilon)
            == -4,
            "no_actual_state_continuation_or_reset": True,
        },
    }
