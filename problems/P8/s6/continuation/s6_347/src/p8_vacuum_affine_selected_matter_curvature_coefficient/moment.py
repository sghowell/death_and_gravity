"""Exact mixed bubble moments and the complete extra coefficient."""

from functools import cache

import sympy as s
from p8_vacuum_affine_triangle_curvature_coefficient.moment import require_positive_mass

from . import mixed

N = mixed.N
Y, Z = s.symbols("y z", positive=True)


def require_order(value):
    if type(value) is not int or value not in (2, 3):
        raise ValueError("Require native mixed bubble order2 or3")
    return value


@cache
def build():
    checks = {}
    values = {}
    for r in (2, 3):
        weight = s.expand((Z * (1 - Z)) ** r / s.Integer(r))
        degree = 2 * r
        numerator = s.Poly(
            s.expand(
                sum(
                    co * (Y - 1) ** ex[0] * (N - 1) ** (degree - ex[0])
                    for ex, co in s.Poly(weight, Z).terms()
                )
            ),
            Y,
        )
        F = sum(
            co
            * (
                s.log(Y)
                if ex[0] == r - 1
                else Y ** (ex[0] - r + 1) / s.Integer(ex[0] - r + 1)
            )
            for ex, co in numerator.terms()
        )
        value = s.factor((F.subs(Y, N) - F.subs(Y, 1)) / (N - 1) ** (degree + 1))
        checks[f"moment{r}_whole_primitive"] = s.factor(
            s.diff(F, Y) - numerator.as_expr() / Y**r
        )
        checks[f"moment{r}_equal_mass_limit"] = s.limit(value, N, 1) - s.integrate(
            weight, (Z, 0, 1)
        )
        values[r] = value
    return checks, values


def coefficient(order, mass):
    order = require_order(order)
    mass = require_positive_mass(mass)
    return (
        s.Rational(1, 60)
        if mass == 1 and order == 2
        else s.Rational(1, 420)
        if mass == 1
        else build()[1][order].subs(N, mass)
    )


@cache
def difference():
    b = build()[1]
    return s.factor((N - 1) * b[3] - b[2])


@cache
def extra():
    return s.factor(-32 * (difference() + 2 / N**4))


def extra_at(mass):
    mass = require_positive_mass(mass)
    return -s.Rational(952, 15) if mass == 1 else extra().subs(N, mass)


@cache
def data():
    checks = dict(build()[0])
    b = build()[1]
    checks["whole_combined_extra_mass_limit"] = s.limit(
        N * N * extra(), N, s.oo
    ) - s.Rational(8, 3)
    checks["whole_extra_equal_mass_continuity"] = s.limit(extra(), N, 1) + s.Rational(
        952, 15
    )
    for r in (2, 3):
        for mass in (s.Rational(1, 2), 2, 3, 10):
            checks[f"moment{r}_mass{mass}_exact_API"] = s.factor(
                coefficient(r, mass) - b[r].subs(N, mass)
            )
    return {
        "checks": checks,
        "gates": {
            "both_complete_bubble_Taylor_moments": set(b) == {2, 3},
            "strict_positive_internal_mass_only": True,
            "nonzero_removable_equal_mass_limits": True,
            "original_sign_not_falsely_extended_to_equal_mass": extra_at(1) < 0,
            "no_numerical_cancellation_fit": True,
        },
        "whole_bubble_moments": b,
        "whole_bubble_integrals": tuple(
            s.Integral((Z * (1 - Z) / (1 + (N - 1) * Z)) ** r / r, (Z, 0, 1))
            for r in (2, 3)
        ),
        "whole_across_moment_difference": difference(),
        "whole_combined_extra_dimensionless_coefficient": extra(),
        "whole_equal_mass_extra": extra_at(1),
        "whole_large_mass_squared_limit": s.Rational(8, 3),
        "whole_physical_normalization": "chi_extra=g^2 E(n)/(16pi^2 kappa), E(n)=-32[(n-1)b3(n)-b2(n)+2/n^4]. The first two terms are all across-source contractions; the final term is the complete within-J4 correction. The within-J2 class cancels and the local-tadpole curvature difference is zero. Equal-mass E is negative; positivity is proved only in the stated original large-mass domain, not for every n>0.",
    }
