"""Only repeated insertions of the computed one-loop kernel, with UV boundary."""

from functools import cache

import sympy as sp


def tails(alpha, order):
    if not isinstance(alpha, sp.Rational) or isinstance(alpha, sp.Float):
        raise TypeError("Require a SymPy exact rational multiplier bound")
    if not 0 < alpha < 1:
        raise ValueError("Require multiplier bound strictly between zero and one")
    if type(order) is not int or not 0 <= order <= 32:
        raise ValueError("Require native truncation order zero through 32")
    return {
        "order": order,
        "single_line_tail": alpha ** (order + 1) / (1 - alpha),
        "two_line_total_order_tail": alpha ** (order + 1)
        * (order + 2 - (order + 1) * alpha)
        / (1 - alpha) ** 2,
    }


@cache
def data():
    a, q, r = sp.symbols(
        "uniform_bound first_line_ratio second_line_ratio", positive=True
    )
    n = sp.Symbol("nonnegative_truncation_order", integer=True, nonnegative=True)
    z = sp.Symbol("formal_insertion_order")
    y, R = sp.symbols(
        "nonnegative_radial_momentum positive_radial_upper", positive=True
    )
    uv = sp.log(1 + R) + 1 / (1 + R) - 1
    checks = {}
    for order in range(5):
        single = sum(a**j for j in range(order + 1))
        two = sum((j + 1) * a**j for j in range(order + 1))
        checks["single_line_geometric_tail_" + str(order)] = sp.factor(
            1 / (1 - a) - single - a ** (order + 1) / (1 - a)
        )
        checks["two_line_total_order_tail_" + str(order)] = sp.factor(
            1 / (1 - a) ** 2
            - two
            - a ** (order + 1) * (order + 2 - (order + 1) * a) / (1 - a) ** 2
        )
        actual = (
            sp.series(1 / ((1 - z * q) * (1 - z * r)), z, 0, order + 1)
            .removeO()
            .coeff(z, order)
        )
        checks["two_distinct_line_order_" + str(order)] = sp.expand(
            actual - sum(q**j * r ** (order - j) for j in range(order + 1))
        )
    checks.update(
        {
            "positive_covariance_relative_shift": sp.factor(
                1 / (1 - a) - 1 - a / (1 - a)
            ),
            "two_line_relative_shift": sp.factor(
                1 / (1 - a) ** 2 - 1 - a * (2 - a) / (1 - a) ** 2
            ),
            "single_beyond_one_insertion_tail": sp.factor(
                1 / (1 - a) - 1 - a - a * a / (1 - a)
            ),
            "two_line_beyond_one_total_insertion_tail": sp.factor(
                1 / (1 - a) ** 2 - 1 - 2 * a - a * a * (3 - 2 * a) / (1 - a) ** 2
            ),
            "unsubtracted_four_dimensional_bubble_primitive": sp.factor(
                sp.diff(uv, R) - R / (1 + R) ** 2
            ),
            "unsubtracted_bubble_primitive_anchor": uv.subs(R, 0),
            "unsubtracted_bubble_logarithmic_asymptote": sp.limit(
                uv - sp.log(R), R, sp.oo
            )
            + 1,
        }
    )
    return {
        "uniform_bound": a,
        "single_line_covariance_factor": 1 / (1 - a),
        "two_line_covariance_factor": 1 / (1 - a) ** 2,
        "single_line_tail_formula": a ** (n + 1) / (1 - a),
        "two_line_total_order_tail_formula": a ** (n + 1)
        * (n + 2 - (n + 1) * a)
        / (1 - a) ** 2,
        "unsubtracted_bubble_positive_integral": sp.Integral(
            y / (1 + y) ** 2, (y, 0, R)
        ),
        "unsubtracted_bubble_closed_primitive": uv,
        "unsubtracted_bubble_diverges_to_positive_infinity": sp.limit(uv, R, sp.oo)
        == sp.oo,
        "scope": "Geometric iterations of one known one-loop self-energy only. Operator/wavepacket bounds or already absolutely integrable test kernels do not bound primitive higher-loop self-energies or renormalized scattering coefficients. The actual local quartic bubble is logarithmically divergent before subtraction.",
        "checks": checks,
    }
