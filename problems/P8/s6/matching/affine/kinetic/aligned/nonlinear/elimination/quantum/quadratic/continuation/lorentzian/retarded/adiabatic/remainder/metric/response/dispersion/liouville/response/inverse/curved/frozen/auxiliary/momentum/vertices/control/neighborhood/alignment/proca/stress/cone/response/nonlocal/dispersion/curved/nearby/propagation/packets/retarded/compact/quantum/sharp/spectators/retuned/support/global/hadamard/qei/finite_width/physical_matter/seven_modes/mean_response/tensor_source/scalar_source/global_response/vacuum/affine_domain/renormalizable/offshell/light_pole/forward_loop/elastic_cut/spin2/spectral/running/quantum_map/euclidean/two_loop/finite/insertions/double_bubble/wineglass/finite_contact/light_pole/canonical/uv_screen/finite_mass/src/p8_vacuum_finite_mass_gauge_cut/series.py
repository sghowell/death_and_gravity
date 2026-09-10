"""All-degree box-tail bound with exact ordered-insertion multiplicities."""

from functools import cache
from itertools import permutations

import sympy as sp


def compositions(total, parts):
    if type(total) is not int or type(parts) is not int:
        raise TypeError("Require native integer composition arguments")
    if total < 0 or parts < 1:
        raise ValueError("Require nonnegative degree and positive part count")
    if parts == 1:
        return [(total,)]
    return [
        (j,) + r for j in range(total + 1) for r in compositions(total - j, parts - 1)
    ]


@cache
def data():
    n = sp.symbols("integer_total_insertions", integer=True, positive=True)
    x = sp.symbols("nonnegative_Neumann_ratio", nonnegative=True)
    t = sp.symbols("positive_radial_variable", positive=True)
    m = sp.symbols("positive_fermion_mass", positive=True)
    a, Y = sp.symbols("positive_gauge_squared positive_Yukawa_squared", positive=True)
    q = 16 * sp.pi**2
    orders = [("Phi1",) + p for p in permutations(("Phi2", "gauge1", "gauge2"))]
    weight = (n + 3) * (n + 1) / (6 * n)
    geometric = x**3 / (1 - x)
    majorant = sp.factor((x * sp.diff(geometric, x) + 5 * geometric) / 6)
    # Primitive of t/(1+t)^(2+n/2), normalized to zero at t=0.
    primitive = (
        4 / (n * (n + 2))
        - 2 * (1 + t) ** (-n / 2) / n
        + 2 * (1 + t) ** (-1 - n / 2) / (n + 2)
    )
    checks = {
        "six_cyclic_box_orders": len(orders) - 6,
        "weak_composition_ratio": sp.expand((n + 3) * (n + 2) * (n + 1) / 6)
        - sp.expand(sp.binomial(n + 3, 3).expand(func=True)),
        "weight_reduction": sp.factor(
            (n + 3) * (n + 2) * (n + 1) / (6 * n * (n + 2)) - weight
        ),
        "weight_majorant_positive_residual": sp.factor(
            (n + 5) / 6 - weight - (n - 3) / (6 * n)
        ),
        "geometric_tail_closed_form": sp.factor(
            majorant - x**3 * (8 - 7 * x) / (6 * (1 - x) ** 2)
        ),
        "primitive_derivative": sp.simplify(
            sp.diff(primitive, t) - t / (1 + t) ** (2 + n / 2)
        ),
        "primitive_zero_anchor": sp.simplify(primitive.subs(t, 0)),
        "primitive_infinity_anchor": sp.limit(primitive, t, sp.oo) - 4 / (n * (n + 2)),
        "tail_majorant_monotonic_derivative": sp.factor(
            sp.diff((8 - 7 * x) / (6 * (1 - x) ** 2), x)
            - (9 - 7 * x) / (6 * (1 - x) ** 3)
        ),
        "radial_integer_three_check": sp.integrate(
            t / (1 + t) ** sp.Rational(7, 2), (t, 0, sp.oo)
        )
        - sp.Rational(4, 15),
        "radial_integer_four_check": sp.integrate(t / (1 + t) ** 4, (t, 0, sp.oo))
        - sp.Rational(1, 6),
        "tail_constant_spin_orders_radial": 4 * 6 * 4 - 96,
        "tail_constant_majorant_and_routing": 96 * 6 * 12**3 - 995328,
        "strict_rounded_constant_margin": 1000000 - 995328 - 4672,
        "half_disc_six_majorant_margin": sp.Rational(6)
        - sp.Rational(16, 3)
        - sp.Rational(2, 3),
    }
    composition_rows = []
    for degree in range(13):
        values = compositions(degree, 4)
        composition_rows.append(
            {
                "degree": degree,
                "enumerated": len(values),
                "formula": sp.binomial(degree + 3, 3),
            }
        )
        checks[f"all_ordered_insertions_degree_{degree}"] = len(values) - sp.binomial(
            degree + 3, 3
        )
        checks[f"distinct_insertions_degree_{degree}"] = len(values) - len(set(values))
    return {
        "cyclic_labelled_box_orders": orders,
        "small_degree_independent_compositions": composition_rows,
        "all_degree_composition_count": sp.binomial(n + 3, 3),
        "all_degree_radial_integral": 4 / (q * n * (n + 2) * m**n),
        "dimensionless_radial_primitive": primitive,
        "weighted_tail_majorant": majorant,
        "uniform_per_color_polarization_box_remainder_upper": 1000000
        * a
        * Y
        / (q * m**3),
        "unrounded_box_remainder_upper": 995328 * a * Y / (q * m**3),
        "normalized_box_epsilon": 1000000 / m,
        "scope": "All Taylor degrees n>=3 are bounded, not just the enumerated examples. The complete regulated degree-zero term and the vanishing degree-one term are removed by Ward/Lorentz arguments; the finite degree-two term is the inherited leading operator. The bound requires mF>=24 and cumulative external Euclidean component norm strictly below 12.",
        "checks": checks,
    }
