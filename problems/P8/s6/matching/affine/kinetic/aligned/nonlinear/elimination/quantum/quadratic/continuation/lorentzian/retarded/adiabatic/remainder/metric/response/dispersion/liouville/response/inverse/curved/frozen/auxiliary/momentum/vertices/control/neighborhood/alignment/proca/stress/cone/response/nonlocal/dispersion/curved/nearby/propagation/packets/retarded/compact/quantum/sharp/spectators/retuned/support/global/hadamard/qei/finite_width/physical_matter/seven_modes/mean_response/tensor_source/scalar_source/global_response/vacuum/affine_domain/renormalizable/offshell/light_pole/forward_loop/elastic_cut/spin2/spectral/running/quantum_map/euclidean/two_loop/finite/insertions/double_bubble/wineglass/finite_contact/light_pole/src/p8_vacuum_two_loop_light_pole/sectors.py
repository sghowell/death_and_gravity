"""Every finite sunset Hepp sector and the twice-differentiated local sunset."""

from functools import cache
from itertools import permutations
from math import prod

import sympy as sp

from . import denominators, graphs


def ordered_bound(cotrees, edge_count, extra_powers=None):
    total = sp.S.Zero
    minimum = edge_count + 1
    patterns = set()
    orders = 0
    for order in permutations(range(edge_count)):
        ranks = {e: i for i, e in enumerate(order)}
        positions = sorted(tuple(sorted((ranks[i], ranks[j]))) for i, j in cotrees)
        chosen = positions[0]
        if not all(chosen[0] <= p[0] and chosen[1] <= p[1] for p in positions):
            raise ValueError("The selected co-tree does not dominate this ordering")
        if extra_powers is None:
            powers = tuple(-2 if i in chosen else 0 for i in range(edge_count))
        else:
            powers = tuple(
                extra_powers[order[i]] - (4 if i in chosen else 0)
                for i in range(edge_count)
            )
        tails = tuple(sum(p + 1 for p in powers[i:]) for i in range(edge_count))
        if min(tails) <= 0:
            raise ValueError("A UV sector is not integrable")
        minimum = min(minimum, min(tails))
        patterns.add((powers, tails))
        total += sp.Rational(1, prod(tails))
        orders += 1
    return {
        "orders": orders,
        "minimum_tail_exponent": minimum,
        "ordered_integral_patterns": tuple(sorted(patterns)),
        "unit_cube_integral_upper": total,
    }


@cache
def data():
    finite = []
    checks = {}
    patterns = set()
    for kind, choice in graphs.cases():
        if graphs.group(kind, choice) != "finite":
            continue
        g = graphs.refine(kind, choice)
        if graphs.uv(g):
            raise ValueError("A selected finite refinement has a UV core")
        d = denominators.point(kind, choice)
        row = ordered_bound(d["cotrees"], len(g["edges"]))
        row["choices"] = choice
        finite.append(row)
        patterns.update(row["ordered_integral_patterns"])
        key = "".join(map(str, choice))
        checks[key + "_all_five_edge_orders"] = row["orders"] - 120
        checks[key + "_positive_minimum_tail"] = row["minimum_tail_exponent"] - 1
        checks[key + "_exact_cube_upper"] = row["unit_cube_integral_upper"] - 14
    local = denominators.point("sunset", (0, 0))
    local_row = ordered_bound(local["cotrees"], 3, (2, 2, 2))
    patterns.update(local_row["ordered_integral_patterns"])
    checks["six_finite_refinements"] = len(finite) - 6
    checks["all_finite_edge_orders"] = sum(row["orders"] for row in finite) - 720
    checks["local_all_three_edge_orders"] = local_row["orders"] - 6
    checks["local_second_derivative_cube_upper"] = (
        local_row["unit_cube_integral_upper"] - 1
    )
    a, b, c = local["parameters"]
    checks["local_preserved_derivative_numerator"] = sp.expand(local["P"] - a * b * c)
    checks["local_stronger_denominator_identity"] = sp.expand(
        local["U"] * (a + b + c)
        - 9 * a * b * c
        - a * (b - c) ** 2
        - b * (c - a) ** 2
        - c * (a - b) ** 2
    )
    elementary = []
    for index, (powers, tails) in enumerate(sorted(patterns)):
        n = len(powers)
        t = sp.symbols("t0:" + str(n), positive=True)
        coordinates = [sp.prod(t[: i + 1]) for i in range(n)]
        jacobian = sp.Matrix(coordinates).jacobian(t).det()
        density = sp.prod(x**p for x, p in zip(coordinates, powers)) * jacobian
        expected = sp.prod(x ** (a - 1) for x, a in zip(t, tails))
        checks["independent_ordered_change_of_variables_" + str(index)] = sp.cancel(
            density - expected
        )
        direct = sp.prod(
            sp.integrate(x ** (a - 1), (x, 0, 1)) for x, a in zip(t, tails)
        )
        checks["independently_anchored_ordered_integral_" + str(index)] = (
            direct - sp.Rational(1, prod(tails))
        )
        elementary.append({"powers": powers, "tails": tails, "integral": direct})
    r = sp.Symbol("positive_box_radius", positive=True)
    checks["layer_cake_radial_integral"] = (
        sp.integrate(r * sp.exp(-r), (r, 0, sp.oo)) - 1
    )
    checks["finite_scaling_factor"] = 4 ** (5 - 4) - 4
    checks["local_differentiated_scaling_factor"] = 4 ** (3 + 6 - 8) - 4
    checks["finite_sunset_group_OS_prefactor"] = 6 * sp.Rational(1, 6) * 4 * 14 / 2 - 28
    checks["local_sunset_group_OS_prefactor"] = sp.Rational(1, 6) * 4 / 2 - sp.Rational(
        1, 3
    )
    return {
        "six_UV_finite_sunset_ordered_bounds": finite,
        "local_twice_differentiated_sunset_bound": local_row,
        "independent_ordered_change_of_variables": elementary,
        "finite_group_OS_majorant_without_g_squared_loop_measure": 28,
        "local_group_OS_majorant_without_L_squared_loop_measure": sp.Rational(1, 3),
        "scope": "All 720 finite-refinement edge orderings and all six local twice-differentiated orderings. The local numerator P^2 is retained to remove proper UV boundary singularities. These bounds integrate all Schwinger parameters and impose no physical cutoff.",
        "checks": checks,
    }
