"""Third-order map identity and independent omitted-interaction controls."""

from functools import cache

import sympy as s


@cache
def data():
    x, r = s.symbols("x r")
    h, K = s.symbols("h K", positive=True)
    F = x + r * x**3
    jac = 1 + 3 * r * x * x
    weight = (
        1
        - r * K * x**4 / h
        + r**2 * (K**2 * x**8 / (2 * h**2) - K * x**6 / (2 * h))
        + r**3 * (K**2 * x**10 / (2 * h**2) - K**3 * x**12 / (6 * h**3))
    )

    def average(expr):
        return s.expand(
            sum(
                c * s.factorial2(n - 1) * (h / K) ** (n // 2)
                for (n,), c in s.Poly(s.expand(expr), x).terms()
                if n % 2 == 0
            )
        )

    checks = {}
    moment_rows = {}
    for n in range(0, 17, 2):
        expanded = s.expand(jac * weight * F**n)
        values = {j: average(expanded.coeff(r, j)) for j in range(4)}
        moment_rows[n] = values
        for j in range(4):
            checks[f"Gaussian_physical_source_n{n}_r{j}"] = s.expand(
                values[j] - (average(x**n) if j == 0 else 0)
            )
    # General even-monomial identities after factoring its nonzero Gaussian
    # moment and h^j/K^j. Odd monomials vanish by parity.
    n = s.Symbol("n", integer=True, nonnegative=True)
    ratio = lambda k: s.prod(n + 2 * j + 1 for j in range(k))
    general = {
        1: (n + 3) * ratio(1) - ratio(2),
        2: n * (n + 5) * ratio(2) / 2
        - (n + s.Rational(7, 2)) * ratio(3)
        + ratio(4) / 2,
        3: n * (n - 1) * (n + 7) * ratio(3) / 6
        - (n * n + 6 * n + 3) * ratio(4) / 2
        + (n + 4) * ratio(5) / 2
        - ratio(6) / 6,
    }
    for j, v in general.items():
        checks[f"arbitrary_even_monomial_map_order_{j}"] = s.expand(v)
    no_jac = average(weight.coeff(r, 1))
    no_source = average(s.expand(jac * weight * x * x).coeff(r, 1))
    no_sextic = average(
        s.expand(jac * (weight + r * r * K * x**6 / (2 * h))).coeff(r, 2)
    )
    q, J, Y, m = s.symbols("quartic J Y mass", positive=True)
    shifted_eighth = sum(
        s.binomial(8, k) * (h * J / K) ** (8 - k) * average(x**k) for k in range(9)
    )
    no_octic = q * s.diff(shifted_eighth, J, 4).subs(J, 0) / (4 * h)
    moments = {
        n: average(x**n)
        - Y * average(x ** (n + 2)) / m**2
        + 2 * Y * r * average(x ** (n + 4)) / m**2
        for n in (0, 2, 4)
    }
    wrong_connected = moments[4] / moments[0] - 3 * (moments[2] / moments[0]) ** 2
    no_Yukawa = s.simplify(s.diff(wrong_connected, Y, r).subs({Y: 0, r: 0}))
    checks.update(
        {
            "Jacobian_omission_is_nonzero": no_jac + 3 * h / K,
            "physical_source_omission_is_nonzero": no_source + 6 * h * h / K**2,
            "generated_sextic_omission_is_nonzero": no_sextic
            - 15 * h * h / (2 * K * K),
            "generated_octic_omission_changes_two_loop_four_point": no_octic
            - 1260 * q * h**5 / K**6,
            "generated_Yukawa_omission_changes_physical_four_point": no_Yukawa
            - 48 * h**4 / (m**2 * K**4),
            "third_map_order_weight_from_literal_Gaussian_action": s.expand(
                s.diff(s.exp(-K * (F**2 - x**2) / (2 * h)), r, 3).subs(r, 0) / 6
                - weight.coeff(r, 3)
            ),
        }
    )
    return {
        "finite_dimensional_map": F,
        "Gaussian_action_weight_through_map_order_three": weight,
        "physical_source_Gaussian_moments": moment_rows,
        "arbitrary_even_monomial_reduced_identities": general,
        "omission_controls": {
            "Jacobian_partition": no_jac,
            "physical_source_two_point": no_source,
            "sextic_partition": no_sextic,
            "octic_four_point_coefficient_per_r_squared": no_octic,
            "Yukawa_four_point_coefficient_per_r_Y": no_Yukawa,
        },
        "checks": checks,
        "scope": "Finite-dimensional Gaussian diagnostics retain their nonzero Jacobian. General monomial recurrences cover arbitrary polynomial parent weights, including scalar counterterms and the opposite-Yukawa determinant pair. They support, but do not replace, the regulated continuum functional Ward proof.",
    }
