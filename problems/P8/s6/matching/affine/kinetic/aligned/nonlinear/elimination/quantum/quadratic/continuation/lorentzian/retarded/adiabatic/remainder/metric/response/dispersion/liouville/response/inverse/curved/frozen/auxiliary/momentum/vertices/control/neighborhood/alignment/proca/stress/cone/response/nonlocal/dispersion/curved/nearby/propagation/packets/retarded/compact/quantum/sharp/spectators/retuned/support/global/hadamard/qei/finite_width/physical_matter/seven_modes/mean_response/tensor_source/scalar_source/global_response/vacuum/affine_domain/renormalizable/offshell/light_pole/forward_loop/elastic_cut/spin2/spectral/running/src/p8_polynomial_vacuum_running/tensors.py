"""Literal two-field tensor beta coefficients and existing-counterterm match."""

from functools import cache
from itertools import product

import sympy as sp
from p8_vacuum_forward_loop import subtraction


@cache
def data():
    L, G, M, m = sp.symbols(
        "polynomial_quartic cubic_coupling heavy_mass_squared light_mass_squared",
        real=True,
    )
    g = G * G

    def lam(a, b, c, d):
        return L if (a, b, c, d) == (0, 0, 0, 0) else sp.S.Zero

    def cubic(a, b, c):
        return G if sorted((a, b, c)) == [0, 0, 1] else sp.S.Zero

    def mass(a, b):
        return (m if a == 0 else M) if a == b else sp.S.Zero

    B4 = {
        indices: sp.expand(
            sum(
                lam(a, b, e, f) * lam(e, f, c, d)
                + lam(a, c, e, f) * lam(e, f, b, d)
                + lam(a, d, e, f) * lam(e, f, b, c)
                for e, f in product(range(2), repeat=2)
            )
        )
        for indices in product(range(2), repeat=4)
        for a, b, c, d in [indices]
    }
    B3 = {
        indices: sp.expand(
            sum(
                lam(a, b, e, f) * cubic(e, f, c)
                + lam(a, c, e, f) * cubic(e, f, b)
                + lam(b, c, e, f) * cubic(e, f, a)
                for e, f in product(range(2), repeat=2)
            )
        )
        for indices in product(range(2), repeat=3)
        for a, b, c in [indices]
    }
    B2 = {
        indices: sp.expand(
            sum(
                mass(e, f) * lam(a, b, e, f) + cubic(a, e, f) * cubic(b, e, f)
                for e, f in product(range(2), repeat=2)
            )
        )
        for indices in product(range(2), repeat=2)
        for a, b in [indices]
    }
    checks = {}
    checks.update(
        {
            "quartic_tensor_" + "".join(map(str, k)): v
            - (3 * L * L if k == (0, 0, 0, 0) else 0)
            for k, v in B4.items()
        }
    )
    checks.update(
        {
            "cubic_tensor_" + "".join(map(str, k)): v
            - (L * G if sorted(k) == [0, 0, 1] else 0)
            for k, v in B3.items()
        }
    )
    checks.update(
        {
            "mass_tensor_" + "".join(map(str, k)): v
            - ({(0, 0): L * m + 2 * g, (1, 1): g}.get(k, 0))
            for k, v in B2.items()
        }
    )
    old = subtraction.data()
    expressions = [
        old[k]
        for k in (
            "delta_polynomial_quartic_UV",
            "delta_cubic_squared_UV",
            "delta_heavy_mass_squared_UV",
        )
    ]
    symbols = {str(s): s for expression in expressions for s in expression.free_symbols}
    mapping = {symbols["quartic_L"]: L, symbols["cubic_squared"]: g}
    I0 = symbols["regulated_zero_bubble"]
    beta = (3 * L * L, 2 * L * g, g)
    for name, expression, target in zip(
        ("quartic", "cubic_squared", "heavy_mass_squared"), expressions, beta
    ):
        checks["actual_parent_" + name + "_counterterm_factor"] = sp.factor(
            32 * sp.pi**2 * expression.subs(mapping) / I0 - target
        )
    return {
        "L": L,
        "G": G,
        "M": M,
        "light_mass_squared": m,
        "literal_quartic_beta_tensor": B4,
        "literal_cubic_beta_tensor": B3,
        "literal_mass_beta_tensor": B2,
        "one_loop_beta_in_log_reference_over_16pi2": {
            "quartic": beta[0],
            "cubic_squared": beta[1],
            "heavy_mass_squared": beta[2],
        },
        "light_mass_convention": "The tensor UV coefficient of the light quadratic term is retained as a diagnostic. The actual light pole mass is kept one by its existing on-shell counterterm, not evolved with a new MS mass prescription.",
        "scope": "Pure scalar one-loop logarithmic coefficients. No higher-loop, gravitational or independent finite-scheme conversion coefficient is imported.",
        "checks": checks,
    }
