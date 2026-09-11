"""Finite-difference Dirac norm bound after the complete local MS reference."""

from functools import cache

import sympy as s


@cache
def data():
    phi, m, y, N, Q = s.symbols("Phi m y N Q", positive=True)
    potential = (
        -N
        * (m + y * phi) ** 4
        * (s.log((m + y * phi) ** 2 / m**2) - s.Rational(3, 2))
        / Q
    )
    constant = s.simplify(-s.diff(potential, phi, 4).subs(phi, 0))
    radial = s.simplify(s.gamma(s.Rational(1, 2)) / s.gamma(s.Rational(5, 2)))
    return {
        "zero_momentum_MS_scattering_vertex": constant,
        "six_cyclic_words": s.factorial(4) / 4,
        "finite_difference_bound": "|A_F(p)-A_F(0)| <=768 N Y^2/[Q m (1-6/m)^4] <=3072 N Y^2/(Q m), m>=24.",
        "physical_route_norm": "For physical s in [4,6], each Wick-continued external Dirac matrix has norm at most |E|+|p|<2. Each unshifted route is a sum of at most three externals, so its perturbation norm is below six.",
        "Dirac_resolvent_bound": "On real Euclidean q, ||S0(q)||=(q^2+m^2)^-1/2. Telescoping the four-propagator product gives 4*6 ||S0||^5/(1-6/m)^4. The six cyclic words and trace dimension four are retained.",
        "checks": {
            "full_fixed_MS_zero_vertex": constant - 64 * N * y**4 / Q,
            "all_six_cyclic_words": s.factorial(4) / 4 - 6,
            "UV_convergent_four_dimensional_radial_integral": radial - s.Rational(4, 3),
            "trace_words_product_and_radial_multiplicity": 6 * 4 * 4 * 6 * radial - 768,
            "route_denominator_bound_at_m24": (1 - s.Rational(6, 24)) ** (-4)
            - s.Rational(256, 81),
            "rational_finite_difference_simplification": 768 * 4 - 3072,
            "local_MS_vertex_not_zero": s.diff(constant, N) - 64 * y**4 / Q,
        },
        "scope": "The full dimensional zero-momentum local term is removed before the convergent four-dimensional difference is bounded. No high-energy momentum expansion, regulator-dependent rational term, second field conversion or inert-flavor Yukawa is discarded.",
    }
