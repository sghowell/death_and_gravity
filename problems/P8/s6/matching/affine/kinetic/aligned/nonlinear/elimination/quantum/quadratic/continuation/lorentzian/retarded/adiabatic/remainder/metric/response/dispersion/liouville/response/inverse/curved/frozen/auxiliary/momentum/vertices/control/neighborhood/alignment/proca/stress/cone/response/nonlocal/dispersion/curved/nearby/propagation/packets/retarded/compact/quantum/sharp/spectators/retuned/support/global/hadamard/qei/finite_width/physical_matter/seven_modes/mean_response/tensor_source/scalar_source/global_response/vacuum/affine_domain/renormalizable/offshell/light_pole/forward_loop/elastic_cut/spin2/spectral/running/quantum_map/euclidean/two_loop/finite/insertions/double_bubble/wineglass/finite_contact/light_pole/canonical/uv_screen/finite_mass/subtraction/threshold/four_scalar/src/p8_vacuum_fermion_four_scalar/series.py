"""All higher momentum degrees of the finite four-scalar box."""

from functools import cache
from itertools import permutations

import sympy as sp
from p8_vacuum_finite_mass_gauge_cut import series as parent_series


@cache
def data():
    n = sp.symbols("integer_total_insertions", integer=True, positive=True)
    x = sp.symbols("nonnegative_Neumann_ratio", nonnegative=True)
    m, Y = sp.symbols("positive_fermion_mass positive_Yukawa_squared", positive=True)
    Q = 16 * sp.pi**2
    orders = [("Phi1",) + p for p in permutations(("Phi2", "Phi3", "Phi4"))]
    weight = (n + 3) * (n + 1) / (6 * n)
    geometric = x**4 / (1 - x)
    majorant = sp.factor((x * sp.diff(geometric, x) + 5 * geometric) / 6)
    ratio = sp.factor(majorant / x**4)
    checks = {
        "six_cyclic_four_scalar_orders": len(orders) - 6,
        "same_propagator_radial_integral": parent_series.data()[
            "all_degree_radial_integral"
        ]
        - 4 / (Q * n * (n + 2) * m**n),
        "weak_composition_weight_reduction": sp.factor(
            sp.binomial(n + 3, 3).expand(func=True) / (n * (n + 2)) - weight
        ),
        "weight_majorant_positive_residual": sp.factor(
            (n + 5) / 6 - weight - (n - 3) / (6 * n)
        ),
        "all_degree_four_and_higher_majorant": sp.factor(
            majorant - x**4 * (9 - 8 * x) / (6 * (1 - x) ** 2)
        ),
        "majorant_ratio_monotonic_derivative": sp.factor(
            sp.diff(ratio, x) - (10 - 8 * x) / (6 * (1 - x) ** 3)
        ),
        "half_ratio_exact_upper": ratio.subs(x, sp.Rational(1, 2)) - sp.Rational(10, 3),
        "strict_four_majorant_margin": 4 - sp.Rational(10, 3) - sp.Rational(2, 3),
        "active_color_spin_cycle_radial_prefactor": 6 * 4 * 6 * 4 - 576,
        "unrounded_full_momentum_tail_constant": 576 * 4 * 18**4 - 241864704,
        "strict_rounded_tail_margin": 300000000 - 241864704 - 58135296,
        "unit_disc_Cauchy_radius_factor": sp.Integer(1) ** 2 - 1,
    }
    rows = []
    for degree in range(11):
        values = parent_series.compositions(degree, 4)
        checks[f"independent_insertions_degree_{degree}"] = len(values) - sp.binomial(
            degree + 3, 3
        )
        checks[f"distinct_insertions_degree_{degree}"] = len(values) - len(set(values))
        rows.append(
            {
                "degree": degree,
                "enumerated_count": len(values),
                "formula_count": sp.binomial(degree + 3, 3),
            }
        )
    return {
        "cyclic_labelled_four_scalar_orders": orders,
        "ordered_insertion_count": sp.binomial(n + 3, 3),
        "independent_small_degree_counts": rows,
        "all_degree_radial_integral": 4 / (Q * n * (n + 2) * m**n),
        "complete_tail_majorant": majorant,
        "unrounded_four_scalar_tail_upper": 241864704 * Y**2 / (Q * m**4),
        "four_scalar_tail_upper": 300000000 * Y**2 / (Q * m**4),
        "complex_tail_domain": "|s-2|<=5 at forward transfer, all external masses one; cumulative component norm <18 and mF>=36.",
        "coefficient_bound": "The degree-zero and degree-two terms are constant on shell. The remaining holomorphic function obeys the displayed sup bound, hence its center second coefficient obeys the same bound by Cauchy on |s-2|=1.",
        "scope": "All momentum degrees n>=4 are bounded, including an overcount of the identically vanishing odd degrees. This is the complete tail of the one-loop fermion box, not a bound on later primitive loop orders.",
        "checks": checks,
    }
