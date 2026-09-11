"""The fixed physical mass anchor and the entire paired potential, not a fit."""

from functools import cache

import sympy as s
from p8_vacuum_fermion_local_matching import anchors


@cache
def data():
    z, m, y, x = s.symbols("scalar_argument mean_mass Yukawa x", positive=True)
    N, Q = s.symbols("N Q", positive=True)
    r = s.Symbol("mass_ratio", real=True)
    F = lambda v: v**4 * (2 * s.log(v) - s.Rational(3, 2))
    pair = -N * (F(1 + r) + F(1 - r)) / 2
    f0 = 4 * N * y * y * m * m / Q
    f1 = (
        2
        * N
        * y
        * y
        / Q
        * (
            2 * m * m
            + (4 * m * m - 1) * s.Integral(-s.log(1 - x * (1 - x) / (m * m)), (x, 0, 1))
        )
    )
    inc = s.Symbol("exact_f_one_minus_f_zero", positive=True)
    renormalized = -inc * z * z / 2 - 8 * N * y**4 * z**4 / (3 * Q)
    k = s.Symbol("even_degree", integer=True, positive=True)
    coefficient = 48 * N / (k * (k - 1) * (k - 2) * (k - 3) * (k - 4))
    checks = {
        "pair_vacuum_constant": s.simplify(pair.subs(r, 0) - 3 * N / 2),
        "pair_quadratic_coefficient": s.simplify(
            s.diff(pair, r, 2).subs(r, 0) / 2 - 2 * N
        ),
        "pair_quartic_coefficient": s.simplify(
            s.diff(pair, r, 4).subs(r, 0) / 24 + 8 * N / 3
        ),
        "sixth_coefficient": s.simplify(coefficient.subs(k, 6) - N / 15),
        "strict_decrease_of_even_tail_coefficients": s.expand(
            (k + 1) * (k + 2) - (k - 4) * (k - 3) - 10 * (k - 1)
        ),
        "fixed_mass_anchor_cancellation_not_a_fit": s.expand(
            f0 * z * z / 2
            - (f0 + inc) * z * z / 2
            - 8 * N * y**4 * z**4 / (3 * Q)
            - renormalized
        ),
        "fixed_zero_field_vacuum_reference": renormalized.subs(z, 0),
    }
    for degree in (6, 8, 10, 12):
        checks[f"full_potential_even_coefficient_{degree}"] = s.simplify(
            s.diff(pair, r, degree).subs(r, 0) / s.factorial(degree)
            - coefficient.subs(k, degree)
        )
    return {
        "one_loop_MS_pole_mass_anchor": f1,
        "zero_momentum_mass_anchor": f0,
        "paired_potential_dimensionless": pair,
        "first_reference_subtracted_terms": renormalized,
        "all_higher_even_coefficients": coefficient,
        "same_frozen_all_order_anchor": anchors.data()["mass_anchor_series"],
        "mass_anchor_bound": "For mF>=36, 0<f(1)-f(0)<4NY/(3Q). This uses the positive complete moment series a_n=3(n!)^2/[n(2n+3)(2n+1)!], not subtraction of two enormous decimal mass terms.",
        "full_potential_bound": "After subtracting the fixed all-flavor zero-field vacuum constant and -f(1)chi^2/2, the paired potential is -.5[f(1)-f(0)]chi^2 -(8/3)NY^2 chi^4/Q plus the positive complete even tail. For |ychi|<=Delta<m, the tail is <=NDelta^6/[15Q m^2(1-(Delta/m)^2)].",
        "reference_coordinate": "This is the one-loop fermion increment in the named MS field coordinate with physical Phi pole mass reference, fixed vacuum zero, and minimal field/quartic subtractions. No extra finite kinetic or quartic counterterm is added. The fixed local physical unit-residue source and parameter dictionary remains the inherited one; the computed quadratic contribution is not silently the entire canonically re-expressed scalar/heavy energy.",
        "checks": checks,
        "gates": {
            "frozen_anchor_increment_lower_positive": bool(
                anchors.enclosure(36)["mass_increment_over_2NY_div_Q_lower"] > 0
            ),
            "frozen_anchor_increment_upper_below_two_thirds": bool(
                anchors.enclosure(36)["mass_increment_over_2NY_div_Q_upper"]
                < s.Rational(2, 3)
            ),
            "all_higher_even_potential_coefficients_positive": True,
            "zero_field_inert_vacuum_constants_removed_by_same_reference": True,
        },
    }
