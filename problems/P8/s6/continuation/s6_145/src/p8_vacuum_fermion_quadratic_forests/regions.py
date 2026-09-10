"""Symmetric middle/high regions with both overlapping vertex subtractions."""

from functools import cache

import sympy as sp


@cache
def data():
    u, x, y, m, R = sp.symbols("u x y m R", positive=True)
    I0 = sp.integrate(1 - u, (u, 0, 1))
    Ilog = sp.integrate(-(1 - u) * sp.log(u), (u, 0, 1))
    selfradial = Ilog + 6 * I0
    middle = 2 * (4 * 2 * I0 + 2 * I0 + 2 * I0)
    high = 2 * (32 * I0 + Ilog)
    anchors = 2 * 2 * I0
    vertex = middle + high + anchors
    base = 32 * 360**4
    return {
        "regions": "Use symmetry q<->l. Half-domain x=q^2<=y=l^2. MIDDLE x<=y<=4(m^2+x); HIGH y>=4(m^2+x). Then multiply by two.",
        "raw_MIDDLE_radius": "sqrt(m^2+x)/360; x<=y",
        "first_subtraction_MIDDLE_radius": "sqrt(m^2+x)/360; K0(l) has no soft momentum",
        "second_subtraction_radius": "sqrt(m^2+y)/360; K0(q) has no soft momentum, in BOTH middle and high regions",
        "paired_HIGH_radius": "sqrt(m^2+x)/360 for raw minus the first subtraction; inherited paired kernel norm <=32 sqrt(m^2+x)/y^(5/2)",
        "self_energy_radial_moment": selfradial,
        "symmetric_middle_constant": middle,
        "symmetric_high_constant": high,
        "two_local_anchor_constant": anchors,
        "complete_vertex_constant": vertex,
        "vertex_base": base,
        "exact_vertex_prefactor": base * vertex,
        "exact_self_prefactor": 2 * 64 * 360**4 * selfradial,
        "checks": {
            "compact_radial_density": sp.factor(
                ((1 - u) / u) / (1 + (1 - u) / u) ** 3 / u**2 - (1 - u)
            ),
            "basic_radial_moment": I0 - sp.Rational(1, 2),
            "log_radial_moment": Ilog - sp.Rational(3, 4),
            "self_radial_moment": selfradial - sp.Rational(15, 4),
            "three_self_propagators_trace": 4 * 2**3 - 32,
            "two_self_marks": 2 * 64 * 360**4 * selfradial - sp.Integer(8062156800000),
            "raw_vertex_trace_to_base_ratio": 4 * 2**4 / (4 * sp.Integer(2) ** 2) - 4,
            "minimum_radius_on_half_domain": (m * m + y) - (m * m + x) - (y - x),
            "middle_log_ratio": sp.factor(
                (m * m + 4 * (m * m + x)) / (m * m + x) - (4 + m * m / (m * m + x))
            ),
            "middle_opposite_log_small_y": (m * m + 4 * m * m) / m**2 - 5,
            "middle_opposite_log_large_y": 4 * (1 + sp.Rational(1, 4)) - 5,
            "symmetric_middle_three_terms": middle - 12,
            "high_paired_inner_integral": sp.integrate(
                y ** (-sp.Rational(3, 2)), (y, 4 * (m * m + x), sp.oo)
            )
            - 1 / sp.sqrt(m * m + x),
            "high_complement_outer_log": sp.integrate(
                1 / (m * m + x), (x, 0, y / 4 - m * m)
            )
            - sp.log(y / (4 * m * m)),
            "symmetric_high_includes_complement": high - sp.Rational(67, 2),
            "both_MS_anchor_insertions": anchors - 2,
            "complete_vertex_radial_factor": vertex - sp.Rational(95, 2),
            "exact_vertex_base": base - 537477120000,
            "exact_vertex_prefactor": base * vertex - 25530163200000,
            "vertex_rounding_slack": 3 * 10**13 - base * vertex - 4469836800000,
            "self_rounding_slack": 10**13
            - 2 * 64 * 360**4 * selfradial
            - 1937843200000,
            "soft_geometric_tail_gap": sp.factor(
                2 / R**4 - 1 / (R**4 * (1 - 1 / R)) - (R - 2) / (R**4 * (R - 1))
            ),
        },
    }
