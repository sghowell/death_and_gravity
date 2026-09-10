"""Explicit Cauchy subtraction, physical boundary germs and center coefficient."""

from functools import cache

import sympy as sp


@cache
def data():
    s, x, T = sp.symbols("channel_s cut_variable cutoff_T", positive=True)
    u = 4 - s
    primitive = x * x / 2 + s * x + s * s * sp.log(x - s)
    regular = lambda w: T * T / 2 + T * w + w * w * (sp.log(T - w) - sp.log(w))
    # These are local analytic continuations of the two physical boundaries.
    # Their imaginary pieces are opposite, including in the crossed channel.
    plus = regular(s) + regular(u) + sp.I * sp.pi * (s * s - u * u)
    minus = regular(s) + regular(u) - sp.I * sp.pi * (s * s - u * u)
    b2 = sp.diff(plus, s, 2).subs(s, 2) / 2
    formula = 2 * sp.log(T - 2) - 2 * sp.log(2) - 8 / (T - 2) - 4 / (T - 2) ** 2 - 3
    # T>4 ensures positive log arguments at the center.
    b2six = 2 * sp.log(2) - sp.Rational(21, 4)
    rho = sp.Function("rho")
    jump = 2 * sp.I * (rho(s) - rho(4 - s))
    c = sp.symbols("regular_matching_coefficient", real=True)
    checks = {
        "Cauchy_primitive_derivative": sp.factor(
            sp.diff(primitive, x) - x * x / (x - s)
        ),
        "leading_Cauchy_quotient": sp.factor((x * x - s * s) / (x - s) - (x + s)),
        "integrated_leading_quotient": sp.integrate(x + s, (x, 0, T))
        - T * T / 2
        - T * s,
        "upper_minus_lower_leading_germs": sp.expand(
            plus - minus - 2 * sp.I * sp.pi * (s * s - u * u)
        ),
        "crossing_exchanges_boundary_germs": sp.simplify(plus.subs(s, u) - minus),
        "general_center_second_coefficient": sp.simplify(b2 - formula),
        "fixed_cutoff_six_coefficient": sp.simplify(formula.subs(T, 6) - b2six),
        "cutoff_derivative_from_endpoint_kernel": sp.factor(
            sp.diff(formula, T) - 2 * T * T / (T - 2) ** 3
        ),
        "kernel_second_coefficient": sp.factor(
            (sp.diff(1 / (x - s) + 1 / (x - u), s, 2) / 2).subs(s, 2) - 2 / (x - 2) ** 3
        ),
        "general_cut_jump_crossing_odd": sp.simplify(jump.subs(s, u) + jump),
        "regular_polynomial_preserves_crossing": sp.expand(
            c * (u - 2) ** 2 - c * (s - 2) ** 2
        ),
        "same_cut_polynomial_changes_b2": sp.diff(c * (s - 2) ** 2, s, 2) / 2 - c,
        "same_cut_polynomial_has_no_boundary_jump": c * (s - 2) ** 2 - c * (s - 2) ** 2,
    }
    for degree in (0, 2, 4, 6):
        checks[f"general_jump_even_center_derivative_{degree}"] = sp.simplify(
            sp.diff(jump, s, degree).subs(s, 2)
        )
    return {
        "s": s,
        "u": u,
        "fixed_upper_cut_invariant": 6,
        "symbolic_cutoff_domain": "T>4; the published subtraction uses T=6",
        "Cauchy_kernel": 1 / (x - s) + 1 / (x - u),
        "finite_cut_subtraction": "F_low(s)=pi^-1 integral_0^6 rho(x)[(x-s)^-1+(x-(4-s))^-1] dx, initially Im(s)!=0",
        "Plemelj_divided_difference": "C_rho(s)=pi^-1[integral_0^6 (rho(x)-rho(s))/(x-s) dx +rho(s)(Log(6-s)-Log(-s))]",
        "leading_primitive_over_K": primitive,
        "leading_upper_boundary_germ_over_K": plus.subs(T, 6),
        "leading_lower_boundary_germ_over_K": minus.subs(T, 6),
        "leading_general_cutoff_b2_over_K": formula,
        "leading_fixed_cutoff_b2_over_K": b2six,
        "general_first_cut_boundary_jump": jump,
        "regular_part_ambiguity_diagnostic": c * (s - 2) ** 2,
        "scope": "The local coefficient belongs to the specified finite-cut subtraction, not the full amplitude. Its two physical boundary germs have the same even derivatives at the crossing center. Adding an analytic crossing-even term preserves the cut but changes b2: this is an insufficiency diagnostic, not permission to retune the fixed renormalizable model.",
        "checks": checks,
    }
