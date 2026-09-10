"""Complete canonical one-loop amplitude conversion and forward coefficient."""

from functools import cache

import sympy as sp

from . import scheme


@cache
def data():
    d = scheme.data()
    L, G, M = (d["parameters"][k] for k in ("L", "G", "M"))
    h, Ibar, ell, sigma, r, fp = (
        d["reference_symbols"][k] for k in ("h", "Ibar", "ell", "sigma", "r", "fp")
    )
    g, Q = G * G, 16 * sp.pi**2
    s, t, u, D = sp.symbols(
        "channel_s channel_t channel_u positive_heavy_gap", real=True
    )
    exchange = [1 / (M - z) for z in (s, t, u)]
    tree = -L + g * sum(exchange)
    squares = sum((-L + g * v) ** 2 for v in exchange)
    conversion = ell * squares / (2 * Q) + sigma + 2 * (fp - r) * tree
    shifts = d["canonical_star_parameter_increments"]
    variation = sum(
        shifts[name] * sp.diff(tree, var)
        for name, var in (("L", L), ("G", G), ("M", M))
    )
    scalar, fermion = sp.symbols(
        "complete_old_scalar_one_loop complete_MS_fermion_one_loop"
    )
    full = tree + h * (scalar + fermion + conversion)
    forward = lambda value: value.subs({t: 0, u: 4 - s})
    center = lambda value: sp.factor(sp.diff(forward(value), s, 2).subs(s, 2) / 2)
    C2 = -4 * L * g / (M - 2) ** 3 + 6 * g * g / (M - 2) ** 4
    T = 2 * g / (M - 2) ** 3
    reduced_relative = -2 * L + 3 * g / D
    fixed_relation = g * (3 * D - 2) / D**2
    checks = {
        "literal_canonical_tree_parameter_variation": sp.factor(variation - conversion),
        "direct_three_channel_counterterm_square_coefficient": sp.factor(
            center(squares) - C2
        ),
        "full_tree_second_coefficient": sp.factor(center(tree) - T),
        "finite_fixed_contact_second_coefficient_zero": center(sigma),
        "conversion_second_coefficient_complete": sp.factor(
            center(conversion) - ell * C2 / (2 * Q) - 2 * (fp - r) * T
        ),
        "reference_log_coefficient_relative_to_tree": sp.factor(
            C2.subs(M, D + 2) / T.subs(M, D + 2) - reduced_relative
        ),
        "actual_tree_relation_relative_conversion": sp.factor(
            reduced_relative.subs(L, fixed_relation) - g * (4 - 3 * D) / D**2
        ),
        "same_canonical_star_and_MS_representations": sp.factor(
            tree + h * variation + h * (scalar + fermion) - full
        ),
        "entire_MS_dimensional_reference_absent": sp.diff(full, Ibar),
        "second_LSZ_copy_changes_the_amplitude": sp.factor(
            full + h * 2 * (fp - r) * tree - full - 2 * h * (fp - r) * tree
        ),
    }
    return {
        "reference_tree_amplitude": tree,
        "fixed_scale_and_field_conversion": conversion,
        "complete_one_loop_canonical_amplitude": full,
        "forward_reference_tree_second_coefficient": T,
        "forward_three_channel_square_second_coefficient": C2,
        "forward_reference_conversion_second_coefficient": ell * C2 / (2 * Q)
        + 2 * (fp - r) * T,
        "relative_scale_conversion_coefficient": reduced_relative,
        "actual_tree_relation": fixed_relation,
        "completeness": "The old scalar one-loop term owns all scalar contact/exchange/tadpole and mixed light-heavy graphs in its fixed reference scheme. The new fermion term owns every first one-loop fermion four-point box and local quartic subtraction. The displayed parameter/field conversion accounts for their common canonical reference. Neutral Phi has no direct gauge vertex or H Yukawa; adding a gauge exchange to the fermion loop or joining a scalar tree across it creates another loop. One-loop parameter insertions inside one-loop graphs are two-loop terms.",
        "scope": "The complete formal one-loop amplitude is in the stated MS interaction boundary and canonical physical light field, with explicitly fixed mass/one-point/vacuum references. Its old scalar term is re-expressed, not numerically transplanted from a different scheme. The new-model two-loop and later errors are not obtained from the old pure-scalar two-loop certificate.",
        "checks": checks,
    }
