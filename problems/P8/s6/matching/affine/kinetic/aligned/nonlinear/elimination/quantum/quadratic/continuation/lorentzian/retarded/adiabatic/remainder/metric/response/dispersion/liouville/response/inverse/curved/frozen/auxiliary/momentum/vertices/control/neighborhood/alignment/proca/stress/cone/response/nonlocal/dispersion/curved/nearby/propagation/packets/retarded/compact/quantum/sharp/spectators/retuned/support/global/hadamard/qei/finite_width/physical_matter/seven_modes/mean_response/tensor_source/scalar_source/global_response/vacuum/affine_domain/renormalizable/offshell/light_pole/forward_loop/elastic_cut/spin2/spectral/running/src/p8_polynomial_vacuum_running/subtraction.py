"""Exact reference-bubble change and full-channel one-loop compensation."""

from functools import cache

import sympy as sp


@cache
def data():
    y, r, nu = sp.symbols(
        "positive_radial_variable positive_reference_mass_squared positive_reference_mass",
        positive=True,
    )
    difference = y / (y + r) ** 2 - y / (y + 1) ** 2
    primitive = sp.log(y + r) + r / (y + r) - sp.log(y + 1) - 1 / (y + 1)
    derivative_integral = sp.integrate(-2 * y / (y + r) ** 3, (y, 0, sp.oo))
    L, g, M = sp.symbols(
        "polynomial_quartic cubic_squared heavy_mass_squared", real=True
    )
    s, t, u = sp.symbols("s t u", real=True)
    h = [1 / (M - z) for z in (s, t, u)]
    A = -L + g * sum(h)
    C = [-L + g * v for v in h]
    beta = (3 * L * L, 2 * L * g, g)
    tree_flow = sum(v * sp.diff(A, p) for v, p in zip(beta, (L, g, M)))
    d = L - 3 * g / M
    dflow = sum(v * sp.diff(d, p) for v, p in zip(beta, (L, g, M)))
    return {
        "reference_bubble_integrand_difference": difference,
        "anchored_reference_bubble_primitive": primitive,
        "reference_bubble_finite_change": -sp.log(r),
        "full_tree_amplitude": A,
        "explicit_one_loop_reference_derivative": sum(v * v for v in C),
        "tree_parameter_flow_derivative": tree_flow,
        "effective_quartic_margin_derivative": dflow,
        "scheme": "Change only the common zero-bubble reference from mass one to nu. The finite potential contact is fixed once at nu=1 and transported, not retuned to impose a fresh zero-potential-loop condition at every nu. Light pole mass and LSZ residue remain on shell.",
        "order": "The cancellation is through one matter loop: parameter derivatives of a one-loop term are two-loop contributions and are not set to zero or bounded here.",
        "checks": {
            "literal_reference_bubble_primitive": sp.factor(
                sp.diff(primitive, y) - difference
            ),
            "reference_bubble_derivative_under_integral": derivative_integral + 1 / r,
            "reference_bubble_zero_anchor": difference.subs(r, 1),
            "reference_bubble_upper_endpoint": sp.limit(primitive, y, sp.oo),
            "reference_bubble_lower_endpoint": sp.simplify(
                primitive.subs(y, 0) - sp.log(r)
            ),
            "reference_mass_logarithm_factor_two": sp.diff(
                -sp.log(r).subs(r, nu * nu), nu
            )
            * nu
            + 2,
            "all_channel_one_loop_reference_compensation": sp.factor(
                tree_flow + sum(v * v for v in C)
            ),
            "effective_quartic_running_is_positive_square": sp.factor(
                dflow - 3 * (L - g / M) ** 2
            ),
            "heavy_mass_flow_cannot_be_omitted": sp.factor(
                (tree_flow - g * sp.diff(A, M))
                + sum(v * v for v in C)
                - g * g * sum(v * v for v in h)
            ),
        },
    }
