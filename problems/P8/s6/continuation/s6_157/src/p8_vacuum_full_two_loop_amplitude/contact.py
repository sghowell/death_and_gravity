"""Pointwise isolated-contact cancellation survives subsequent field matching."""

from functools import cache

import sympy as s


@cache
def data():
    h, L, g, k, c, u, v, w = s.symbols("h L g k c u v w")
    sigma = -k * L + c
    star = (L - h * c) / (1 - h * k)
    f = u * L**2 + v * L * g + w * g**2
    Fstar = f.subs(L, star)
    sigstar = sigma.subs(L, star)
    derivative = s.diff(f, L).subs(L, star)
    old = -star + h * (Fstar - sigstar) + h**2 * sigstar * derivative
    direct = -L + h * f
    coeffs = [s.factor(s.diff(old - direct, h, n).subs(h, 0)) for n in range(3)]
    # No commutation assumption is needed. Even noncommuting directions
    # act on a pointwise zero remainder through the stated order.
    field = lambda value: 2 * L * s.diff(value, L) + 2 * g * s.diff(value, g)
    sigma_g = -g * L + g**2
    comm = s.expand(field(sigma_g) - 2 * sigma_g)
    return {
        "regulated_affine_inverse": star,
        "isolated_finite_contact": sigma,
        "nonlocal_first_loop_representative": f,
        "old_and_direct_functionals_through_order_two": {"old": old, "direct": direct},
        "composition_rule": "An analytic pointwise R(h,q)=O(h^3) remains O(h^3) under q=q0+O(h), before removing the common regulator. No commutation of the isolated contact and field directions is presumed.",
        "nonzero_direction_commutator_example": comm,
        "checks": {
            "exact_contact_inverse": s.factor(star + h * sigstar - L),
            "zero_tree_difference": coeffs[0],
            "zero_one_loop_difference": coeffs[1],
            "zero_two_loop_nonlocal_difference": coeffs[2],
            "noncommuting_example_evaluated": comm - (-2 * g * L + 2 * g**2),
            "no_new_sigma_insertion_after_direct_MS_conversion": s.diff(coeffs[2], g),
        },
    }
