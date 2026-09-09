"""Exact boost family, normalized compact shapes and two-null-width scaling."""

from functools import cache

import sympy as sp

from . import flat


@cache
def data():
    b, Hp, Hpp, Jp, Jpp = sp.symbols("boost H1 H2 J1 J2", positive=True)
    xi = sp.Symbol("xi", real=True)
    # Product profiles have vanishing first/second mixed cross integrals,
    # while <f_++,f_-->=H1*J1. These are exact compact IBP identities.
    gram = sp.Matrix([[Hpp, 0, Hp * Jp], [0, Hp * Jp, 0], [Hp * Jp, 0, Jpp]])
    U2 = sp.Matrix([b * b, 2, b**-2])
    UE = sp.Matrix([b * b, 0, -(b**-2)])
    UU = (U2.T * gram * U2)[0]
    EE = (UE.T * gram * UE)[0]
    cross = (U2.T * gram * UE)[0]
    form = flat.data()["proper_derivative_quadratic_form"].subs(flat.data()["xi"], xi)
    cost = sp.expand(
        b**-2 * (form[0, 0] * UU + 2 * form[0, 1] * cross + form[1, 1] * EE)
    )
    dp, dm, r = sp.symbols("delta_plus delta_minus ratio", positive=True)
    scaled = sp.factor(
        cost.subs(
            {Hp: Hp / dp**2, Hpp: Hpp / dp**4, Jp: Jp / dm**2, Jpp: Jpp / dm**4},
            simultaneous=True,
        ).subs(b, sp.sqrt(r * dp / dm))
    )
    s = sp.Symbol("shape_coordinate", real=True)
    shape = sp.sqrt(sp.Rational(315, 256)) * (1 - s * s) ** 2
    h1 = sp.integrate(sp.diff(shape, s) ** 2, (s, -1, 1))
    h2 = sp.integrate(sp.diff(shape, s, 2) ** 2, (s, -1, 1))
    selected = sp.factor(scaled * dp**3 * dm).subs(
        {xi: sp.Rational(1, 6), Hp: h1, Hpp: h2, Jp: h1, Jpp: h2}
    )
    selected = sp.simplify(selected)
    r2 = (12 + 9 * sp.sqrt(29)) / 35
    # Evaluate the actual polynomial Gram independently of the formal IBP rows.
    pairings = sp.Matrix(
        3,
        3,
        lambda i, j: sp.integrate(
            sp.diff(shape, s, i) * sp.diff(shape, s, j), (s, -1, 1)
        ),
    )
    orders = [(2, 0), (1, 1), (0, 2)]
    actual_gram = sp.Matrix(
        3,
        3,
        lambda i, j: (
            pairings[orders[i][0], orders[j][0]] * pairings[orders[i][1], orders[j][1]]
        ),
    )
    return {
        "product_derivative_Gram": gram,
        "boosted_UU_norm": UU,
        "boosted_UE_norm": EE,
        "boosted_cross": cross,
        "full_boost_cost_without_hbar_over_8pi_squared": cost,
        "normalized_compact_shape": shape,
        "shape_first_norm_squared": h1,
        "shape_second_norm_squared": h2,
        "conformal_product_cost_polynomial": selected,
        "exact_boost_minimizer_squared": r2,
        "rational_boost_cost": selected.subs(r, sp.Rational(4, 3)),
        "actual_shape_derivative_pairings": pairings,
        "actual_product_derivative_Gram": actual_gram,
        "checks": {
            "complete_product_UU_norm": sp.expand(
                UU - b**4 * Hpp - 6 * Hp * Jp - b**-4 * Jpp
            ),
            "complete_product_UE_norm": sp.expand(
                EE - b**4 * Hpp + 2 * Hp * Jp - b**-4 * Jpp
            ),
            "complete_product_cross": sp.expand(cross - b**4 * Hpp + b**-4 * Jpp),
            "complete_nonminimal_boost_family": sp.expand(
                cost
                - sp.Rational(2, 3) * (1 + 4 * xi) * b * b * Hpp
                - 4 * (1 - 2 * xi) * b**-2 * Hp * Jp
                - sp.Rational(2, 3) * b**-6 * Jpp
            ),
            "both_null_widths_retained": sp.factor(
                scaled * dp**3 * dm
                - (
                    sp.Rational(2, 3) * (1 + 4 * xi) * r * Hpp
                    + 4 * (1 - 2 * xi) * Hp * Jp / r
                    + sp.Rational(2, 3) * Jpp / r**3
                )
            ),
            "actual_compact_shape_normalization": sp.integrate(
                shape * shape, (s, -1, 1)
            )
            - 1,
            "actual_compact_shape_first_norm": h1 - 3,
            "actual_compact_shape_second_norm": h2 - sp.Rational(63, 2),
            "actual_polynomial_product_Gram_matches_IBP": actual_gram
            - gram.subs({Hp: h1, Hpp: h2, Jp: h1, Jpp: h2}),
            "compact_shape_outer_value_traces": shape.subs(s, 1) ** 2
            + shape.subs(s, -1) ** 2,
            "compact_shape_outer_first_derivative_traces": sp.diff(shape, s).subs(s, 1)
            ** 2
            + sp.diff(shape, s).subs(s, -1) ** 2,
            "conformal_product_exact_cost": sp.factor(
                selected - 35 * r - 24 / r - 21 / r**3
            ),
            "conformal_rational_boost_cost": sp.factor(
                selected.subs(r, sp.Rational(4, 3)) - sp.Rational(14117, 192)
            ),
            "conformal_boost_optimizer_equation": sp.simplify(
                35 * r2 * r2 - 24 * r2 - 63
            ),
            "strict_boost_convexity_formula": sp.factor(
                sp.diff(selected, r, 2) - 48 / r**3 - 252 / r**5
            ),
            "state_cost_is_boost_independent": sp.factor(
                2 * xi * b**-2 * (4 * b * b * Hp) - 8 * xi * Hp
            ),
            "actual_conformal_shape_state_cost": 8 * sp.Rational(1, 6) * h1 - 4,
        },
    }
