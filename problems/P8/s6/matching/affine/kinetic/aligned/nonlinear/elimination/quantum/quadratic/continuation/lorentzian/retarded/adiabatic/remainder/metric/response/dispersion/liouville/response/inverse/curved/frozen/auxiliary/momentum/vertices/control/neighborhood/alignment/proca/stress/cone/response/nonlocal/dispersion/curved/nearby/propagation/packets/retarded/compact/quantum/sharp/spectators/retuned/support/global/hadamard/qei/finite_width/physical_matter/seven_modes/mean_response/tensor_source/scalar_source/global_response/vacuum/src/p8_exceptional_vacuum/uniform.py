"""Uniform clock-tube control along the interaction-preserving even-order family."""

from functools import cache

import sympy as sp

from . import analytic, family


def localized_rational_norm(value):
    u, X = family.u, family.X
    if isinstance(value, bool) or not isinstance(value, (int, sp.Basic)):
        raise TypeError("The envelope requires an exact rational expression")
    value = sp.sympify(value)
    if value.has(sp.Float, sp.oo, -sp.oo, sp.zoo, sp.nan) or value.free_symbols - {
        u,
        X,
    }:
        raise ValueError(
            "Only finite exact rational functions of u and X are supported"
        )
    if not value.is_rational_function(u, X):
        raise ValueError("The envelope requires a rational function of u and X")
    rows = {}
    for a in range(analytic.JET_ORDER + 1):
        for b in range(analytic.JET_ORDER + 1 - a):
            num, den = sp.fraction(sp.cancel((1 + u * u) ** 6 * sp.diff(value, X, b)))
            D = sp.Poly(den, u, domain=sp.QQ)
            if D.nth(0) <= 0 or any(
                power[0] % 2 or coefficient < 0 for power, coefficient in D.terms()
            ):
                raise ValueError(
                    "A positive even denominator is required for this global envelope"
                )
            P = sp.Poly(num, u, X, domain=sp.QQ)
            rows[(a, b)] = sum(
                abs(c) * sp.Rational(11, 10) ** p[1] * sp.factorial((p[0] + 3) // 4)
                for p, c in P.terms()
            ) / (D.nth(0) * sp.factorial(a) * sp.factorial(b))
        value = sp.cancel(sp.diff(value, u) - 4 * u**3 * value)
    return rows


@cache
def data():
    u, X = family.u, family.X
    Ftree = family.data()["original_retuned_tree_scalar"]
    F0 = Ftree.subs({u: 0, X: 0})
    linear = (analytic.VACUUM_LAMBDA_BAR / analytic.FIXED_GAMMA - F0) * X**2 - u**4 / 3
    constant = X / 2 - u**2 / 2 - F0 * u**4 - Ftree
    n = sp.Symbol("even_switch_order_at_least_1024", positive=True, integer=True)
    Glin = localized_rational_norm(linear)
    Gconst = localized_rational_norm(constant)
    Gu = analytic.ORDER * sum(Glin.values()) + sum(Gconst.values())
    scalar_upper = analytic.coefficient_error_bounds()["S_jet_norm"] * Gu
    # For 0<=j<=6 and every factor of (n)_{j-b}, increasing n by 2
    # multiplies it by at most L=1+2/(1024-5). The extra linear n
    # in the scalar envelope contributes one further L.
    L = 1 + sp.Rational(2, analytic.ORDER - 5)
    ratio = analytic.TUBE_Q_UPPER**2 * L**7
    # Negative-X and two positive-X regions, valid for each n separately.
    floors = {
        "small_negative_X": 1 - sp.Rational(5, 64) / n,
        "small_nonnegative_X": 1 - sp.Rational(1, 4) / n,
        "remaining_nonnegative_X": 1 / (2 * n),
    }
    return {
        "family_order": "every even integer n>=1024",
        "family_action_normalization": n / analytic.FIXED_GAMMA,
        "localized_scalar_affine_order_coefficient": linear,
        "localized_scalar_order_independent_coefficient": constant,
        "linear_and_constant_envelopes": {"linear": Glin, "constant": Gconst},
        "uniform_scalar_error_at_first_family_order": scalar_upper,
        "successive_even_order_majorant_ratio_upper": ratio,
        "uniform_curvature_errors": analytic.coefficient_error_bounds()[
            "h_weighted_curvature_and_DHOST_jet_error_bounds"
        ],
        "per_order_tensor_lower_bounds": floors,
        "positive_domain": "all real u; -1/(4n)<X<6/5",
        "uniform_clock_tube_tensor_floor": sp.Rational(4, 5),
        "checks": {
            "scalar_affine_in_order_identity": sp.factor(
                analytic.vacuum_lower_function(n / analytic.FIXED_GAMMA, n)
                - Ftree
                - n * linear
                - constant
            ),
            "selected_family_normalization": (n / analytic.FIXED_GAMMA).subs(
                n, analytic.ORDER
            )
            - analytic.KAPPA,
            "localizer_fourth_derivative_exponent": sp.diff(-(u**4), u) + 4 * u**3,
        },
        "bounds": {
            "even_step_all_majorants_contract_below_one_twentieth": ratio
            < sp.Rational(1, 20),
            "first_uniform_scalar_majorant_below_one_e_minus_400": scalar_upper
            < sp.Rational(1, 10**400),
            "negative_region_lower_above_one_half": floors["small_negative_X"].subs(
                n, analytic.ORDER
            )
            > sp.Rational(1, 2),
            "small_positive_region_lower_above_three_quarters": floors[
                "small_nonnegative_X"
            ].subs(n, analytic.ORDER)
            > sp.Rational(3, 4),
            "upper_X_domain_keeps_absolute_one_minus_X_squared_below_one": abs(
                1 - sp.Rational(6, 5) ** 2
            )
            < 1,
            "minimal_order_exceeds_required_jet_order": analytic.ORDER
            > analytic.JET_ORDER + 2,
        },
    }
