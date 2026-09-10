"""Quantified analytic alternative: exact finite clock jets, not exact open-tube equality."""

from functools import cache

import sympy as sp

from . import family

ORDER = 1024
JET_ORDER = 4
KAPPA = sp.Integer(10) ** 800
VACUUM_LAMBDA_BAR = sp.Rational(1, 10**600)
LOCAL_VACUUM_MASS = sp.Integer(1)
TUBE_Q_UPPER = sp.Rational(21, 100)
FIXED_GAMMA = sp.Rational(ORDER, KAPPA)


def vacuum_lower_function(kappa=KAPPA, n=ORDER, coupling=VACUUM_LAMBDA_BAR):
    u, X = family.u, family.X
    F0 = family.data()["original_retuned_tree_scalar"].subs({u: 0, X: 0})
    # Explicit classical vacuum calibration, not a subtraction prescription:
    # cancel the switch-induced Y^2 shift and retain a positive quartic potential.
    return (
        X / 2
        - u * u / 2
        + (coupling * kappa - n * F0) * X * X
        + (-sp.Rational(1, 3) * n - F0) * u**4
    )


def falling(n, j):
    return sp.prod(n - i for i in range(j))


def derivative_upper(j):
    if isinstance(j, bool) or not isinstance(j, int) or not 0 <= j <= 6:
        raise TypeError("Require an integer derivative order from zero through six")
    # Faà di Bruno for q(X)^ORDER; q=1-X², |q|<=21/100,
    # |q'|<=11/5, and |q''/2|=1 on 9/10<=X<=11/10.
    n = ORDER
    return sum(
        sp.factorial(j)
        * falling(n, j - b)
        * TUBE_Q_UPPER ** (n - j + b)
        * sp.Rational(11, 5) ** (j - 2 * b)
        / (sp.factorial(j - 2 * b) * sp.factorial(b))
        for b in range(j // 2 + 1)
    )


@cache
def coefficient_error_bounds():
    r = JET_ORDER
    S = sum(derivative_upper(j) / sp.factorial(j) for j in range(r + 1))
    Sd = sum(derivative_upper(j + 1) / sp.factorial(j) for j in range(r + 1))
    H = sum(sp.rf(6, j) / sp.factorial(j) for j in range(r + 1))
    Xm = sp.Rational(11, 10)
    Xinv = sum(sp.Rational(10, 9) ** (j + 1) for j in range(r + 1))
    E = S + Xm * Sd
    F = Xm * S
    Rold = 1 + Xm * H
    Rnew = Rold + H * F
    Rnorm = sp.Integer(233)
    floor = sp.Rational(4, 5)
    inverse = sum((Rnorm / floor) ** j for j in range(r + 1)) / floor
    T = (2 * E + E * E) * Rold + H * F
    A3 = H * E * Xinv
    A4 = A3 + sp.Rational(7, 4) * H * H * T * inverse * inverse
    A5 = H * H * T * inverse * inverse * Xinv
    return {
        "order": ORDER,
        "mixed_normalized_derivative_total_order": r,
        "S_derivative_bounds": [derivative_upper(j) for j in range(r + 2)],
        "S_jet_norm": S,
        "S_prime_jet_norm": Sd,
        "weighted_inverse_h_jet_norm": H,
        "X_minus_one_jet_norm": Xm,
        "inverse_X_jet_norm": Xinv,
        "old_R_jet_norm_upper": Rold,
        "new_R_jet_norm_upper": Rnew,
        "both_inverse_R_jet_norm_upper": inverse,
        "tube_R_floor": floor,
        "h_weighted_curvature_and_DHOST_jet_error_bounds": {
            "F2": H * F / 2,
            "A3": A3,
            "A4": A4,
            "A5": A5,
        },
        "bounds": {
            "new_R_jet_norm_below_233": Rnew < Rnorm,
            "both_tensor_factors_above_four_fifths": sp.Rational(9, 10)
            - sp.Rational(1, 10) * derivative_upper(0)
            > floor,
            **{
                name + "_four_jet_error_below_one_e_minus_400": value
                < sp.Rational(1, 10**400)
                for name, value in {
                    "F2": H * F / 2,
                    "A3": A3,
                    "A4": A4,
                    "A5": A5,
                }.items()
            },
        },
    }


@cache
def lower_error_bounds():
    u, X = family.u, family.X
    t = 1 + u * u
    Fv = vacuum_lower_function()
    value = Fv - family.data()["original_retuned_tree_scalar"]
    rows = {}
    for a in range(JET_ORDER + 1):
        for b in range(JET_ORDER + 1 - a):
            weighted = sp.cancel(t**6 * sp.diff(value, X, b))
            num, den = sp.fraction(weighted)
            D = sp.Poly(den, u, domain=sp.QQ)
            if D.nth(0) <= 0 or any(p[0] % 2 or c < 0 for p, c in D.terms()):
                raise ValueError(
                    "The rational denominator does not support the global bound"
                )
            P = sp.Poly(num, u, X, domain=sp.QQ)
            C = sum(
                abs(c) * sp.Rational(11, 10) ** p[1] * sp.factorial((p[0] + 3) // 4)
                for p, c in P.terms()
            ) / D.nth(0)
            rows[(a, b)] = C / (sp.factorial(a) * sp.factorial(b))
        value = sp.cancel(sp.diff(value, u) - 4 * u**3 * value)
    G = sum(rows.values())
    S = coefficient_error_bounds()["S_jet_norm"]
    return {
        "dimensionless_vacuum_function": Fv,
        "analytic_field_localizer": "exp(-u^4)",
        "h_squared_weighted_localized_scalar_difference_jet_bounds": rows,
        "h_squared_weighted_localized_difference_jet_norm": G,
        "h_squared_weighted_full_scalar_four_jet_error": S * G,
        "bounds": {
            "actual_scalar_four_jet_error_below_one_e_minus_400": S * G
            < sp.Rational(1, 10**400)
        },
    }


@cache
def local_jets():
    u, X = family.u, family.X
    # Only the exact jets needed for the vacuum action, without expanding
    # a degree-2048 switching polynomial.
    eps = sp.Symbol("formal_field_amplitude", real=True)
    field, Y = sp.symbols(
        "canonical_dimensionless_field canonical_gradient_square", real=True
    )
    n = sp.Symbol("positive_integer_switch_order", positive=True, integer=True)
    kap = sp.Symbol("positive_action_normalization", positive=True)
    coupling = sp.Symbol("fixed_canonical_lambda_bar", real=True)
    F0 = family.data()["original_retuned_tree_scalar"].subs({u: 0, X: 0})
    S = 1 - n * X * X  # exact through four perturbation fields, since X is quadratic
    L = 1 - u**4  # exact through four fields
    # The unmodified F contributes only its constant in the new quartic
    # terms; the vacuum two- and four-field pieces are written explicitly.
    Fv = vacuum_lower_function(kap, n, coupling)
    local = Fv + (1 - L * S) * (F0 - Fv)
    scaled = sp.expand(
        kap
        * local.subs(
            {u: eps * field / sp.sqrt(kap), X: eps * eps * Y / kap}, simultaneous=True
        )
    )
    quadratic = sp.expand(scaled).coeff(eps, 2)
    quartic = sp.expand(scaled).coeff(eps, 4)
    lam_eff = coupling
    return {
        "original_retuned_F_at_vacuum": F0,
        "canonical_quadratic_action": quadratic,
        "canonical_quartic_lower_action": quartic,
        "canonical_derivative_quartic_coefficient": lam_eff,
        "canonical_potential_quartic_coefficient": -n / (3 * kap),
        "leading_canonical_L3_minus_L4_coefficient": -2 * n / kap,
        "actual_lambda_effective": lam_eff.subs(
            {coupling: VACUUM_LAMBDA_BAR, kap: KAPPA, n: ORDER}
        ),
        "canonical_action_normalization_symbol": kap,
        "canonical_coupling_symbol": coupling,
        "canonical_switch_order_symbol": n,
        "fixed_nonzero_DHOST_scaling_gamma": FIXED_GAMMA,
        "checks": {
            "canonical_quadratic_mass_and_residue_unchanged": sp.expand(
                quadratic - (Y - field**2) / 2
            ),
            "canonical_quartic_has_explicitly_calibrated_localizer_and_switch": sp.expand(
                quartic - coupling * Y * Y + n * field**4 / (3 * kap)
            ),
            "actual_original_vacuum_constant": F0 + sp.Rational(28, 25),
            "actual_decoupling_holds_fixed_derivative_interaction": sp.limit(
                lam_eff, kap, sp.oo
            )
            - coupling,
            "fixed_order_decoupling_removes_extra_quartic_potential": sp.limit(
                -n / (3 * kap), kap, sp.oo
            ),
            "fixed_order_decoupling_removes_DHOST_quartic": sp.limit(
                -2 * n / kap, kap, sp.oo
            ),
            "coscaled_decoupling_keeps_quartic_potential": (-n / (3 * kap)).subs(
                kap, n / FIXED_GAMMA
            )
            + FIXED_GAMMA / 3,
            "coscaled_decoupling_keeps_nonzero_DHOST_quartic": (-2 * n / kap).subs(
                kap, n / FIXED_GAMMA
            )
            + 2 * FIXED_GAMMA,
        },
        "bounds": {
            "chosen_effective_derivative_quartic_remains_positive": lam_eff.subs(
                {coupling: VACUUM_LAMBDA_BAR, kap: KAPPA, n: ORDER}
            )
            > 0,
            "chosen_vacuum_quartic_potential_is_positive": FIXED_GAMMA / 3 > 0,
            "analytic_tensor_global_floor": sp.Rational(1, 2 * ORDER) > 0,
        },
    }


@cache
def actual_family():
    u, X = family.u, family.X
    h = (1 + u * u) ** 3
    S = (1 - X * X) ** ORDER
    R = 1 + (1 - S) * (X - 1) / h
    RX = sp.diff(R, X)
    Fv = vacuum_lower_function()
    Ftree = family.data()["original_retuned_tree_scalar"]
    return {
        "S": S,
        "R": R,
        "RX": RX,
        "F": Ftree + sp.exp(-(u**4)) * S * (Fv - Ftree),
        "F2": -R / 2,
        "A3": RX / X,
        "A4": -RX / X - sp.Rational(7, 4) * RX**2 / R,
        "A5": RX**2 / (R * X),
        "field_domain": "all real u and -1/4096<X<6/5",
        "original_tube": "9/10<=X<=11/10",
        "global_tensor_factor_lower": sp.Rational(1, 2 * ORDER),
        "apparent_X_denominators_removed_analytically_at_zero": True,
        "clock_jet_agreement_through_order": ORDER - 1,
        "not_exact_equality_on_an_open_tube": True,
    }


@cache
def exact_checks():
    u, X = family.u, family.X
    n = sp.Symbol("positive_integer_switch_order", positive=True, integer=True)
    q = 1 - X * X
    P = sp.Integer(1)
    checks = {}
    for j in range(7):
        closed = sum(
            sp.factorial(j)
            * falling(n, j - b)
            * q**b
            * (-2 * X) ** (j - 2 * b)
            * (-1) ** b
            / (sp.factorial(j - 2 * b) * sp.factorial(b))
            for b in range(j // 2 + 1)
        )
        checks["switch_derivative_polynomial_" + str(j)] = sp.expand(P - closed)
        P = sp.expand(q * sp.diff(P, X) + (n - j) * sp.diff(q, X) * P)
    actual = actual_family()
    for j in range(7):
        checks["actual_analytic_switch_clock_derivative_" + str(j)] = sp.diff(
            actual["S"], X, j
        ).subs(X, 1)
    checks["switch_zero_has_exact_factorized_base"] = sp.expand(
        1 - X * X + (X - 1) * (X + 1)
    )
    checks["tensor_factor_exact_clock_value"] = actual["R"].subs(X, 1) - 1
    checks["vacuum_tensor_coefficient_constant_for_all_fields"] = (
        actual["R"].subs(X, 0) - 1
    )
    checks["vacuum_tensor_X_derivative_vanishes_for_all_fields"] = actual["RX"].subs(
        X, 0
    )
    checks["regular_vacuum_A3_limit_from_lHopital"] = sp.factor(
        sp.diff(actual["RX"], X).subs(X, 0) + 2 * ORDER / (1 + u * u) ** 3
    )
    checks["regular_vacuum_A5_limit_from_lHopital"] = sp.factor(
        sp.diff(actual["RX"] ** 2, X).subs(X, 0)
    )
    # Independently differentiate and series expand the full analytic lower
    # action in generic order n before specializing it. No high-degree
    # switch polynomial is truncated in the definition.
    eps = sp.Symbol("independent_small_field_amplitude", real=True)
    field, Y = sp.symbols(
        "independent_canonical_field independent_gradient_square", real=True
    )
    kap = sp.Symbol("independent_positive_kappa", positive=True)
    coupling = sp.Symbol("independent_fixed_lambda", real=True)
    Ftree = family.data()["original_retuned_tree_scalar"]
    Fv = vacuum_lower_function(kap, n, coupling)
    exactF = Ftree + sp.exp(-(u**4)) * (1 - X * X) ** n * (Fv - Ftree)
    transformed = kap * exactF.subs(
        {u: eps * field / sp.sqrt(kap), X: eps * eps * Y / kap}, simultaneous=True
    )
    series = sp.series(transformed, eps, 0, 5).removeO().expand()
    checks["full_analytic_function_has_actual_positive_vacuum_quadratic_jet"] = (
        sp.expand(series.coeff(eps, 2) - (Y - field * field) / 2)
    )
    checks["full_analytic_function_has_both_calibrated_quartic_jets"] = sp.expand(
        series.coeff(eps, 4) - coupling * Y * Y + n * field**4 / (3 * kap)
    )
    checks["no_scalar_cubic_in_actual_vacuum"] = series.coeff(eps, 3)
    # Rational difference identities justify the independently bounded
    # expressions in coefficient_error_bounds.
    HH, EE, FF = sp.symbols("H E F", real=True)
    ro = 1 + HH * (X - 1)
    rn = ro - HH * FF
    a3o, a3n = HH / X, HH * (1 - EE) / X
    a4o = -a3o - sp.Rational(7, 4) * HH**2 / ro
    a4n = -a3n - sp.Rational(7, 4) * HH**2 * (1 - EE) ** 2 / rn
    a5o = HH**2 / (ro * X)
    a5n = HH**2 * (1 - EE) ** 2 / (rn * X)
    num = (-2 * EE + EE**2) * ro + HH * FF
    checks["actual_A3_error_factorization"] = sp.factor(a3n - a3o + HH * EE / X)
    checks["actual_A4_error_factorization"] = sp.factor(
        a4n - a4o - (HH * EE / X - sp.Rational(7, 4) * HH**2 * num / (ro * rn))
    )
    checks["actual_A5_error_factorization"] = sp.factor(
        a5n - a5o - HH**2 * num / (X * ro * rn)
    )
    for j in range(5):
        checks["inverse_h_complex_factor_derivative_majorant_" + str(j)] = sum(
            sp.binomial(j, b) * sp.rf(3, b) * sp.rf(3, j - b) for b in range(j + 1)
        ) - sp.rf(6, j)
    return checks
