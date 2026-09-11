"""A common complex field-amplitude circle for the full analytic target."""

from functools import cache

import sympy as s
from p8_exceptional_vacuum import analytic, family


@cache
def data():
    n, kappa, lam = analytic.ORDER, analytic.KAPPA, analytic.VACUUM_LAMBDA_BAR
    rho = s.Integer(10) ** 300
    U = rho / s.sqrt(kappa)
    Xcap = 4 * rho * rho / kappa
    u, X = family.u, family.X
    tree = family.data()["original_retuned_tree_scalar"]
    num, den = s.fraction(s.factor(tree))
    numerator_l1 = sum(abs(c) for c in s.Poly(num, u, X).coeffs())
    treeup = numerator_l1 / (200 * s.Rational(15, 16) ** 12)
    F0 = abs(tree.subs({u: 0, X: 0}))
    Vup = (
        (Xcap + U * U) / 2
        + (lam * kappa + n * F0) * Xcap**2
        + (s.Rational(n, 3) + F0) * U**4
    )
    Tup = s.Rational(2, 3**n)
    Wup = 2 * n * Xcap**2
    Sup = (1 + Tup) * (1 + Wup)
    Fup = treeup + 3 * Sup * (Vup + treeup)
    Bup = Tup + (1 + Tup) * Wup
    Rdef = (1 + Xcap) * Bup / (1 - U * U) ** 3
    TXup = s.Rational(8, 3**n)
    TpXup = s.Rational(100 * n, 3**n)
    WXup = 2 * n * Xcap
    WpXup = 8 * n
    BXup = TXup + (1 + Tup) * WXup
    BpXup = TpXup * (1 + Wup) + (1 + Tup) * WpXup
    A3up = ((1 + Xcap) * BpXup + BXup) / (1 - U * U) ** 3
    RXup = Xcap * A3up
    A4up = A3up + s.Rational(7, 4) * RXup**2 / (1 - Rdef)
    A5up = A3up**2 * Xcap / (1 - Rdef)
    density = (
        kappa * Fup
        + 64 * rho**4 * (A3up + A4up) / kappa
        + 256 * rho**6 * A5up / kappa**2
    )
    M = s.Integer(10) ** 811
    return {
        "fixed_n": n,
        "fixed_kappa": kappa,
        "fixed_lambda": lam,
        "complex_field_amplitude_radius": rho,
        "dimensionless_field_modulus_upper": U,
        "dimensionless_gradient_invariant_modulus_upper": Xcap,
        "retuned_scalar_rational_coefficient_bound": treeup,
        "full_analytic_scalar_coefficient_upper": Fup,
        "R_distance_from_one_upper": Rdef,
        "A3_absolute_upper": A3up,
        "A4_absolute_upper": A4up,
        "A5_absolute_upper": A5up,
        "full_canonical_flat_density_circle_upper": M,
        "bounds": {
            "dimensionless_complex_field_and_X_inside_quarter": bool(
                U < s.Rational(1, 4) and Xcap < s.Rational(1, 4)
            ),
            "rational_step_denominator_has_strict_gap": bool(
                Tup < s.Rational(1, 10**450)
            ),
            "Gaussian_bump_exponential_argument_below_one_half": bool(
                n * Xcap**2 < s.Rational(1, 2)
            ),
            "retuned_tree_bound_below_one_e_8": bool(treeup < 10**8),
            "calibrated_vacuum_lower_function_modulus_below_one": bool(Vup < 1),
            "full_switch_modulus_below_four": bool(Sup < 4),
            "full_scalar_coefficient_modulus_below_one_e_10": bool(Fup < 10**10),
            "full_R_holomorphic_inverse_gap": bool(Rdef < s.Rational(1, 2)),
            "full_A3_modulus_below_one_e_6": bool(A3up < 10**6),
            "full_A4_modulus_below_two_e_6": bool(A4up < 2 * 10**6),
            "full_A5_modulus_below_eight_e_minus_188": bool(
                A5up < s.Rational(8, 10**188)
            ),
            "full_density_modulus_below_declared_circle_upper": bool(density < M),
        },
        "checks": {
            "literal_retuned_rational_denominator": s.expand(
                den - 200 * (1 + u * u) ** 12
            ),
            "same_actual_canonical_normalization": kappa - 10**800,
            "same_actual_switch_order": n - 1024,
            "same_actual_fixed_gamma": n / kappa - analytic.FIXED_GAMMA,
            "complex_field_scaling": U * s.sqrt(kappa) - rho,
            "four_gradient_components_in_X_bound": Xcap * kappa - 4 * rho * rho,
            "full_density_all_scalar_and_DHOST_groups": density
            - kappa * Fup
            - 64 * rho**4 * (A3up + A4up) / kappa
            - 256 * rho**6 * A5up / kappa**2,
        },
        "scope": "This is a complex field-amplitude Cauchy circle at fixed jets and fixed flat metric, not a momentum cutoff, quantum loop bound, global field-map inverse or rolling-background domain. The rational step is kept in every full-function modulus bound.",
    }
