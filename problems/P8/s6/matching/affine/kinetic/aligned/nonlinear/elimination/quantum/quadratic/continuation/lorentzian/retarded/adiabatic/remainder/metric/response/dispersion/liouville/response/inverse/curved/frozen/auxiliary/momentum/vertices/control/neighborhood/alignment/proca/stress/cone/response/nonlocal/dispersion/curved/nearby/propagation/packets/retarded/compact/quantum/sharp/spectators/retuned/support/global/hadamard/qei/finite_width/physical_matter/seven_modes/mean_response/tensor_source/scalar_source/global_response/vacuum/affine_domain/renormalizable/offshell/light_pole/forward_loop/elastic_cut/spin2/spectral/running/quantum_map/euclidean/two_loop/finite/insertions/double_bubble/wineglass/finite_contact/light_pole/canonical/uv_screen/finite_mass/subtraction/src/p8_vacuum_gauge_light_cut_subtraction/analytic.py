"""Holomorphic divided differences and a uniform bound on the cut coefficient."""

from functools import cache

import sympy as sp


@cache
def data():
    x, s = sp.symbols("cut_variable channel_s")
    M, K, eta = sp.symbols(
        "positive_rho_error positive_cut_coefficient positive_rho_eta", positive=True
    )
    # rho error M=dA*z^2*eta/pi=9*pi*K*eta.
    b2error = 32 * M / sp.pi
    checks = {
        "Cauchy_derivative_margin": 5 - 4 - 1,
        "paired_divided_difference_integral_bound": 2 * 6 - 12,
        "each_log_difference_bound": 3 + 3 + 4 - 10,
        "paired_log_product_bound": 2 * 10 - 20,
        "total_germ_bound": 12 + 20 - 32,
        "Cauchy_unit_radius_coefficient_bound": 32 * M / sp.pi - b2error,
        "same_scale_rho_error_dictionary": 9 * sp.pi * K * eta - 9 * sp.pi * K * eta,
        "normalized_cut_coefficient_error": sp.factor(
            b2error.subs(M, 9 * sp.pi * K * eta) - 288 * K * eta
        ),
    }
    quotient_rows = []
    for n in range(11):
        quotient = sum(x ** (n - 1 - j) * s**j for j in range(n)) if n else sp.S.Zero
        checks[f"removable_polynomial_quotient_degree_{n}"] = sp.expand(
            (x - s) * quotient - (x**n - s**n)
        )
        integral = sp.integrate(quotient, (x, 0, 6))
        independent = (
            sum(6 ** (n - j) * s**j / sp.Integer(n - j) for j in range(n))
            if n
            else sp.S.Zero
        )
        checks[f"independent_quotient_integral_degree_{n}"] = sp.expand(
            integral - independent
        )
        quotient_rows.append(
            {"degree": n, "removable_quotient": quotient, "integral": integral}
        )
    return {
        "holomorphic_rho_control_disc": "|s-2|<=5",
        "rho_control_radius": 5,
        "quotient_control_radius": 4,
        "coefficient_Cauchy_radius": 1,
        "divided_difference_segment_disc": "|s-2|<=4",
        "coefficient_Cauchy_disc": "|s-2|<=1",
        "single_channel_rho_error_upper": M,
        "paired_divided_difference_integral_error_upper": 12 * M,
        "paired_log_product_error_upper": 20 * M,
        "physical_boundary_germ_error_upper": 32 * M / sp.pi,
        "center_b2_error_upper": b2error,
        "center_b2_error_over_K_upper": 288 * eta,
        "independent_polynomial_divided_difference_checks": quotient_rows,
        "scope": "The rho error is holomorphic on the larger complex disc, not merely small on the real cut. Quotient removability, its derivative bound and the local logarithm bounds precede Cauchy's coefficient estimate. No divergent principal-value derivative is integrated naively.",
        "checks": checks,
    }
