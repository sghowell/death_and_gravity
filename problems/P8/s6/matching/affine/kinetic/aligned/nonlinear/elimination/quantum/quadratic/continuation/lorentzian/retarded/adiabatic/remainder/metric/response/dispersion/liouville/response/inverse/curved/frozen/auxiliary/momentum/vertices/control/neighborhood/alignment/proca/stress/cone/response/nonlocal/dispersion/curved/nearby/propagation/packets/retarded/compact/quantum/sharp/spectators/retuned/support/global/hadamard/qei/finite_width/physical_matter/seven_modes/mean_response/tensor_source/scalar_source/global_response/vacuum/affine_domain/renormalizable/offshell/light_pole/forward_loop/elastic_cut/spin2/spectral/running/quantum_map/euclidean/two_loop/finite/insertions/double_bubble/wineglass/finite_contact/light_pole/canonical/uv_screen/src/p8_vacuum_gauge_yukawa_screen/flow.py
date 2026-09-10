"""One-loop marginal fixed-flow screen, derived from matrices and tensors."""

from functools import cache

import sympy as sp

from . import model


@cache
def data():
    a, y, L, t, a0 = sp.symbols(
        "positive_gauge_squared positive_Yukawa quartic_L nonnegative_log_scale positive_initial_gauge_squared",
        positive=True,
    )
    b = sp.Rational(5, 3)
    cf = model.data()["fundamental_Casimir"]
    Y = sp.diag(y, y, y, -y, -y, -y)
    YY = Y * Y
    beta_matrix = (
        (YY * Y + Y * YY) / 2
        + 2 * Y * Y * Y
        + 2 * sp.trace(YY) * Y
        - 3 * a * (cf * Y + Y * cf)
    )
    beta_y = sp.factor(beta_matrix[0, 0])
    phi = sp.symbols("scalar_Phi", real=True)
    potential = L * phi**4 / 24
    beta_potential = (
        sp.diff(potential, phi, 2) ** 2 / 2
        + 2 * phi * sp.trace(YY) * sp.diff(potential, phi)
        - 2 * sp.trace(Y**4) * phi**4
    )
    beta_L = sp.factor(sp.diff(beta_potential, phi, 4))
    ry = sp.Rational(19, 45)
    x = sp.symbols("quartic_to_gauge_ratio", real=True)
    polynomial = 3 * x * x + (48 * ry + 2 * b) * x - 288 * ry * ry
    rx = (sp.sqrt(65985) - 177) / 45
    q = 16 * sp.pi**2
    running_a = a0 / (1 + 2 * b * a0 * t / q)
    running_y2 = ry * running_a
    running_L = rx * running_a
    checks = {
        "gauge_flavor_coefficient": sp.Rational(11, 3) * 3 - sp.Rational(2, 3) * 14 - b,
        "matrix_Yukawa_beta": beta_matrix - (15 * y * y - 8 * a) * Y,
        "direct_quartic_beta": beta_L - (3 * L * L + 48 * y * y * L - 288 * y**4),
        "quartic_dictionary_to_L_over_24": 6
        * (18 * (L / 6) ** 2 + 48 * y * y * (L / 6) - 48 * y**4)
        - beta_L,
        "Yukawa_fixed_ratio": 15 * ry - 8 + b,
        "positive_quartic_root": sp.simplify(polynomial.subs(x, rx)),
        "integer_quartic_polynomial": sp.expand(
            225 * polynomial - (675 * x * x + 5310 * x - 11552)
        ),
        "running_gauge_ODE": sp.factor(
            q * sp.diff(running_a, t) + 2 * b * running_a**2
        ),
        "running_Yukawa_squared_ODE": sp.factor(
            q * sp.diff(running_y2, t)
            - 2 * running_y2 * (15 * running_y2 - 8 * running_a)
        ),
        "running_quartic_ODE": sp.simplify(
            q * sp.diff(running_L, t)
            - (3 * running_L**2 + 48 * running_y2 * running_L - 288 * running_y2**2)
        ),
        "gauge_initial_value": running_a.subs(t, 0) - a0,
        "Yukawa_initial_value": running_y2.subs(t, 0) - ry * a0,
        "quartic_initial_value": running_L.subs(t, 0) - rx * a0,
    }
    return {
        "loop_denominator": q,
        "beta_gauge_squared_times_loop_denominator": -2 * b * a * a,
        "beta_Yukawa_times_loop_denominator": beta_y,
        "beta_quartic_times_loop_denominator": beta_L,
        "matrix_Yukawa_derivation": beta_matrix,
        "potential_tensor_derivation": beta_potential,
        "Yukawa_squared_to_gauge_squared": ry,
        "quartic_to_gauge_squared": rx,
        "quartic_ratio_polynomial": polynomial,
        "quartic_root_lower_test": polynomial.subs(x, 1),
        "quartic_root_upper_test": polynomial.subs(x, 2),
        "quartic_polynomial_derivative": sp.diff(polynomial, x),
        "one_loop_running": {
            "gauge_squared": running_a,
            "Yukawa_squared": running_y2,
            "quartic": running_L,
        },
        "scope": "Exact solution of the displayed one-loop marginal equations for nonnegative log scale. The positive ray and its decay do not establish relevant-coupling stability, a nonperturbative UV theory, a mass gap or P8 positivity applicability. The superrenormalizable H Phi^2 coupling does not enter these one-loop marginal logarithms.",
        "checks": checks,
    }
