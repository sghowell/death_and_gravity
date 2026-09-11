"""Full affine trace source, regular lower dictionary, and exact new alignment."""

from functools import cache

import sympy as s
from p8_affine import connection
from p8_affine_retuned import degeneracy


@cache
def data():
    p, v = connection.P, connection.S
    R, RX, Ru, q, H = s.symbols("R R_X R_u q_regular H_clock", real=True)
    X = s.symbols("positive_X", positive=True)
    delta = R - 1
    old = degeneracy.source()
    substitutions = {
        connection.F3: (q + Ru / 2) / (4 * p * X),
        connection.PP: Ru / (8 * p),
        connection.PX: -RX / (8 * p),
        v**2: X,
    }
    expected = {
        "alpha": s.Rational(3, 2) * (q + Ru / 2) + 3 * delta * Ru / (4 * R),
        "c": -delta / X,
        "d": -delta / X**2 + 3 * delta * RX / (2 * R * X),
    }
    checks = {}
    for key, value in expected.items():
        actual = old[key].subs(substitutions, simultaneous=True)
        checks["literal_all64_" + key] = s.factor(
            s.expand(actual).subs(p**2, R / 4) - value
        )
    B = 3 * H * delta + s.Rational(3, 2) * (q + Ru / 2)
    Box, Z = s.symbols("Box_P8_u uHu_P8", real=True)
    source_scalar = delta * (
        -3 * H + 3 * Ru / (4 * R) + Box / X + (-1 / X**2 + 3 * RX / (2 * R * X)) * Z
    )
    H00, trace, K = s.symbols("H00 spatial_Hessian_trace Khat", real=True)
    original_scalar = (
        expected["alpha"] + expected["c"] * (-H00 + trace) + expected["d"] * X * H00 - B
    )
    transformed_trace = (
        -v * K + 3 * v * v * Ru / (4 * R) + 3 * v * v * RX * H00 / (2 * R)
    )
    actual_normal = s.factor(
        v * original_scalar.subs(trace, transformed_trace).subs(X, v * v)
    )
    desired = (R - 1) * (K - 3 * H * v)
    checks["full_covariant_source_to_clock_chart"] = s.factor(actual_normal - desired)
    checks["source_has_no_lapse_velocity"] = s.diff(actual_normal, H00)
    checks["regular_lower_q_cancels_from_shifted_source"] = s.diff(source_scalar, q)
    checks["source_signature_dictionary"] = s.factor(
        expected["alpha"] - B - expected["c"] * Box + expected["d"] * Z - source_scalar
    )
    eps, n, k1, h = s.symbols("epsilon lapse_perturbation Khat_first h", real=True)
    vjet = 1 - eps * n
    deltajet = (vjet * vjet - 1) / h
    expression = s.expand(deltajet * (3 * H + eps * k1 - 3 * H * vjet))
    second = -2 * n * (k1 + 3 * H * n) / h
    checks["clock_source_zero"] = expression.coeff(eps, 0)
    checks["clock_source_first_zero"] = expression.coeff(eps, 1)
    checks["clock_source_second_nontrivial"] = s.factor(
        expression.coeff(eps, 2) - second
    )
    a, b = s.symbols("a b", real=True)
    x = s.symbols("x", real=True)
    germ_delta = a * x * x + b * x**3
    regular_box = s.cancel(germ_delta / x)
    regular_Z = s.cancel(-germ_delta / x**2)
    checks["regular_null_gradient_box_coefficient"] = regular_box.subs(x, 0)
    checks["regular_null_gradient_Hessian_coefficient"] = regular_Z.subs(x, 0) + a
    return {
        "unshifted_trace_coefficients_source_signature": expected,
        "new_local_B": B,
        "regular_lower_dictionary": "The unique S6.109 vacuum-regular analytic q is retained. Its nonzero clock boundary is not reset to the old singular choice.",
        "new_retained_vector": "W_mu=T_mu-B(u,X)*partial_mu u, T=V_trace+U_trace of the unrestricted projective quotient",
        "shifted_source_P8": source_scalar,
        "source_one_form": "S_mu=partial_mu u times shifted_source_P8; Box and uHu use the prescribed +--- physical metric",
        "normal_source_timelike_hat_chart": desired,
        "new_source_second_variation_coefficient": second,
        "regularity": "R-1 is divisible by X^2; all X denominators in S and B have analytic removable limits, including nonzero null gradients. S is generally nonzero away from the clock.",
        "checks": checks,
    }
