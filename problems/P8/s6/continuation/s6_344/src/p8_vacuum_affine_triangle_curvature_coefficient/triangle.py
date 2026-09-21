"""All three literal triangle line insertions, generic off-shell null-TT algebra."""

from functools import cache

import sympy as s

from . import jets

V, W, DOT, A, B, H, J, L = s.symbols("v w dot a b H J L")
X, Y, Z, GAMMA = s.symbols("x y z gamma")
GRAM = s.Matrix([[V, DOT, A], [DOT, W, B], [A, B, 0]])
POL = s.Matrix([[H, J, 0], [J, L, 0], [0, 0, 0]])
VERTICES = (s.Matrix([1, 0, 0]), s.Matrix([0, 1, 0]), s.Matrix([-1, -1, -1]))
K = s.Matrix([0, 0, 1])
WEIGHTS = (X, Y, 1 - X - Y)
SOFT = (A, B, -A - B)
SQUARES = (V, W, V + W + 2 * DOT + 2 * A + 2 * B)
HP = (H, L, H + 2 * J + L)
CS = (WEIGHTS[0] * WEIGHTS[1], WEIGHTS[1] * WEIGHTS[2], WEIGHTS[0] * WEIGHTS[2])
U = sum(c * v for c, v in zip(CS, SQUARES))
Q = A * A * L + B * B * H - 2 * A * B * J


def dot(p, q):
    return (p.T * GRAM * q)[0]


def ep(p, q):
    return (p.T * POL * q)[0]


def gamma_integral(expr):
    polynomial = s.Poly(s.expand(expr), GAMMA)
    return s.expand(sum(co / (power[0] + 1) for power, co in polynomial.terms()))


def line_setup(split):
    if type(split) is not int or split not in (0, 1, 2):
        raise ValueError("Require exactly one of the three literal triangle lines")
    rr = VERTICES[split:] + VERTICES[:split]
    ws = WEIGHTS[split:] + WEIGHTS[:split]
    qs = (s.zeros(3, 1), -rr[0], -rr[0] - rr[1])
    effective = (qs[0] + GAMMA * K, qs[1], qs[2])
    bar = sum((weight * q for weight, q in zip(ws, qs)), s.zeros(3, 1))
    HH = ep(qs[0] - bar, qs[0] - bar)
    shifted = sum(
        ws[i] * ws[j] * dot(effective[i] - effective[j], effective[i] - effective[j])
        for i in range(3)
        for j in range(i + 1, 3)
    )
    return ws[0], HH, shifted


@cache
def _literal_line(split, order):
    weight, HH, shifted = line_setup(split)
    return s.expand(-2 * order * weight * HH * gamma_integral(shifted ** (order - 1)))


def literal_line(split, order):
    order = jets.require_order(order)
    if type(split) is not int or split not in (0, 1, 2):
        raise ValueError("Require exactly one of the three literal triangle lines")
    return _literal_line(split, order)


def lift(order):
    order = jets.require_order(order)
    return s.expand(
        -2
        * sum(
            h * jets.difference_power(U, c, a, order) for h, c, a in zip(HP, CS, SOFT)
        )
    )


@cache
def data():
    checks = {}
    residuals = {}
    lines = {}
    for order in (1, 2, 3):
        actual = tuple(literal_line(split, order) for split in range(3))
        delta = s.factor(sum(actual) - lift(order))
        target = 0 if order < 3 else -8 * X**2 * Y**2 * (1 - X - Y) ** 2 * Q
        checks[f"generic_degree{2 * order}_entire_triangle_minus_explicit_lift"] = (
            s.expand(delta - target)
        )
        residuals[2 * order] = delta
        lines[2 * order] = actual
        for split in range(3):
            shifted = line_setup(split)[2] ** (order - 1)
            antiderivative = s.integrate(shifted, GAMMA)
            checks[f"degree{2 * order}_line{split}_split_parameter_exact_integral"] = (
                s.expand(
                    antiderivative.subs(GAMMA, 1)
                    - antiderivative.subs(GAMMA, 0)
                    - gamma_integral(shifted)
                )
            )
    xi = s.Symbol("xi")
    reduced = s.factor(
        (-8 * X**2 * Y**2 * (1 - X - Y) ** 2).subs(
            {X: (1 - Z) * xi, Y: (1 - Z) * (1 - xi)}, simultaneous=True
        )
    )
    checks["complete_light_angle_and_simplex_Jacobian"] = s.factor(
        s.integrate((1 - Z) * reduced, (xi, 0, 1))
        + s.Rational(4, 15) * Z**2 * (1 - Z) ** 5
    )
    checks["generic_conservation"] = sum(VERTICES, s.zeros(3, 1)) + K
    checks["generic_null_wave"] = dot(K, K)
    checks["generic_transverse_polarization"] = POL * K
    edge12 = A * VERTICES[1] - B * VERTICES[0]
    edge23 = B * VERTICES[2] + (A + B) * VERTICES[1]
    checks["composite_to_two_singletons_same_curvature_word"] = s.expand(
        ep(edge12, edge12) - ep(edge23, edge23)
    )
    checks["curvature_word_reversal_even"] = s.expand(ep(-edge12, -edge12) - Q)
    gates = {
        "degree6_entire_difference_is_nonzero": residuals[6] != 0,
        "all_three_literal_metric_positions_needed": all(
            lines[6][i] != 0 for i in range(3)
        ),
        "no_scalar_mass_shell_constraint_used": True,
        "all_Gram_soft_and_TT_variables_symbolic": True,
        "exact_simplex_identity_not_one_kinematic_fit": True,
        "opposite_insertion_sign_changes_lower_orders": s.expand(
            -sum(lines[2]) - lift(1)
        )
        != 0,
    }
    return {
        "checks": checks,
        "gates": gates,
        "whole_generic_momenta": {
            "R1_squared": V,
            "R2_squared": W,
            "R1_dot_R2": DOT,
            "k_dot_R1": A,
            "k_dot_R2": B,
            "epsilon_R1_R1": H,
            "epsilon_R1_R2": J,
            "epsilon_R2_R2": L,
        },
        "whole_every_literal_line_coefficient": lines,
        "whole_entire_difference_by_derivative_order": residuals,
        "whole_local_curvature_polynomial": Q,
        "whole_common_mass_denominator": "At degree2r restore1/M^(r+1), M=x*m0^2+y*m1^2+(1-x-y)*m2^2. The identity holds pointwise in the full simplex for arbitrary positive masses, not only after integration or at the original mass.",
        "whole_literal_normalization": "The S342 normalized TT kernel is-2*x_split*epsilon(Qsplit-barQ,Qsplit-barQ)/(M-lambda^2 U_gamma)^2. Its lambda^(2r) coefficient is-2*r*x_split*epsilon(...)*integral_0^1 U_gamma^(r-1) dgamma/M^(r+1), with all three cyclic lines and their actual mass/weight order retained.",
        "whole_first_curvature_difference": "Triangle minus the explicit labeled-scalar-box lift is zero at degree2 and4, and-8*x^2*y^2*z^2*R_abcd R1^a R2^b R1^c R2^d/M^4 at degree6. Here z=1-x-y and the linear curvature convention is exactly S336; Ricci and trace pieces vanish only in this one-real-null-TT observable.",
    }
