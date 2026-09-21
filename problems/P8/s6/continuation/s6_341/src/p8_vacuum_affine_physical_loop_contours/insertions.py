"""Literal minimal TT line insertion and explicit weighted absolute budgets."""

from functools import cache

import sympy as s
from p8_vacuum_affine_massive_scalar_graviton_bremsstrahlung import tree

from . import parameters, source

eta = s.diag(1, -1, -1, -1)


def dot(a, b):
    return (a.T * eta * b)[0]


def avg(expr, loop):
    expr = s.Poly(s.expand(expr), *loop)
    out = 0
    sigma = s.Symbol("isotropic_radial_second_moment")
    for powers, coeff in expr.terms():
        degree = sum(powers)
        if degree == 0:
            out += coeff
        elif degree == 2 and max(powers) == 2:
            out += coeff * sigma * eta[powers.index(2), powers.index(2)]
        elif degree > 2:
            raise AssertionError("unexpected loop numerator degree")
    return s.expand(out)


BUDGETS = {
    "triangle_light": s.Integer(8388608),
    "triangle_heavy": s.Integer(65536000),
    "box_light": s.Integer(1073741824),
    "box_heavy": s.Integer(16777216000),
}


def line_bound(kind, line_type, mass=source.HEAVY_MASS2, numerator_norm2=10000):
    kind = parameters.require_kind(kind)
    mass = parameters.heavy_mass(mass)
    if not isinstance(line_type, str) or line_type not in ("light", "heavy"):
        raise ValueError("Require a literal light or heavy internal line")
    bound = parameters.rational(numerator_norm2)
    if not 0 <= bound <= 10000:
        raise ValueError("Require the proved routing-vector squared-norm budget")
    topology = "triangle" if kind == "triangle" else "box"
    power = (1 if kind == "triangle" else 2) + (line_type == "heavy")
    return BUDGETS[topology + "_" + line_type] * bound / mass**power


@cache
def data():
    checks = {}
    gates = {}
    R = s.Rational
    k = s.Matrix([1, 0, 0, 1])
    Qbar = s.Matrix(s.symbols("barQ0:4", real=True))
    Qa = s.Matrix(s.symbols("lineQ0:4", real=True))
    mass, y, xa, gamma = s.symbols("mass y xa gamma", real=True)
    loop = s.Matrix(s.symbols("loop0:4", real=True))
    Rvec = Qa - Qbar
    shifted = loop + Rvec - y * k
    literal = tree.tensor(shifted, shifted + k, mass, eta)
    pols = [s.diag(0, 1, -1, 0), s.zeros(4)]
    pols[1][1, 2] = pols[1][2, 1] = 1
    for i, eps in enumerate(pols):
        full = sum(eps[a, b] * literal[a, b] for a in range(4) for b in range(4))
        checks[f"literal_TT_vertex_isotropic_average_{i}"] = s.expand(
            avg(full, loop) - 2 * (Rvec.T * eps * Rvec)[0]
        )
        checks[f"effective_line_shift_TT_invariance_{i}"] = s.expand(
            ((Rvec + (gamma - y) * k).T * eps * (Rvec + (gamma - y) * k))[0]
            - (Rvec.T * eps * Rvec)[0]
        )
    checks["general_split_denominator_change"] = s.expand(
        -2 * y * dot(Qa, k)
        + dot(Qbar + y * k, Qbar + y * k)
        - dot(Qbar, Qbar)
        - 2 * y * dot(Qbar - Qa, k)
    )
    checks["effective_single_line_denominator"] = s.expand(
        -xa * (dot(Qa + gamma * k, Qa + gamma * k) - dot(Qa, Qa))
        + dot(Qbar + xa * gamma * k, Qbar + xa * gamma * k)
        - dot(Qbar, Qbar)
        - 2 * xa * gamma * dot(Qbar - Qa, k)
    )
    checks["split_parameter_Jacobian"] = (
        s.Matrix([xa - y, y]).jacobian([xa, y]).det() - 1
    )
    checks["unit_gamma_interval_weight"] = s.diff(xa * gamma, gamma) - xa
    for N in (3, 4):
        eps = s.Symbol("dimensional_epsilon")
        D = 4 - 2 * eps
        checks[f"N{N}_full_radial_Gamma_argument"] = s.expand(
            N + 1 - D / 2 - (N - 1 + eps)
        )
        checks[f"N{N}_D4_radial_Gamma_factor"] = s.gamma(N - 1) - s.factorial(N - 2)
        aa, dd, gg = s.symbols("a Delta gamma")
        p = N - 2
        # Inner split integration is a divided difference, including coincidence.
        anti = -((dd + 2 * aa * gg) ** (-p)) / (2 * aa * p)
        checks[f"N{N}_split_exact_antiderivative"] = s.factor(
            s.diff(anti, gg) - (dd + 2 * aa * gg) ** (-(p + 1))
        )
        checks[f"N{N}_split_coincident_limit"] = s.limit(
            (anti.subs(gg, 1) - anti.subs(gg, 0)), aa, 0
        ) - dd ** (-(p + 1))

    # Scalar contour bounds also cover the entire subthreshold interval[-12,4/3].
    gates["subthreshold_light_gap"] = 1 - R(4, 3) / 4 == R(2, 3) > R(1, 32)
    gates["subthreshold_light_upper"] = 1 + 12 / 4 == 4
    gates["complex_simplex_total_absolute_weight"] = 4 * (
        1 - s.Symbol("t")
    ) + 2 * s.Symbol("t") == 4 - 2 * s.Symbol("t")
    gates["routing_vector_budget"] = (12 + 4 * 12) ** 2 < 10000
    budgets = {
        "triangle_light": 8 * 4 * 2 * 64**2 * 32,
        "triangle_heavy": 8 * 2 * 2 * 64**2 * 500,
        "box_light": 16 * 4 * 2 * 2 * 64**3 * 16,
        "box_heavy": 16 * 2 * 2 * 2 * 64**3 * 500,
    }
    for key, expected in (
        ("triangle_light", 8388608),
        ("triangle_heavy", 65536000),
        ("box_light", 1073741824),
        ("box_heavy", 16777216000),
    ):
        checks[key + "_weighted_absolute_budget"] = s.Integer(budgets[key]) - expected
    gates["all_triangle_lines_with_R2_budget"] = (
        2 * budgets["triangle_light"] + R(budgets["triangle_heavy"], 10**6)
    ) * 10000 < 10**12
    gates["all_box_lines_with_R2_budget"] = (
        2 * budgets["box_light"] + R(2 * budgets["box_heavy"], 10**6)
    ) * 10000 < 10**14

    # Two required deformations: each omitted-variable control really fails.
    z = R(1, 10**6 + 2)
    n = (1 - z) ** 2 / z
    checks["omitted_heavy_contour_actual_pole"] = -((1 - z) ** 2) + n * z
    gates["omitted_heavy_contour_mass_in_theorem_domain"] = 10**6 <= n < 10**198
    checks["omitted_light_contour_endpoint_pole"] = 1 - R(25, 4) * R(1, 5) * R(4, 5)
    xx = R(1, 5)
    tt = R(1, 3)
    X = xx - s.I * xx * (1 - xx) * (1 - 2 * xx)
    Z = tt + s.I * tt * (1 - tt)
    bad = (1 - Z) ** 2 * (1 - R(25, 4) * X * (1 - X)) + 10**6 * Z
    gates["wrong_sign_contour_violates_Feynman_prescription"] = s.im(bad) > 0

    affineP = s.Matrix(s.symbols("affineP0:4", real=True))
    checks["null_shift_invariant_affine_in_gamma"] = s.expand(
        dot(affineP + gamma * k, affineP + gamma * k)
        - (1 - gamma) * dot(affineP, affineP)
        - gamma * dot(affineP + k, affineP + k)
    )
    return {
        "checks": {key: s.factor(value) for key, value in checks.items()},
        "gates": {key: bool(value) for key, value in gates.items()},
        "whole_literal_TT_numerator": 2 * Rvec * Rvec.T,
        "whole_split_denominator": "Delta_y=Delta+2y(barQ-Qa).k; y=xa*gamma. This equals the original N-line denominator with Qa replaced by Qa+gamma*k, since k^2=0.",
        "whole_full_TT_line_integral": "The literal scalar stress tensor, isotropically averaged in exact D, gives2*epsilon(Qa-barQ,Qa-barQ), independent of y. The loop measure is xa*d_gamma times the original simplex measure, with Gamma(N-1+epsilon)/Delta_gamma^(N-1+epsilon) for D=4-2epsilon. At D4 the factors are1 for triangles and2 for boxes. Overall Wick/vertex signs are not needed for this absolute bound and are not used to assemble a complete amplitude.",
        "whole_weighted_per_line_budgets": BUDGETS,
        "whole_routing_budget": "Original cyclic Qj are partial sums of hard vectors with Euclidean norm<=12. On the complex simplex sum|xj|<=4, so|Qa-barQ|<=60 and its squared norm<10000. Physical canonical TT Frobenius projection is contractive.",
        "whole_complete_internal_triangle_bound": "Sum of its two light and one heavy insertions has norm<10^12/n, excluding couplings, the common loop1/(16pi^2), and1/sqrt(kappa).",
        "whole_complete_internal_box_bound": "Sum of its two light and two heavy insertions has norm<10^14/n^2, with the same excluded factors. Ordered mass assignments are never exchanged.",
        "whole_scope": "These are line-insertion kernels and absolute bounds, not the multiplicity/sign/source/outer-heavy-branch/finite-counterterm completed four-hard radiative amplitude. No Ward-only transverse matching inference.",
    }
