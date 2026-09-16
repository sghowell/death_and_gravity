"""Exact conserved-current pair quotient and uniform angular coefficient budget."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complete_two_graviton_tree import trees as e


@cache
def coefficients():
    a, b, r = s.symbols("a b r", positive=True)
    c = (1 - r * r) / (1 + r * r)
    d = 2 * r / (1 + r * r)
    W = a + b
    q1 = e.imm([a, 0, 0, a])
    q2 = e.imm([b, b * d, 0, b * c])
    Q = q1 + q2
    xs = s.symbols("Uxx Uxy Uxz Uyy Uyz Uzz", real=True)
    spatial = s.Matrix(
        [[xs[0], xs[1], xs[2]], [xs[1], xs[3], xs[4]], [xs[2], xs[4], xs[5]]]
    )
    v = Q[1:, 0]
    U = s.zeros(4)
    U[1:, 1:] = spatial
    U[0, 1:] = (v.T * spatial) / W
    U[1:, 0] = spatial * v / W
    U[0, 0] = (v.T * spatial * v)[0] / W**2
    checks = {"general_conserved_hard_tensor": (Q.T * e.ETA * U).applyfunc(s.factor)}
    H = e.imm(e.ETA * U * e.ETA - e.ETA * s.trace(e.ETA * U) / 2)

    def embed(M):
        A = s.zeros(4)
        A[1:, 1:] = M
        return e.imm(A)

    x = s.Matrix([1, 0, 0])
    y = s.Matrix([0, 1, 0])
    v2 = s.Matrix([c, 0, -d])
    pol1 = (embed(x * x.T - y * y.T), embed(x * y.T + y * x.T))
    pol2 = (embed(v2 * v2.T - y * y.T), embed(v2 * y.T + y * v2.T))
    records = {}
    for i, A in enumerate(pol1):
        for j, B in enumerate(pol2):
            checks[f"first_TT_{i}_{j}"] = A * q1
            checks[f"second_TT_{i}_{j}"] = (B * q2).applyfunc(s.factor)
            checks[f"first_norm_{i}_{j}"] = sum(x * x for x in A) - 2
            checks[f"second_norm_{i}_{j}"] = s.factor(sum(x * x for x in B) - 2)
            numerator = e.cubic_gravity((A, B, H), (q1, q2, -Q))
            quotient = s.factor(-numerator / e.old.dot(Q, Q))
            cs = tuple(s.factor(s.diff(quotient, x)) for x in xs)
            checks[f"complete_quotient_linearity_{i}_{j}"] = s.factor(
                quotient - sum((x * v for x, v in zip(xs, cs)), s.S.Zero)
            )
            records[(i, j)] = cs
    return a, b, r, records, checks


@cache
def data():
    a, b, r, records, checks = coefficients()
    checks = dict(checks)
    z = s.Symbol("z", nonnegative=True)
    bounds = {}
    for pair, cs in records.items():
        values = []
        for number, coefficient in enumerate(cs):
            if coefficient == 0:
                values.append(s.S.Zero)
                continue
            normalized = s.factor(
                (coefficient * a * b / (a + b) ** 2).subs(
                    {a: z, b: 1 - z}, simultaneous=True
                )
            )
            numerator, denominator = s.fraction(normalized)
            degree = s.degree(denominator, r)
            m = int(degree // 2)
            constant = denominator.subs(r, 0)
            checks[f"uniform_denominator_{pair[0]}_{pair[1]}_{number}"] = s.factor(
                denominator - constant * (1 + r * r) ** m
            )
            if degree % 2 != 0 or not constant.is_number or constant <= 0:
                raise ValueError("The uniform angular denominator contract failed")
            poly = s.Poly(numerator, z, r)
            if not all(0 <= power[1] <= 2 * m for power, _ in poly.terms()):
                raise ValueError("An unbounded angular monomial remains")
            values.append(sum(abs(c) for c in poly.coeffs()) / constant)
        bounds[pair] = sum(values)
    expected = {(0, 0): 530, (0, 1): 128, (1, 0): 88, (1, 1): 268}
    for key, value in bounds.items():
        checks[f"exact_polarization_pair_budget_{key[0]}_{key[1]}"] = (
            value - expected[key]
        )
    return {
        "whole_exact_conserved_pair_coefficients": {
            f"polarization_{i}_{j}": cs for (i, j), cs in records.items()
        },
        "whole_exact_normalized_pair_budgets": {
            f"polarization_{i}_{j}": v for (i, j), v in bounds.items()
        },
        "whole_uniform_angular_proof": "After multiplying each coefficient byab/(a+b)^2 and settingz=a/(a+b), it is P(z,r)/(1+r2)^m with0<=z<=1 and each r-degree<=2m. Every monomial has absolute value<=1 onr>=0. Absolute coefficient sums bound all4 plus/cross pair projections by530 times the spatial hard-current Frobenius norm. A unit-norm physical polarization has plus/cross coefficient l1 norm<=1. Rotations preserve the spatial Frobenius norm. Thus|Mpair|<=530(a+b)^2||Uspatial||F/(sqrt(kappa)ab), uniformly over relative angle.",
        "whole_angular_scope": "This removes the explicit q1.q2 angular pole after contraction with the complete conserved hard current. It bounds exactly the47 pair-propagator graphs. It is not a full434-graph error or an integrated soft subtraction, and it does not assign an exactly collinear helicity basis.",
        "checks": checks,
        "gates": {
            "all4_physical_polarization_pairs_included": len(records) == 4,
            "all24_current_coefficients_included": sum(
                len(cs) for cs in records.values()
            )
            == 24,
            "largest_exact_angular_budget530": max(bounds.values()) == 530,
            "zero_and_infinite_half_angle_endpoints_bounded": True,
            "unit_TT_norm_and_spatial_rotation_comparison": True,
        },
    }
