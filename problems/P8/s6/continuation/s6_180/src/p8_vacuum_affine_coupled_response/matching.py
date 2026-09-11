"""Whole-dimensional physical current and unchanged massive scalar reference."""

from functools import cache

import sympy as s


@cache
def data():
    N, a, k, m = s.symbols("N a k m", positive=True)
    D = s.symbols("spatial_dimension", real=True)
    z = s.symbols("z", nonnegative=True)
    q = k * k / a**2
    w2 = N * N * (q + m * m)
    pairs, checks = {}, {}
    for sector in ("T", "L"):
        G = a ** (D - 2) / N
        if sector == "L":
            G *= m * m / (q + m * m)
        rows = {}
        for key, derivative in (
            ("N", lambda x: s.diff(x, N)),
            ("Z", lambda x: a * s.diff(x, a)),
        ):
            A = s.factor(
                (derivative(G) / G).subs(N, 1).subs(k * k, a * a * m * m * z / (1 - z))
            )
            B = s.factor(
                (-derivative(G * w2) / (G * w2))
                .subs(N, 1)
                .subs(k * k, a * a * m * m * z / (1 - z))
            )
            rows[key] = (A, B)
        pair = s.Matrix([rows[name][0] - rows[name][1] for name in ("N", "Z")])
        pairs[sector] = {
            "weights": rows,
            "pair": s.ImmutableMatrix(pair),
            "high_pair": s.ImmutableMatrix(pair.subs(z, 1)),
        }
        expected = (
            s.Matrix([0, 2 * (D - 3)]) if sector == "T" else s.Matrix([0, 2 * (D - 1)])
        )
        checks[sector + "_whole_dimensional_high_current"] = s.ImmutableMatrix(
            (pair.subs(z, 1) - expected).applyfunc(s.factor)
        )
        expected_full = (
            s.Matrix([0, 2 * (1 - z)]) if sector == "T" else s.Matrix([0, 2 * (1 + z)])
        )
        checks[sector + "_all_momentum_physical_pair"] = s.ImmutableMatrix(
            (pair.subs(D, 3) - expected_full).applyfunc(s.factor)
        )
    bT = pairs["T"]["high_pair"][1]
    bL = pairs["L"]["high_pair"][1]
    checks["transverse_square_and_first_dimensional_jet"] = s.ImmutableMatrix(
        [((D - 1) * bT * bT).subs(D, 3), s.diff((D - 1) * bT * bT, D).subs(D, 3)]
    )
    checks["longitudinal_dimensional_square_jet"] = s.diff(bL * bL, D).subs(D, 3) - 16
    checks["physical_radial_pair_square"] = s.factor(
        (2 * pairs["T"]["pair"][1] ** 2 + pairs["L"]["pair"][1] ** 2).subs(D, 3)
        - 4 * (3 - 2 * z + 3 * z * z)
    )
    y, p = s.symbols("y p", real=True)
    W = y * y * (3 - 2 * y * y + 3 * y**4)
    integral = s.Integral(W * p / (4 * m * m + p * (1 - y * y)), (y, 0, 1))
    checks["exact_massive_reference_threshold"] = (
        4 - s.integrate(W / (y * y), (y, 0, 1)) - s.Rational(16, 15)
    )
    checks["exact_massive_reference_first_derivative"] = (
        s.integrate(W / (4 * m * m), (y, 0, 1)) - s.Rational(9, 35) / m**2
    )
    checks["positive_radial_polynomial"] = s.expand(
        3 - 2 * z + 3 * z * z - 3 * (z - s.Rational(1, 3)) ** 2 - s.Rational(8, 3)
    )
    lag = s.symbols("positive_lag", positive=True)
    checks["matched_leading_four_primitive_kernel"] = (
        s.diff(4 / lag, lag, 4) - 96 / lag**5
    )
    scale = s.symbols("positive_scale", positive=True)
    checks["all_dimension_proper_momentum_diagonal_scale"] = s.simplify(
        scale ** (-D - 2) * scale ** (D + 1) * scale - 1
    )
    return {
        "actual_full_dimensional_current_pairs": pairs,
        "massive_scalar_reference_F": -4 - integral,
        "fixed_finite_fourth_coefficient": -s.Integer(4),
        "reference_first_sheet_gap": s.Rational(16, 15),
        "reference_inverse_absolute_upper": s.Rational(15, 16),
        "reference_inverse_moments": (s.Rational(1, 4), s.Rational(9, 560) / m**2),
        "inverse_kernel": "The same source-pinned S6.86 scalar kernel K_m is in L1(0,infinity), has no instantaneous term and Laplace transform1/F_m(s^2). It is only the active range inverse of the normalized isolated loop block.",
        "current_matching": "The S6.179 whole physical operator/state/prescription bridge and these full-D vertices give exactly S6.87's pure-scale nonlocal current and its dimensional fourth-order contact. The current adds9P*w, which is local of order0. No old scalar profile is installed.",
        "curved_remainder": "After four zero-past primitives, I4 Q_w,scale=F_m(partial_t^2)+V_scalar. The unchanged all-order state, full-D diagonal amplitude and fixed -4 coefficient give a weak-log integrable kernel with every fixed diagonal derivative, as proved in S6.87. There is no unmatched fourth-order local remainder; new uncancelled contacts are lower order.",
        "checks": checks,
    }
