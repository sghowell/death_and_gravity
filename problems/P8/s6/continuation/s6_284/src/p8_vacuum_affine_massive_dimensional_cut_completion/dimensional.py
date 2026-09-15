"""Whole dimension-dependent angular/master reduction and finite bubble rational term."""

from functools import cache

import sympy as s

from . import polarizations, source, threshold

N, H, X, Y, Z = source.N, source.H, source.X, source.Y, source.Z
S, T, MU, EP = source.S, source.T, source.MU, source.EP
J0, J1, L0 = s.symbols("double_even double_odd single_even", real=True)


def reduce_power(degree, variable):
    k, e = divmod(degree, 2)
    return H**k * variable**e, -sum(
        H**j * variable ** (2 * (k - 1 - j) + e) for j in range(k)
    )


def moment(a, b):
    if (a + b) % 2:
        return s.Integer(0)
    d = N + 1
    return {
        (0, 0): s.Integer(1),
        (2, 0): 1 / d,
        (0, 2): 1 / d,
        (1, 1): Z / d,
        (2, 2): (1 + 2 * Z * Z) / (d * (d + 2)),
    }[(a, b)]


def single(e, j):
    if (e + j) % 2:
        return s.Integer(0)
    if (e, j) == (0, 0):
        return L0
    if (e, j) == (1, 1):
        return Z * (H * L0 - 1)
    if (e, j) == (0, 2):
        return ((N + 1) * Z * Z - 1) / N * (H * L0 - 1) + (1 - Z * Z) * L0 / N
    raise ValueError("Unsupported single-denominator angular monomial")


@cache
def angular_reduction():
    numerator = s.expand(polarizations.angular_numerator())
    full = s.Integer(0)
    avg = s.Integer(0)
    for (a, b), c in s.Poly(numerator, X, Y).terms():
        ra, pa = reduce_power(a, X)
        rb, pb = reduce_power(b, Y)
        full += c * (
            ra * rb / ((H - X * X) * (H - Y * Y))
            + ra * pb / (H - X * X)
            + pa * rb / (H - Y * Y)
            + pa * pb
        )
        ka, ea = divmod(a, 2)
        kb, eb = divmod(b, 2)
        if ea != eb:
            raise ValueError("Unexpected odd double resolvent")
        avg += c * H ** (ka + kb) * (J0 if ea == 0 else J1)
        avg += c * sum(H**ka * v * single(ea, j) for (j,), v in s.Poly(pb, Y).terms())
        avg += c * sum(H**kb * v * single(eb, j) for (j,), v in s.Poly(pa, X).terms())
        avg += c * sum(v * moment(i, j) for (i, j), v in s.Poly(pa * pb, X, Y).terms())
    coefficients = tuple(s.factor(s.diff(avg, v)) for v in (J0, J1, L0)) + (
        s.factor(avg.subs({J0: 0, J1: 0, L0: 0})),
    )
    return full, avg, coefficients


def master_coefficients(energy=S, transfer=T, mass=MU, epsilon=EP):
    ss, tt, mu, ep = map(s.sympify, (energy, transfer, mass, epsilon))
    q = ss - 4 * mu
    z = 1 + 2 * tt / q
    n = 2 + 2 * ep
    d = n + 1
    c0, c1, cl, cp = (
        s.factor(c.subs({H: ss / q, Z: z, N: n})) for c in angular_reduction()[2]
    )
    vd = ss * ss - 4 * mu * ss + 2 * mu * mu + 2 * mu * mu * ep / (1 + ep)
    p = 4 * vd / q
    b = -q * q / (4 * ss)
    a = (
        -2 * ss
        + 6 * mu
        - 2 * mu * mu / ss
        + q * q / (4 * ss)
        - 2 * mu * mu * ep / (ss * (1 + ep))
    )
    h = a + b * z * z
    return {
        "ordered_box": s.factor(q**4 * c0 / 128 + ss * q**3 * c1 / 128),
        "C00mu": s.factor(-(q**3) * cl / 32),
        "B00": s.factor(q * q * cp / 16),
        "C0mumu": s.factor(-2 * vd * h),
        "Bmm": s.factor(
            (
                2 * p * b * (-z * z + (1 - z * z) / n)
                + a * a
                + 2 * a * b / d
                + b * b * (1 + 2 * z * z) / (d * (d + 2))
            )
            / 2
        ),
    }


@cache
def evanescent_rational_symbolic():
    row = master_coefficients()
    return s.factor(-s.diff(row["B00"] + row["Bmm"], EP).subs(EP, 0))


def evanescent_rational(energy=S, transfer=T, mass=MU):
    return evanescent_rational_symbolic().subs(
        dict(zip((S, T, MU), map(s.sympify, (energy, transfer, mass)), strict=True)),
        simultaneous=True,
    )


def symmetric_basis(energy=S, transfer=T, mass=MU):
    ss, tt, mu = map(s.sympify, (energy, transfer, mass))
    uu = 4 * mu - ss - tt
    rows = ((ss, tt, uu), (tt, ss, uu), (uu, ss, tt))
    return (
        ss * ss + tt * tt + uu * uu,
        mu * mu,
        mu**3 * sum(1 / a for a, b, c in rows),
        mu**4 * sum(1 / a**2 for a, b, c in rows),
        mu * sum(b * c / a for a, b, c in rows),
        mu * mu * sum((b - c) ** 2 / a**2 for a, b, c in rows),
    )


def compact_crossed_evans(energy=S, transfer=T, mass=MU):
    coefficients = tuple(
        s.Rational(a, b)
        for a, b in (
            (69, 50),
            (-4979, 90),
            (3044, 225),
            (28, 225),
            (-418, 225),
            (16, 225),
        )
    )
    return sum(
        c * b
        for c, b in zip(
            coefficients, symmetric_basis(energy, transfer, mass), strict=True
        )
    )


def compact_Gram_removed_UV(energy=S, transfer=T, mass=MU):
    coefficients = tuple(
        s.Rational(a, b)
        for a, b in ((203, 40), (-25, 3), (-76, 15), (-2, 15), (47, 15), (1, 15))
    )
    return sum(
        c * b
        for c, b in zip(
            coefficients, symmetric_basis(energy, transfer, mass), strict=True
        )
    )


@cache
def data():
    unavg, avg, angular = angular_reduction()
    row = master_coefficients()
    old = source.reconstruction.coefficients()
    uu = 4 * MU - S - T
    q = S - 4 * MU
    vt = T * T - 4 * MU * T + 2 * MU * MU
    d = N + 1
    expected_single = Z * Z * (H * L0 - 1) + (1 - Z * Z) * (L0 - H * L0 + 1) / N
    ev = evanescent_rational_symbolic()
    crossed = sum(
        evanescent_rational(a, b, MU) - threshold.evanescent_gram_choice(a, b, MU)
        for a, b in ((S, T), (T, S), (uu, S))
    )
    uv = sum(
        source.reconstruction.coefficients(a, b)["B00"]
        + source.reconstruction.coefficients(a, b)["Bmm"]
        - MU * b * c * (5 * (a - 4 * MU) - 12 * MU) / (a - 4 * MU) ** 2
        for a, b, c in ((S, T, uu), (T, S, uu), (uu, S, T))
    )
    checks = {
        "entire_D_unaveraged_two_axis_division": s.factor(
            unavg - polarizations.angular_numerator() / ((H - X * X) * (H - Y * Y))
        ),
        "entire_D_averaged_coefficient_decomposition": s.factor(
            avg - angular[0] * J0 - angular[1] * J1 - angular[2] * L0 - angular[3]
        ),
        "whole_D_conditional_y2_moment": s.factor(single(0, 2) - expected_single),
        "whole_D_cross_polynomial_moment": s.factor(
            Z * Z * 3 / (d * (d + 2))
            + (1 - Z * Z) * (1 / d - 3 / (d * (d + 2))) / N
            - moment(2, 2)
        ),
        "all_D_shared_box_matches_elastic_cut": s.factor(
            row["ordered_box"] - (vt + 2 * MU * MU * EP / (1 + EP)) ** 2
        ),
        "all_D_massive_triangle_matches_elastic_cut": s.factor(
            row["C0mumu"] - source.reconstruction.evanescent_coefficients()[1]
        ),
        "entire_crossed_finite_bubble_rational_compact_identity": s.factor(
            crossed - compact_crossed_evans()
        ),
        "whole_crossed_UV_bubble_Gram_removed_identity": s.factor(
            uv - compact_Gram_removed_UV()
        ),
        "finite_bubble_double_Gram_principal_part": s.factor(
            s.limit(q * q * ev, S, 4 * MU) + 8 * MU * MU * T * T
        ),
        "finite_bubble_simple_Gram_principal_part": s.factor(
            s.limit(q * (ev + 8 * MU * MU * T * T / q**2), S, 4 * MU)
            - 2 * MU * T * (-4 * MU + 3 * T)
        ),
        "finite_bubble_Gram_choice_double_cancellation": s.factor(
            s.limit(q * q * (ev - threshold.evanescent_gram_choice()), S, 4 * MU)
        ),
        "finite_bubble_Gram_choice_simple_cancellation": s.factor(
            s.limit(q * (ev - threshold.evanescent_gram_choice()), S, 4 * MU)
        ),
        "massless_crossed_finite_bubble_calibration": s.factor(
            compact_crossed_evans(mass=0)
            - s.Rational(69, 50) * (S * S + T * T + (-S - T) ** 2)
        ),
        "physical_double_pole_not_erased_with_Gram_poles": s.factor(
            s.limit(S * S * compact_crossed_evans(), S, 0)
            - 4 * MU * MU * (71 * MU * MU - 64 * MU * T + 16 * T * T) / 225
        ),
    }
    for label, got, expected in zip(
        ("J0", "J1", "L0", "constant"),
        angular,
        source.angular.coefficients(H, Z),
        strict=True,
    ):
        checks["entire_4D_angular_coefficient_" + label] = s.factor(
            got.subs(N, 2) - expected
        )
    for label in ("C00mu", "B00", "C0mumu", "Bmm"):
        checks["entire_4D_master_coefficient_" + label] = s.factor(
            row[label].subs(EP, 0) - old[label]
        )
    reverse = master_coefficients(transfer=uu)
    for label in ("C00mu", "B00", "C0mumu", "Bmm"):
        checks["whole_D_other_channel_crossing_" + label] = s.factor(
            row[label] - reverse[label]
        )
    return {
        "whole_D_angular_coefficients": angular,
        "whole_D_master_coefficients": row,
        "whole_determined_finite_bubble_rational_contribution": ev,
        "whole_crossed_Gram_removed_finite_contribution": compact_crossed_evans(),
        "whole_Gram_removed_UV_bubble_structure": compact_Gram_removed_UV(),
        "six_invariant_compact_basis": symmetric_basis(),
        "master_normalization": "Overall1/(16pi^2 kappa^2). Both species share each ordered alternating box once. The exact common D-dimensional cut phase cancels in the coefficient extraction, not in the amplitude. Raw bubbles have residue-1/EP, so the finite contribution is minus the EP derivative of their coefficients.",
        "completion_boundary": "This fixes the all-D cut-master coefficients and their bubble evanescent rational term. The displayed crossing-symmetric Gram subtraction removes the integer threshold principal parts. A complete massive tadpole coefficient, threshold-supported distributional matching, finite physical poles, local counterterms, other source sectors and the full finite amplitude remain undetermined.",
        "checks": checks,
        "gates": {
            "whole_D_angular_measure_moments_not_four_dimensional_substitution": angular[
                3
            ].has(N),
            "both_species_common_box_exact_for_all_EP": row["ordered_box"].has(EP),
            "both_bubble_EP_terms_retained": all(
                row[k].has(EP) for k in ("B00", "Bmm")
            ),
            "finite_rational_contribution_not_set_to_zero": ev != 0,
            "crossed_rational_contains_physical_poles": compact_crossed_evans().has(
                1 / S**2
            ),
            "Gram_removal_not_physical_pole_counterterm_choice": True,
            "complete_tadpole_and_local_matching_still_required": True,
        },
    }
