"""Read-only comparison of additional RATE4 parent-input prescriptions.

RATE4-D8 and RATE5-D10 are PROPOSALS, not adopted or matched parents.
The finite-dimensional algebra and conditional tail/error bounds below
do not supply their physical inputs, a full V/G functional, or bounce data.
"""

import hashlib
import json
from fractions import Fraction as F
from functools import cache
from itertools import product
from math import prod

import sympy as s
from p8_match1_audit import POINTS, fraction_row, solve_fraction
from p8_match1_rate_input import REPO, exact
from p8_rate4_candidate import dot, original_borns, rational
from p8_rate4_remainder import contact, local_words
from p8_rate4_remainder import inputs as remainder_inputs

REMAINDER_SHA256 = "14920b563db7f640a65acb31c859ead22de6d44ea1158bfc26a0157bb3eae341"
FIVE_POINTS = POINTS + ((10, -2, -4),)
EIGHT_POINTS = POINTS + ((10, -1, -5), (12, -1, -7), (12, -2, -6), (16, -1, -11))
FIVE_TARGET = (0, 2, 0, 32, 0)


def inputs():
    manifest = remainder_inputs()
    path = REPO / "scripts/p8_rate4_remainder.py"
    assert hashlib.sha256(path.read_bytes()).hexdigest() == REMAINDER_SHA256
    manifest[str(path.relative_to(REPO))] = REMAINDER_SHA256
    return manifest


def samples(values, degree_limit):
    local_case(degree_limit)
    result = tuple(map(exact, values))
    if len(result) != degree_limit:
        raise ValueError("Require one exact input for each declared matching point")
    return result


def nonnegative(value):
    result = exact(value)
    if result < 0:
        raise ValueError("An absolute error bound cannot be negative")
    return result


@cache
def five_weights():
    rows = [local_words(point) for point in FIVE_POINTS]
    return tuple(solve_fraction(list(zip(*rows, strict=True)), FIVE_TARGET))


@cache
def local_case(degree_limit):
    if type(degree_limit) is not int or degree_limit not in (4, 5):
        raise ValueError("Declare channel degree 4 (D8) or 5 (D10)")
    if degree_limit == 5:
        return FIVE_POINTS, five_weights()
    rows = [local_words(point)[:4] for point in POINTS]
    return POINTS, tuple(solve_fraction(list(zip(*rows, strict=True)), FIVE_TARGET[:4]))


def polynomial_fit(*, degree_limit, residuals):
    """Fit COMPLETE supplied residuals in the explicitly declared truncation.

    This is not a computation of the physical residuals or their provenance.
    Unknown hard terms, tail values and rate conversions have no zero default.
    """
    points, _ = local_case(degree_limit)
    return tuple(
        solve_fraction(
            [local_words(point)[:degree_limit] for point in points],
            samples(residuals, degree_limit),
        )
    )


@cache
def tail_constant(degree_limit):
    """Projection bound per unit of N64, beyond the declared channel degree.

    This norm is an ADDITIONAL parent input, not inferred from power counting.
    Known nonlocal loops, heavy exchange and the Newton shape are excluded.
    """
    forward = F(384, 64**6)
    points, weights = local_case(degree_limit)
    samples = []
    for point in points:
        x, y = F(sum(a * a for a in point), 64**2), F(prod(point), 64**3)
        assert 0 <= x < 1 and 0 <= y < 1
        samples.append(max(x**3, (x if degree_limit == 4 else x**2) * y, y**2))
    return forward + dot(tuple(map(abs, weights)), samples)


def nuisance_bound(
    *, degree_limit, mass_ratio, residue_abs, mass_insertion_abs, newton_abs
):
    """Conditional first-order Q/J/T projection error, NOT an all-loop bound.

    The amplitude coordinates are r*Q+eta*J+D_N*T. In the linear mass
    displacement convention eta=-g^2*delta_n. A physical interpretation
    requires a common reference and properly matched pole/source data.
    """
    n = exact(mass_ratio)
    if n < 32:
        raise ValueError("Require n>=32")
    residue, mass, newton = map(
        nonnegative, (residue_abs, mass_insertion_abs, newton_abs)
    )
    ratio = 16 / n
    power = degree_limit + 1
    points, weights = local_case(degree_limit)
    q_tail = 3 * F(16**power) / (n**power * (n - 16))
    j_tail = 3 / n**2 * ratio**power * (power + 1 - power * ratio) / (1 - ratio) ** 2
    factor = 1 + sum(map(abs, weights))
    # L(T)=0 only in the independently required pole-subtraction convention.
    newton_projection = -dot(weights, [fraction_row(q, n)[3] for q in points])
    return factor * (residue * q_tail + mass * j_tail) + abs(newton_projection) * newton


def conditional_budget(
    *,
    degree_limit,
    mass_ratio,
    sample_errors,
    known_L_error,
    higher_L_error,
    tail_norm64,
    residue_abs,
    mass_insertion_abs,
    newton_abs,
):
    """Combine complete SUPPLIED errors in the declared D8 or D10 proposal.

    sample_errors includes rate/hard/IR/radiation/omitted-order SAMPLE
    uncertainties. higher_L_error separately covers omitted orders in L.
    No physical input, existence or admissibility of the full L is certified.
    """
    _, weights = local_case(degree_limit)
    errors = samples(sample_errors, degree_limit)
    if any(error < 0 for error in errors):
        raise ValueError("Require nonnegative sample errors")
    known, higher, tail = map(nonnegative, (known_L_error, higher_L_error, tail_norm64))
    return (
        dot(tuple(map(abs, weights)), errors)
        + known
        + higher
        + tail_constant(degree_limit) * tail
        + nuisance_bound(
            degree_limit=degree_limit,
            mass_ratio=mass_ratio,
            residue_abs=residue_abs,
            mass_insertion_abs=mass_insertion_abs,
            newton_abs=newton_abs,
        )
    )


def positive_after(expression, variable, lower):
    z = s.Symbol("nonnegative_z")
    for polynomial in s.fraction(s.factor(expression)):
        assert all(
            c > 0 for c in s.Poly(polynomial.subs(variable, z + lower), z).all_coeffs()
        )


@cache
def eight_symbolic():
    """Exact rank-6 constant system plus a 2x2 heavy Schur complement."""
    n = s.Symbol("n", positive=True)
    rows = [
        list(local_words(q)) + [rational(fraction_row(q, F(32))[3])]
        for q in EIGHT_POINTS
    ]
    constant = s.Matrix(rows).T
    assert constant[:, :6].det() == -s.Rational(79793142996096, 5)
    base = s.Matrix(
        list(constant[:, :6].inv() * s.Matrix([0, 2, 0, 32, 0, 0])) + [0, 0]
    )
    null = constant.nullspace()
    assert len(null) == 2
    heavy = s.Matrix(
        [sum((s.Integer(a) + 2) / (n - a) for a in q) for q in EIGHT_POINTS]
    )
    mass = s.Matrix([sum(1 / (n - a) ** 2 for a in q) for q in EIGHT_POINTS])
    h = [s.factor(vector.dot(heavy)) for vector in null]
    j = [s.factor(vector.dot(mass)) for vector in null]
    determinant = s.factor(h[0] * j[1] - h[1] * j[0])
    rh = s.factor(2 * (n + 2) / (n - 2) ** 3 - base.dot(heavy))
    rj = s.factor(6 / (n - 2) ** 4 - base.dot(mass))
    c0 = s.factor((rh * j[1] - h[1] * rj) / determinant)
    c1 = s.factor((h[0] * rj - rh * j[0]) / determinant)
    weights = s.Matrix(
        [s.factor(base[i] + c0 * null[0][i] + c1 * null[1][i]) for i in range(8)]
    )
    assert (constant * weights - s.Matrix([0, 2, 0, 32, 0, 0])).applyfunc(
        s.factor
    ) == s.zeros(6, 1)
    assert s.factor(weights.dot(heavy) - 2 * (n + 2) / (n - 2) ** 3) == 0
    assert s.factor(weights.dot(mass) - 6 / (n - 2) ** 4) == 0
    return n, determinant, weights


def eight_rows(n):
    return [
        list(local_words(q))
        + [
            fraction_row(q, n)[3],
            fraction_row(q, n)[1],
            sum(1 / (n - a) ** 2 for a in q),
        ]
        for q in EIGHT_POINTS
    ]


def audit():
    before = inputs()
    n0, g, kappa, lam, _ = original_borns()
    weights = five_weights()
    matrix = s.Matrix([local_words(q) for q in FIVE_POINTS])
    assert matrix.det() != 0
    symbolic = (s.Matrix([FIVE_TARGET]) * matrix.inv()).tolist()[0]
    comparisons = 0
    for direct, independent in zip(symbolic, weights, strict=True):
        assert direct == rational(independent)
        comparisons += 1
    assert sum(map(abs, weights)) == F(630347, 773752) < 1
    assert -dot(weights, [contact(q) for q in FIVE_POINTS]) == 902256
    coordinate_norms = []
    for index in range(5):
        row = solve_fraction(
            list(zip(*[local_words(q) for q in FIVE_POINTS], strict=True)),
            [int(j == index) for j in range(5)],
        )
        coordinate_norms.append(sum(map(abs, row)))
    assert max(coordinate_norms) == F(1906116, 96719) < 20
    fixture_checks = 0
    for coefficients in ((F(1, 3), -2, F(5, 7), 11, F(-3, 2)), (3, 7, -13, 2, 17)):
        residuals = [dot(local_words(q), coefficients) for q in FIVE_POINTS]
        assert polynomial_fit(degree_limit=5, residuals=residuals) == tuple(
            coefficients
        )
        fixture_checks += 5

    # Both exact rational heavy expansions through channel degree five are in span.
    v, x = s.symbols("v x")
    a, heavy_mass = s.symbols("a heavy_mass")
    q_remainder = a**6 / (heavy_mass**6 * (heavy_mass - a))
    assert (
        s.cancel(
            1 / (heavy_mass - a)
            - sum(a**degree / heavy_mass ** (degree + 1) for degree in range(6))
            - q_remainder
        )
        == 0
    )
    assert (
        s.cancel(
            1 / (heavy_mass - a) ** 2
            - sum(
                (degree + 1) * a**degree / heavy_mass ** (degree + 2)
                for degree in range(6)
            )
            + s.diff(q_remainder, heavy_mass)
        )
        == 0
    )
    for degree in range(6):
        pk = (2 + v) ** degree + (2 - v) ** degree + s.Integer(0) ** degree
        expected = s.expand(pk).coeff(v, 2)
        assert (
            rational(dot(weights, [sum(a**degree for a in q) for q in FIVE_POINTS]))
            == expected
        )
    assert (
        s.cancel(
            1 / (1 - x) ** 2
            - sum((k + 1) * x**k for k in range(6))
            - x**6 * (7 - 6 * x) / (1 - x) ** 2
        )
        == 0
    )
    assert tail_constant(5) == F(217375713, 6490702217216) < F(1, 20000)
    tail_checks = 0
    for i, j in product(range(13), repeat=2):
        degree = 2 * i + 3 * j
        if degree <= 5:
            continue
        derivative = 2 * i * 8 ** (i - 1) if j == 0 else 0
        sample = [sum(a * a for a in q) ** i * prod(q) ** j for q in FIVE_POINTS]
        projected = (derivative - dot(weights, sample)) / 64**degree
        assert abs(projected) <= tail_constant(5)
        tail_checks += 1
    # i=3 maximizes the normalized forward derivative among all i>=3.
    assert F(4, 3) * F(8, 64**2) == F(1, 384) < 1

    # The fifth point remains in the same massive physical nonforward window.
    q = FIVE_POINTS[-1]
    d = n0 - 2
    ctree = -(g**2) * (3 / d - 2 / d**2)
    born5 = ctree + g**2 * sum(1 / (n0 - a) for a in q) - fraction_row(q, n0)[3] / kappa
    assert 232 * lam < born5 < 233 * lam
    assert all(
        sum(q) == 4 and F(25, 4) <= q[0] <= 16 and min(-q[1], -q[2]) >= 1
        for q in FIVE_POINTS + EIGHT_POINTS
    )
    nt = -dot(weights, [fraction_row(q, n0)[3] for q in FIVE_POINTS])
    assert nt == F(26901281, 79585920)
    nuisance = nuisance_bound(
        degree_limit=5,
        mass_ratio=n0,
        residue_abs=g**2 / 4,
        mass_insertion_abs=g**2 * n0 / 4,
        newton_abs=F(1, 4 * kappa),
    )
    assert nuisance < lam / 10**200
    # Sufficient allocation ONLY: the required parent data are not supplied.
    sample_error = F(650, 2 * 10000) * lam + lam / 100
    allocation = {
        "degree_limit": 5,
        "mass_ratio": n0,
        "sample_errors": (sample_error,) * 5,
        "known_L_error": lam / 10,
        "higher_L_error": lam / 4,
        "tail_norm64": 2000 * lam,
        "residue_abs": g**2 / 4,
        "mass_insertion_abs": g**2 * n0 / 4,
        "newton_abs": F(1, 4 * kappa),
    }
    assert conditional_budget(**allocation) < F(987, 2000) * lam < lam / 2

    # A minimum-change alternative retains all four targets but bounds SU in the tail.
    points4, weights4 = local_case(4)
    matrix4 = s.Matrix([local_words(q)[:4] for q in points4])
    assert matrix4.det() == 10011648
    symbolic4 = (s.Matrix([FIVE_TARGET[:4]]) * matrix4.inv()).tolist()[0]
    for direct, independent in zip(symbolic4, weights4, strict=True):
        assert direct == rational(independent)
        comparisons += 1
    assert sum(map(abs, weights4)) == F(3544, 3259)
    coordinate_norms4 = [
        sum(abs(value) for value in matrix4.inv().row(i)) for i in range(4)
    ]
    assert max(coordinate_norms4) == s.Rational(85231, 3259) < 27
    assert tail_constant(4) == F(332221287, 6998649208832) < F(1, 20000)
    for coefficients in ((F(1, 3), -2, F(5, 7), 11), (3, 7, -13, 2)):
        residuals = [dot(local_words(q)[:4], coefficients) for q in points4]
        assert polynomial_fit(degree_limit=4, residuals=residuals) == tuple(
            coefficients
        )
        fixture_checks += 4
    for degree in range(5):
        pk = (2 + v) ** degree + (2 - v) ** degree + s.Integer(0) ** degree
        assert rational(
            dot(weights4, [sum(a**degree for a in q) for q in points4])
        ) == s.expand(pk).coeff(v, 2)
    assert (
        s.cancel(
            1 / (heavy_mass - a)
            - sum(a**degree / heavy_mass ** (degree + 1) for degree in range(5))
            - a**5 / (heavy_mass**5 * (heavy_mass - a))
        )
        == 0
    )
    assert (
        s.cancel(
            1 / (1 - x) ** 2
            - sum((k + 1) * x**k for k in range(5))
            - x**5 * (6 - 5 * x) / (1 - x) ** 2
        )
        == 0
    )
    for i, j in product(range(13), repeat=2):
        degree = 2 * i + 3 * j
        if degree <= 4:
            continue
        derivative = 2 * i * 8 ** (i - 1) if j == 0 else 0
        sampled = [sum(a * a for a in q) ** i * prod(q) ** j for q in points4]
        assert abs((derivative - dot(weights4, sampled)) / 64**degree) <= tail_constant(
            4
        )
        tail_checks += 1
    assert -dot(weights4, [fraction_row(q, n0)[3] for q in points4]) == F(
        536987, 2346480
    )
    allocation4 = dict(allocation, degree_limit=4, sample_errors=(sample_error,) * 4)
    nuisance4 = nuisance_bound(
        degree_limit=4,
        mass_ratio=n0,
        residue_abs=g**2 / 4,
        mass_insertion_abs=g**2 * n0 / 4,
        newton_abs=F(1, 4 * kappa),
    )
    assert nuisance4 < lam / 10**200
    assert conditional_budget(**allocation4) < F(1620429, 3259000) * lam < lam / 2

    # Eight low-energy rates: stable coefficient projection, unstable coordinates.
    n, determinant, eight_weights = eight_symbolic()
    positive_after(determinant, n, 32)
    signs = (-1, -1, -1, -1, 1, -1, 1, 1)
    for sign, value in zip(signs, eight_weights, strict=True):
        positive_after(sign * value, n, 1024)
    positive_after(
        1 - sum(sign * value for sign, value in zip(signs, eight_weights, strict=True)),
        n,
        1024,
    )
    for mass in (F(1024), F(2048), n0):
        target = (
            0,
            2,
            0,
            32,
            0,
            0,
            2 * (mass + 2) / (mass - 2) ** 3,
            6 / (mass - 2) ** 4,
        )
        direct = solve_fraction(list(zip(*eight_rows(mass), strict=True)), target)
        for numeric, symbolic in zip(direct, eight_weights, strict=True):
            assert rational(numeric) == symbolic.subs(n, rational(mass))
            comparisons += 1
        assert sum(map(abs, direct)) < 1
    conditioning = {}
    for index, label, exponent in ((6, "H", 774), (7, "J", 1168)):
        row = solve_fraction(
            list(zip(*eight_rows(n0), strict=True)), [int(j == index) for j in range(8)]
        )
        amplification = sum(map(abs, row))
        radius = amplification * lam / 100
        assert F(10) ** exponent < radius < F(10) ** (exponent + 1)
        # Opposite box corners attain the full coordinate radius exactly.
        signs = tuple(1 if value >= 0 else -1 for value in row)
        assert dot(row, signs) == amplification
        conditioning[label] = (
            f"1e{exponent} < radius < 1e{exponent + 1} for independent sample errors lambda/100"
        )
        if label == "J":
            precision = (g**2 * n0 / 4) / amplification / lam
            assert F(1, 10**982) < precision < F(1, 10**981)

    corners = 0
    for signs in product((-1, 1), repeat=5):
        error = dot(weights, signs)
        assert abs(error) <= sum(map(abs, weights))
        corners += 1
    for signs in product((-1, 1), repeat=4):
        assert abs(dot(weights4, signs)) <= sum(map(abs, weights4))
        corners += 1
    malformed = [
        {key: value for key, value in allocation.items() if key != "tail_norm64"},
        {key: value for key, value in allocation.items() if key != "higher_L_error"},
        dict(allocation, sample_errors=(0,) * 4),
        dict(allocation, sample_errors=(0, 0, 0, 0, -1)),
        dict(allocation, tail_norm64=None),
        dict(allocation, tail_norm64=-1),
        dict(allocation, known_L_error=0.0),
        dict(allocation, higher_L_error=True),
        dict(allocation, mass_ratio=16),
        dict(allocation, residue_abs=-1),
        dict(allocation, mass_insertion_abs=-1),
        dict(allocation, newton_abs=-1),
        {key: value for key, value in allocation.items() if key != "degree_limit"},
        dict(allocation, degree_limit=3),
    ]
    rejected = 0
    for packet in malformed:
        try:
            conditional_budget(**packet)
        except (TypeError, ValueError):
            rejected += 1
        else:
            raise AssertionError("Incomplete or malformed parent budget accepted")
    assert inputs() == before
    return {
        "milestone": "RATE4.PARENT",
        "outcome": "INPUT_DESIGN_COMPARED_NO_COMPLETE_PARENT_SUPPLIED",
        "recommendation": "RATE4-D8 parent-bounded first; RATE5-D10 if SU should be fitted rather than included in the tail bound",
        "proposal_adopted": False,
        "physical_matching_established": False,
        "protected_input_files": len(before),
        "five_points": FIVE_POINTS,
        "four_projector_norm": str(sum(map(abs, weights4))),
        "four_maximum_coordinate_error_amplification": str(max(coordinate_norms4)),
        "four_tail_projection_per_unit_N64": str(tail_constant(4)),
        "four_sufficient_conditional_error_over_lambda": "<1620429/3259000<1/2",
        "five_projector_norm": str(sum(map(abs, weights))),
        "five_maximum_coordinate_error_amplification": str(max(coordinate_norms)),
        "five_tail_projection_per_unit_N64": str(tail_constant(5)),
        "five_sufficient_conditional_error_over_lambda": "<987/2000<1/2",
        "eight_projection_norm": "<1 for n>=1024; not a bound on individual heavy coordinates",
        "eight_coordinate_conditioning": conditioning,
        "eight_J_cap_required_relative_sample_precision": "between 1e-982 and 1e-981",
        "symbolic_fraction_comparisons": comparisons,
        "synthetic_polynomial_fit_coordinates": fixture_checks,
        "tail_monomial_checks": tail_checks,
        "sample_error_corners": corners,
        "rejected_incomplete_or_invalid_inputs": rejected,
        "still_required": "actual parent bounds/data, complete hard evaluation, higher-order errors, admissible V/G functional, independent curved/source/state control; no P8 closure",
    }


if __name__ == "__main__":
    print(json.dumps(audit(), indent=2))
