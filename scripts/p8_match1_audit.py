"""Read-only MATCH-1 diagnostic; not a physical matching/positivity certificate.

Pins the cited scientific inputs, checks a normalization-preserving finite
counterterm direction and a four-value matching projector. SymPy and an
independent Fraction implementation cross-check the finite algebra. No parent
package is imported, no physical coefficient is chosen, and no file is written.
"""

import hashlib
import json
from fractions import Fraction as F
from itertools import permutations, product
from math import factorial
from pathlib import Path

import sympy as s

REPO = Path(__file__).resolve().parents[1]
PINS = {
    233: "31716f2b701c0120c328dd8a2506f677eb30edeca956d531c57975e38d8d25cb",
    238: "934b6c035fb5de0eb89ee1e1a03c2624c78a33985c0503bce883c41240b05905",
    239: "a09852126d79ca223af06555e10283328cb144a14945c127abee73b7ce8b9c05",
    240: "d4a6081d0482b2168255299cdd0080c90a1716eeffbc70f2645e6e2452e78eda",
    285: "608a9394ac4d06f0f2e5c22a52e830f023364a906ea8e811bd2b1c7415e0c019",
    289: "0a09a6d2e69c5618c6f05f943fc5a8a5da9072e75e897b0a8f430b02e2c98aa1",
    290: "191329017c78af9fcb017a3bf77df87bebb2fb47c231948f37dc1e743200c9b7",
    294: "365fbe4a68a52fd00725a6e43341a1c863638289bea1d7b12955067bdc4a22a0",
    302: "499d71347342e1d332367faf397757b617386dbf9faf4ebf17eb2be07b973780",
    336: "f8f9cd57c0240fe6008efd743d399c27dfc393a14d50501e8362289aae204db0",
    347: "df492821c0697cca4fedc5cd6bd3e46e165ba76305083c488f7b0b62b5a7b4e8",
}
POINTS = ((8, -2, -2), (10, -3, -3), (12, -4, -4), (16, -6, -6))


def frozen_inputs():
    manifest = {}
    for number, expected in PINS.items():
        root = REPO / "problems/P8/s6/continuation" / f"s6_{number}"
        reports = sorted((root / "certificates").glob("*.json"))
        assert len(reports) == 1, number
        report = reports[0]
        assert hashlib.sha256(report.read_bytes()).hexdigest() == expected, report
        manifest[str(report.relative_to(REPO))] = expected
        record = json.loads(report.read_text())
        for name, digest in record["source_sha256"].items():
            path = (root / name).resolve()
            assert path.is_relative_to(root.resolve()), path
            assert hashlib.sha256(path.read_bytes()).hexdigest() == digest, path
            manifest[str(path.relative_to(REPO))] = digest
    return manifest


def add(*polys):
    """Independent Q[S,T,mu] coefficient arithmetic, without SymPy."""
    out = {}
    for poly in polys:
        for power, value in poly.items():
            out[power] = out.get(power, F(0)) + value
    return {p: v for p, v in out.items() if v}


def scale(value, poly):
    return {p: F(value) * v for p, v in poly.items() if value * v}


def mul(left, right):
    out = {}
    for p, a in left.items():
        for q, b in right.items():
            key = tuple(x + y for x, y in zip(p, q, strict=True))
            out[key] = out.get(key, F(0)) + a * b
    return {p: v for p, v in out.items() if v}


def fraction_vertex():
    ss, tt, mu = ({(1, 0, 0): F(1)}, {(0, 1, 0): F(1)}, {(0, 0, 1): F(1)})
    uu = add(scale(4, mu), scale(-1, ss), scale(-1, tt))
    grams = {}
    for channel, pairs in zip(
        (ss, tt, uu),
        (((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))),
        strict=True,
    ):
        for i, j in pairs:
            grams[i, j] = grams[j, i] = scale(F(1, 2), add(channel, scale(-2, mu)))
    literal = add(
        *(mul(grams[i, j], grams[k, ell]) for i, j, k, ell in permutations(range(4)))
    )
    mu2 = mul(mu, mu)
    deformation = add(literal, scale(-F(24, 9), mu2))
    expected = scale(
        2, add(mul(ss, ss), mul(tt, tt), mul(uu, uu), scale(-F(16, 3), mu2))
    )
    assert deformation == expected
    return deformation


def solve_fraction(matrix, rhs):
    """Exact independent Gauss-Jordan solve; inputs are diagnostic rational data."""
    size = len(rhs)
    rows = [[F(x) for x in row] + [F(y)] for row, y in zip(matrix, rhs, strict=True)]
    assert len(rows) == size and all(len(row) == size + 1 for row in rows)
    for col in range(size):
        pivot = next(i for i in range(col, size) if rows[i][col])
        rows[col], rows[pivot] = rows[pivot], rows[col]
        divisor = rows[col][col]
        rows[col] = [v / divisor for v in rows[col]]
        for i in range(size):
            if i != col:
                factor = rows[i][col]
                rows[i] = [
                    a - factor * b for a, b in zip(rows[i], rows[col], strict=True)
                ]
    return [row[-1] for row in rows]


def fraction_row(point, x):
    sigma = sum(a * a for a in point)
    heavy = sum((F(a) + 2) / (x - a) for a in point)
    newton = sum(
        F(2 - 2 * a - point[(i + 1) % 3] * point[(i + 2) % 3], a)
        for i, a in enumerate(point)
    )
    return list(map(F, (sigma, heavy, 1, newton)))


def positive_on_halfline(expression, x):
    """Exact sufficient polynomial proof for every x>=32, not point sampling."""
    z = s.Symbol("nonnegative_z")
    numerator, denominator = s.fraction(s.factor(expression))
    for value in (numerator, denominator):
        poly = s.Poly(value.subs(x, z + 32), z)
        assert all(c > 0 for c in poly.all_coeffs()), poly


def audit():
    before = frozen_inputs()
    ss, tt, mu, v, kappa, eta = s.symbols("ss tt mu v kappa eta")
    uu = 4 * mu - ss - tt
    grams = {}
    for channel, pairs in zip(
        (ss, tt, uu),
        (((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))),
        strict=True,
    ):
        for i, j in pairs:
            grams[i, j] = grams[j, i] = (channel - 2 * mu) / 2
    literal = sum(
        grams[i, j] * grams[k, ell] for i, j, k, ell in permutations(range(4))
    )
    polynomial = s.expand(literal - s.Rational(24, 9) * mu**2)
    assert (
        s.expand(polynomial - 2 * (ss**2 + tt**2 + uu**2 - s.Rational(16, 3) * mu**2))
        == 0
    )
    sympy_coefficients = {
        powers: F(int(value.p), int(value.q))
        for powers, value in s.Poly(polynomial, ss, tt, mu).terms()
    }
    assert sympy_coefficients == fraction_vertex()
    amplitude = eta * polynomial / kappa**2
    assert s.factor(amplitude.subs({ss: 4 * mu / 3, tt: 4 * mu / 3})) == 0
    crossed = amplitude.subs(ss, 2 * mu - tt / 2 + v)
    assert s.factor(s.diff(crossed, v, 2) / 2 - 4 * eta / kappa**2) == 0
    mass, coupling, curvature = s.symbols("mass coupling curvature")
    heavy_amplitude = (
        2
        * coupling
        * curvature
        / kappa
        * sum((a + 2 * mu) / (mass - a) for a in (ss, tt, uu))
    )
    heavy_forward = heavy_amplitude.subs(ss, 2 * mu - tt / 2 + v)
    heavy_coefficient = s.diff(heavy_forward, v, 2).subs({tt: 0, v: 0}) / 2
    assert (
        s.factor(
            heavy_coefficient
            - 4
            * coupling
            * curvature
            * (mass + 2 * mu)
            / (kappa * (mass - 2 * mu) ** 3)
        )
        == 0
    )

    # Exact factor identities used in the written 1024-fold vacuum/clock proof.
    p, q = s.symbols("p q")
    assert s.cancel(q / (p + q) - 1 + p / (p + q)) == 0
    assert s.cancel(q / (p + q) * (p + q) - q) == 0

    k0 = 10**800
    lam = F(1, 10**600)
    d = F(10**200, 512)
    n = d + 2
    cubic = F(1, 8192)
    eta_bound = 2 * 10**1000
    assert 2 * cubic**2 / d**3 == 4 * lam
    assert F(4 * eta_bound, k0**2) == 8 * lam
    bare_quartic_margin = cubic**2 * (d - 1) / (6 * d**2 * (d + 2))
    assert F(eta_bound, 9 * k0**2) < bare_quartic_margin / 4
    supremum = 6 * F(eta_bound, k0) * F(9, 55) ** 1024
    for i in range(5):
        for j in range(5 - i):
            assert supremum * factorial(i) * factorial(j) * 64 ** (i + j) < F(
                1, 10**590
            )
    assert 4 * lam * (1 + F(1, 10**203)) - 8 * lam < 0
    assert 4 * lam + 8 * lam > 0

    x = s.Symbol("heavy_to_light_mass_squared", positive=True)
    rows = []
    for point in POINTS:
        assert sum(point) == 4 and point[0] >= 8 and point[1] == point[2] < 0
        sigma = sum(a * a for a in point)
        heavy = sum((s.Integer(a) + 2) / (x - a) for a in point)
        newton = sum(
            (2 - 2 * a - point[(i + 1) % 3] * point[(i + 2) % 3]) / s.Integer(a)
            for i, a in enumerate(point)
        )
        rows.append([sigma, heavy, 1, newton])
    matrix = s.Matrix(rows)
    target = s.Matrix([[2, 2 * (x + 2) / (x - 2) ** 3, 0, 0]])
    weights = [s.factor(value) for value in target * matrix.inv()]
    assert all(s.factor(value) == 0 for value in s.Matrix([weights]) * matrix - target)
    determinant_polynomial = (
        3659 * x**4 + 309616 * x**3 + 815324 * x**2 + 30852016 * x - 93987840
    )
    determinant = (
        -2
        * determinant_polynomial
        / (5 * (x - 16) * (x - 12) * (x - 10) * (x - 8) * (x + 3) * (x + 4) * (x + 6))
    )
    assert s.factor(matrix.det() - determinant) == 0
    positive_on_halfline(-determinant, x)
    signs = (1, -1, 1, -1)
    for sign, weight in zip(signs, weights, strict=True):
        positive_on_halfline(sign * weight, x)
    norm = sum(sign * weight for sign, weight in zip(signs, weights, strict=True))
    positive_on_halfline(15 - norm, x)
    limits = (
        s.Rational(42447, 14636),
        -s.Rational(24900, 3659),
        s.Rational(65923, 14636),
        -s.Rational(4385, 7318),
    )
    assert tuple(s.limit(weight, x, s.oo) for weight in weights) == limits

    comparisons = 0
    for mass_ratio in (F(32), F(64), F(1000), n):
        numeric = [fraction_row(point, mass_ratio) for point in POINTS]
        target_f = [F(2), 2 * (mass_ratio + 2) / (mass_ratio - 2) ** 3, F(0), F(0)]
        fraction_weights = solve_fraction(list(zip(*numeric, strict=True)), target_f)
        for expected, symbolic in zip(fraction_weights, weights, strict=True):
            actual = symbolic.subs(
                x, s.Rational(mass_ratio.numerator, mass_ratio.denominator)
            )
            assert expected == F(int(actual.p), int(actual.q))
            comparisons += 1
        assert sum(abs(w) for w in fraction_weights) < 15
        # All four columns are exercised: nuisance constant and Newton are not zeroed.
        for coefficients in (
            (F(1, 3), F(-7, 5), F(11), F(13, 7)),
            (F(-2), F(9), F(-4), F(6)),
        ):
            values = [
                sum(a * b for a, b in zip(row, coefficients, strict=True))
                for row in numeric
            ]
            recovered = sum(
                a * b for a, b in zip(fraction_weights, values, strict=True)
            )
            expected = sum(a * b for a, b in zip(target_f, coefficients, strict=True))
            assert recovered == expected
        for errors in product((-1, 1), repeat=4):
            assert (
                abs(sum(w * e for w, e in zip(fraction_weights, errors, strict=True)))
                < 15
            )

    after = frozen_inputs()
    assert before == after
    return {
        "route_verdict": "CONDITIONAL_ONLY_FROM_RETAINED_DATA_NOT_A_PHYSICAL_EXCLUSION",
        "pinned_reports": len(PINS),
        "pinned_report_and_source_files": len(before),
        "literal_vertex": "SymPy and independent Fraction polynomial ring agree",
        "comparison": "OS4 value unchanged; clock jets through order1023 zero; four-jet tube bound <1e-590",
        "selected_reference_sign_control": "eta=+-2e1000 gives delta_b20=+-8lambda; not full physical b20",
        "matching_points_in_units_of_mu": POINTS,
        "projector_columns": ["sigma2", "RH_shape", "constant", "Newton_shape"],
        "projector_error": "<15*epsilon/mu^2 for n>=32mu, conditional on four residual error bounds",
        "independent_exact_weight_comparisons": comparisons,
        "not_established": "physical matching values, total residual bounds, UV realization, gravitational positivity or quantum bounce",
    }


if __name__ == "__main__":
    print(json.dumps(audit(), indent=2))
