"""Read-only conditional MATCH-1 rate-input and identifiability diagnostic.

No physical rate data or quantum-parent error bounds are supplied here.
Synthetic fixtures check the interval arithmetic, not physical matching.
"""

import hashlib
import json
from fractions import Fraction as F
from itertools import product
from pathlib import Path

import sympy as s
from p8_match1_audit import (
    POINTS,
    fraction_row,
    frozen_inputs,
    positive_on_halfline,
    solve_fraction,
)

REPO = Path(__file__).resolve().parents[1]
AUDIT_SHA256 = "bbd12fdd3b3d4f6cdf871249a3537596c65613b41f9ba09dd5e6bd50f7ff609a"
EXTRA_PINS = {
    296: "40c04cadfa1a1c542c16eed48f6df12bcfe0605ad69b39dbb70fcefaefe6b3ed",
    297: "6206f5d56959f01a8467e2c2b74ff43a2fad3cc75df0cee67e9bb8ce82f39564",
    335: "f1b5b057f0877a8722d63cb560cbe7507ca3e9e639455d121f574cd6034b3b62",
}


def inputs():
    source = REPO / "scripts/p8_match1_audit.py"
    assert hashlib.sha256(source.read_bytes()).hexdigest() == AUDIT_SHA256
    manifest = frozen_inputs()
    manifest[str(source.relative_to(REPO))] = AUDIT_SHA256
    for number, expected in EXTRA_PINS.items():
        root = REPO / "problems/P8/s6/continuation" / f"s6_{number}"
        reports = list((root / "certificates").glob("*.json"))
        assert len(reports) == 1
        report = reports[0]
        assert hashlib.sha256(report.read_bytes()).hexdigest() == expected
        manifest[str(report.relative_to(REPO))] = expected
        for name, digest in json.loads(report.read_text())["source_sha256"].items():
            path = (root / name).resolve()
            assert path.is_relative_to(root.resolve())
            assert hashlib.sha256(path.read_bytes()).hexdigest() == digest
            manifest[str(path.relative_to(REPO))] = digest
    return manifest


def exact(value):
    if type(value) not in (int, F):
        raise TypeError("Expected an exact int or Fraction; no floats or booleans")
    return F(value)


def interval(value):
    lo, hi = map(exact, value)
    if lo > hi:
        raise ValueError("Reversed interval")
    return lo, hi


def rate_residual_interval(
    *, born, corrected_rate, nonlinear_error, known_hard, shape_error
):
    """Enclose one real four-shape residual, conditional on the supplied bounds.

    corrected_rate has the selected known soft/radiation terms subtracted.
    nonlinear_error bounds ALL remaining rate terms, not only a loop square.
    known_hard encloses the known real hard amplitude in the SAME reference.
    shape_error bounds all remaining amplitude structures at this point.
    Nothing defaults to zero; this does not verify physical provenance.
    """
    born, nonlinear_error, shape_error = map(
        exact, (born, nonlinear_error, shape_error)
    )
    if born <= 0 or nonlinear_error < 0 or shape_error < 0:
        raise ValueError("Require positive Born and nonnegative complete error bounds")
    rate_lo, rate_hi = interval(corrected_rate)
    hard_lo, hard_hi = interval(known_hard)
    return (
        born * (rate_lo - 1 - nonlinear_error) / 2 - hard_hi - shape_error,
        born * (rate_hi - 1 + nonlinear_error) / 2 - hard_lo + shape_error,
    )


def projection_interval(*, mass_ratio, light_mass_squared, residuals):
    """Exact interval image of the two-direction MATCH-1 projection only."""
    x, mu = map(exact, (mass_ratio, light_mass_squared))
    if x < 32 or mu <= 0 or len(residuals) != 4:
        raise ValueError("Require x>=32, mu>0 and exactly four residual intervals")
    matrix = [fraction_row(point, x) for point in POINTS]
    weights = solve_fraction(
        list(zip(*matrix, strict=True)), [F(2), 2 * (x + 2) / (x - 2) ** 3, F(0), F(0)]
    )
    lo, hi = F(0), F(0)
    for weight, values in zip(weights, residuals, strict=True):
        a, b = interval(values)
        low, high = sorted((weight * a, weight * b))
        lo += low / mu**2
        hi += high / mu**2
    return lo, hi


def audit():
    before = inputs()
    kappa, lam = 10**800, F(1, 10**600)
    n, g = F(10**200, 512) + 2, F(1, 8192)
    d = n - 2
    quartic = -(g**2) * (3 / d - 2 / d**2)
    borns = []
    born_intervals = []
    for point in POINTS:
        matter = quartic + g**2 * sum(1 / (n - a) for a in point)
        gravity = -fraction_row(point, n)[3] / kappa
        born = matter + gravity
        leading = 2 * sum((a - 2) ** 2 for a in point)
        assert leading * lam < born < (leading + 1) * lam < 650 * lam
        # Separate SymPy evaluation cross-checks each exact Born value.
        ns, gs, ds = (s.Rational(v.numerator, v.denominator) for v in (n, g, d))
        matter_s = -(gs**2) * (3 / ds - 2 / ds**2) + gs**2 * sum(
            1 / (ns - a) for a in point
        )
        ts = sum(
            (2 - 2 * a - point[(i + 1) % 3] * point[(i + 2) % 3]) / s.Integer(a)
            for i, a in enumerate(point)
        )
        assert matter_s - ts / kappa == s.Rational(born.numerator, born.denominator)
        borns.append(born)
        born_intervals.append([leading, leading + 1])

    # An explicit sufficient error allocation, NOT bounds on actual parent errors.
    total_amplitude_error = F(650, 2) * (F(1, 10000) + F(1, 40000)) + F(1, 100)
    assert 15 * total_amplitude_error == F(243, 320) < 1
    assert F(1, 256) < F(1, 192)
    assert F(10**13, kappa) == F(1, 10**787) < F(1, 10000)

    x = s.Symbol("x", positive=True)
    matrix = s.Matrix(
        [
            [
                sum(a * a for a in point),
                sum((s.Integer(a) + 2) / (x - a) for a in point),
                1,
                sum(
                    (2 - 2 * a - point[(i + 1) % 3] * point[(i + 2) % 3]) / s.Integer(a)
                    for i, a in enumerate(point)
                ),
            ]
            for point in POINTS
        ]
    )
    p = 3659 * x**4 + 309616 * x**3 + 815324 * x**2 + 30852016 * x - 93987840
    common = s.prod(x - a for a in (16, 12, 10, 8, -3, -4, -6)) / p
    positive_on_halfline(common, x)
    constants = (
        -s.Rational(7681, 24),
        s.Rational(4595, 6),
        -s.Rational(12457, 24),
        s.Rational(293, 4),
    )
    b_weights = [constant * common for constant in constants]
    assert all(
        s.factor(v) == 0
        for v in s.Matrix([b_weights]) * matrix - s.Matrix([[0, 1, 0, 0]])
    )
    b_norm = s.Rational(10069, 6) * common
    assert sum(abs(c) for c in constants) == s.Rational(10069, 6)
    assert s.limit(b_norm / x**3, x, s.oo) == s.Rational(10069, 21954)

    comparisons = 0
    for ratio in (F(32), F(64), F(1000), n):
        numeric = [fraction_row(point, ratio) for point in POINTS]
        independent = solve_fraction(
            list(zip(*numeric, strict=True)), [F(0), F(1), F(0), F(0)]
        )
        for expected, formula in zip(independent, b_weights, strict=True):
            actual = formula.subs(x, s.Rational(ratio.numerator, ratio.denominator))
            assert expected == F(int(actual.p), int(actual.q))
            comparisons += 1
    original_b_weights = independent
    b_error_norm = sum(abs(v) for v in original_b_weights)
    signs = [1 if value > 0 else -1 for value in original_b_weights]
    assert (
        sum(v * sign for v, sign in zip(original_b_weights, signs, strict=True))
        == b_error_norm
    )
    curvature_radius = F(kappa, 2 * g) * lam / 15 * b_error_norm
    assert 9 * 10**793 < curvature_radius < 10**794
    canonical_unit_epsilon = 2 * g / (10**400 * b_error_norm)
    assert F(1, 10**996) < canonical_unit_epsilon < F(1, 10**995)

    # Synthetic data only: the complete rate identity contains nonzero nuisance
    # coordinates, an imaginary hard correction and a nonzero rate remainder.
    corners = 0
    for ratio in (F(32), F(64), n):
        mat = [fraction_row(point, ratio) for point in POINTS]
        coefficients = (F(1, 1000), F(-1, 3000), F(2, 1000), F(1, 2000))
        residuals = [
            sum(a * b for a, b in zip(row, coefficients, strict=True)) for row in mat
        ]
        intervals = []
        for i, residual in enumerate(residuals):
            born = F(10 + i)
            known = F(i - 1, 10000)
            imag = F(i + 1, 10000)
            radiation = F((-1) ** i, 100000)
            shape = F((-1) ** i, 20000)
            rate = (
                (born + known + residual + shape) ** 2 + imag**2
            ) / born**2 + radiation
            exact_square = ((known + residual + shape) ** 2 + imag**2) / born**2
            values = rate_residual_interval(
                born=born,
                corrected_rate=(rate - F(1, 10**7), rate + F(1, 10**7)),
                nonlinear_error=exact_square + abs(radiation),
                known_hard=(known - F(1, 10**6), known + F(1, 10**6)),
                shape_error=abs(shape),
            )
            assert values[0] <= residual <= values[1]
            intervals.append(values)
        enclosure = projection_interval(
            mass_ratio=ratio, light_mass_squared=F(1), residuals=intervals
        )
        expected = (
            2 * coefficients[0] + 2 * (ratio + 2) / (ratio - 2) ** 3 * coefficients[1]
        )
        assert enclosure[0] <= expected <= enclosure[1]
        numeric_weights = solve_fraction(
            list(zip(*mat, strict=True)),
            [F(2), 2 * (ratio + 2) / (ratio - 2) ** 3, F(0), F(0)],
        )
        endpoints = []
        for choices in product((0, 1), repeat=4):
            value = sum(
                w * interval_[j]
                for w, interval_, j in zip(
                    numeric_weights, intervals, choices, strict=True
                )
            )
            assert enclosure[0] <= value <= enclosure[1]
            endpoints.append(value)
            corners += 1
        assert enclosure == (min(endpoints), max(endpoints))
        assert projection_interval(
            mass_ratio=ratio, light_mass_squared=2, residuals=intervals
        ) == tuple(value / 4 for value in enclosure)

    rejected = 0
    for bad_call in (
        lambda: exact(0.1),
        lambda: exact(True),
        lambda: interval((F(1), F(0))),
        lambda: projection_interval(
            mass_ratio=31, light_mass_squared=1, residuals=[(0, 0)] * 4
        ),
        lambda: projection_interval(
            mass_ratio=32, light_mass_squared=0, residuals=[(0, 0)] * 4
        ),
        lambda: projection_interval(
            mass_ratio=32, light_mass_squared=1, residuals=[(0, 0)] * 3
        ),
        lambda: rate_residual_interval(
            born=1, corrected_rate=(1, 1), known_hard=(0, 0), shape_error=0
        ),
        lambda: rate_residual_interval(
            born=0,
            corrected_rate=(1, 1),
            nonlinear_error=0,
            known_hard=(0, 0),
            shape_error=0,
        ),
        lambda: rate_residual_interval(
            born=1,
            corrected_rate=(1, 1),
            nonlinear_error=-1,
            known_hard=(0, 0),
            shape_error=0,
        ),
        lambda: rate_residual_interval(
            born=1,
            corrected_rate=(1, 1),
            nonlinear_error=0,
            known_hard=(0, 0),
            shape_error=-1,
        ),
    ):
        try:
            bad_call()
        except (TypeError, ValueError):
            rejected += 1
        else:
            raise AssertionError("Invalid input accepted")
    assert inputs() == before
    return {
        "status": "CONDITIONAL_RATE_INPUT_CONTRACT_NOT_PHYSICAL_MATCHING",
        "protected_input_files": len(before),
        "born_over_lambda_open_intervals": born_intervals,
        "sufficient_conditional_Cpos_error_over_lambda": "243/320 < 1",
        "selected_real_remainder_bound": "<1e-787, not the total rate error",
        "curvature_fit_radius_at_previous_amplitude_error": "9e793 < delta_c_RH_radius < 1e794",
        "additional_exact_B_weight_comparisons": comparisons,
        "synthetic_interval_corners": corners,
        "rejected_invalid_inputs": rejected,
        "not_established": "physical rates, full IR conversion, complete error bounds, curved matching or a UV parent",
    }


if __name__ == "__main__":
    print(json.dumps(audit(), indent=2))
