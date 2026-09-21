"""Read-only first-order construction of the separately normalized RATE4 EFT.

The four physical rate targets are new defining conditions, not measurements.
The solver requires the complete non-matter hard and rate-conversion inputs;
it does not calculate them, certify their provenance or set them to zero.
Synthetic fixtures verify normalization algebra, not physical error bounds.
"""

import hashlib
import json
from fractions import Fraction as F
from itertools import permutations, product
from math import factorial, prod

import sympy as s
from p8_match1_audit import POINTS, fraction_row, positive_on_halfline, solve_fraction
from p8_match1_rate_input import REPO, exact, interval, projection_interval
from p8_match1_rate_input import inputs as previous_inputs

RATE_INPUT_SHA256 = "33f482aa8c3dadfe9046982805889e132e3bbe7aa3c4681b87f084ed1c0ee765"
EXTRA_PINS = {
    293: "ad67bca9916c2ec0a094f28b1783bd05e070c7fddd345ffead3d6c61bc7a64d8",
    305: "22b7a327ff4d4a2c39fc6ef4fe7435e84549def254317347e7f8abc01a362011",
}


def inputs():
    manifest = previous_inputs()
    path = REPO / "scripts/p8_match1_rate_input.py"
    assert hashlib.sha256(path.read_bytes()).hexdigest() == RATE_INPUT_SHA256
    manifest[str(path.relative_to(REPO))] = RATE_INPUT_SHA256
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


def four(values):
    result = tuple(map(exact, values))
    if len(result) != 4:
        raise ValueError("Exactly four inputs in the declared POINTS order required")
    return result


def born_vector(values):
    result = four(values)
    if any(value <= 0 for value in result):
        raise ValueError("Require the four positive unexpanded Born amplitudes")
    return result


def matrix_at(mass_ratio):
    ratio = exact(mass_ratio)
    if ratio < 32:
        raise ValueError("Require heavy/light squared-mass ratio >=32")
    return [fraction_row(point, ratio) for point in POINTS]


def target_slopes(*, born, matter_hard_real):
    """rho_i for the defining target 1+hbar*rho_i (original mass-one units).

    S239's complete matter coefficient, not a selected subset, is intended.
    This exact-arithmetic interface also accepts synthetic inputs for testing.
    It does not authenticate that the caller supplied the physical coefficient.
    """
    return tuple(
        2 * value / amplitude
        for amplitude, value in zip(
            born_vector(born), four(matter_hard_real), strict=True
        )
    )


def rate4_coefficients(*, mass_ratio, born, other_hard_real, rate_conversion):
    """Solve the four first-order RATE4 conditions for fixed complementary data.

    other_hard_real is ALL of Re(F), excluding the S239 matter reference and
    the four adjustable directions. rate_conversion is ALL of the finite
    inclusive-rate conversion from that same hard reference, at x_det=1/256.
    The return order is (A, B, C, D_N). No uncomputed input defaults to zero.
    Numerical evaluation requires inputs not yet supplied by this project.
    """
    amplitudes = born_vector(born)
    hard, conversion = four(other_hard_real), four(rate_conversion)
    rhs = [
        -value - amplitude * rate / 2
        for amplitude, value, rate in zip(amplitudes, hard, conversion, strict=True)
    ]
    return tuple(solve_fraction(matrix_at(mass_ratio), rhs))


def sample_subtraction(*, mass_ratio, samples):
    """Coordinates of the unique four-shape interpolant to real sample values."""
    return tuple(solve_fraction(matrix_at(mass_ratio), four(samples)))


def rate4_intervals(
    *, mass_ratio, born, other_hard_real_intervals, rate_conversion_intervals
):
    """Enclose the first-order solution from supplied complete rational bounds.

    Transcendental hard/conversion values can be enclosed, not rounded and
    treated as exact. Return each coordinate's marginal interval AND the
    stable direct matching projection. Recombining the four marginal boxes
    discards their correlations and can lose the projection's sharp bound.
    No omitted-order or physical provenance validation is supplied here.
    """
    amplitudes = born_vector(born)
    rows = matrix_at(mass_ratio)
    hard = tuple(map(interval, other_hard_real_intervals))
    conversion = tuple(map(interval, rate_conversion_intervals))
    if len(hard) != 4 or len(conversion) != 4:
        raise ValueError("Exactly four complete hard and conversion intervals required")
    residuals = tuple(
        (-hi - amplitude * dhi / 2, -lo - amplitude * dlo / 2)
        for amplitude, (lo, hi), (dlo, dhi) in zip(
            amplitudes, hard, conversion, strict=True
        )
    )
    coordinates = []
    for axis in range(4):
        weights = solve_fraction(
            list(zip(*rows, strict=True)), [int(j == axis) for j in range(4)]
        )
        endpoints = [
            tuple(sorted((weight * lo, weight * hi)))
            for weight, (lo, hi) in zip(weights, residuals, strict=True)
        ]
        coordinates.append(
            (sum(lo for lo, _ in endpoints), sum(hi for _, hi in endpoints))
        )
    return {
        "sample_residuals": residuals,
        "coordinate_intervals": tuple(coordinates),
        "matching_projection": projection_interval(
            mass_ratio=mass_ratio, light_mass_squared=1, residuals=residuals
        ),
    }


def first_order_conversion_bound(*, kappa):
    """Bound the full retained Born-tree conversion at x_det=1/256, delta=1.

    This combines the same-D virtual soft pole, S296's finite angular term,
    the explicit resolution change from analytic x=1, and S335's full
    47-graph real-minus-soft remainder. It is NOT an all-order rate bound
    or a bound on unspecified additional classical interactions.
    """
    coupling = exact(kappa)
    if coupling <= 0:
        raise ValueError("Require kappa>0")
    resolution, remainder_constant = F(1, 256), 4 * 10**9
    # pi^2>9, |log(1/256)|<8*(7/10), |K0|<=800, |Delta_soft|<112/kappa.
    soft = 112 + F(800, 36) * 8 * F(7, 10)
    real = (
        F(361, 324)
        * (
            128 * remainder_constant * resolution
            + remainder_constant**2 * resolution**2 / 2
            + 12288 * resolution
        )
        / 36
    )
    return (soft + real) / coupling


def dot(left, right):
    return sum(a * b for a, b in zip(left, right, strict=True))


def rational(value):
    value = F(value)
    return s.Rational(value.numerator, value.denominator)


def original_borns():
    kappa, lam, n, g = 10**800, F(1, 10**600), F(10**200, 512) + 2, F(1, 8192)
    d = n - 2
    quartic = -(g**2) * (3 / d - 2 / d**2)
    gamma = g**2 / d**4
    assert gamma == F(1024, kappa)
    values = tuple(
        quartic
        + g**2 * sum(1 / (n - a) for a in point)
        - fraction_row(point, n)[3] / kappa
        for point in POINTS
    )
    for point, amplitude, low in zip(POINTS, values, (136, 228, 344, 648), strict=True):
        assert low * lam < amplitude < (low + 1) * lam < 650 * lam
        original = (
            2 * lam * sum((a - 2) ** 2 for a in point)
            + 3 * gamma * prod(point)
            - 8 * gamma
        )
        remainder = gamma * sum(F(a - 2) ** 4 / (d - (a - 2)) for a in point)
        gravity = -fraction_row(point, n)[3] / kappa
        assert amplitude == original + remainder + gravity
        assert original > 0 and remainder > 0 and gravity > 0
    return n, g, kappa, lam, values


def audit():
    before = inputs()
    n_original, g_original, k_original, lam, borns_original = original_borns()

    n = s.Symbol("n", positive=True)
    matrix = s.Matrix(
        [
            [
                sum(a * a for a in point),
                sum((s.Integer(a) + 2) / (n - a) for a in point),
                1,
                rational(fraction_row(point, F(32))[3]),
            ]
            for point in POINTS
        ]
    )
    p = 3659 * n**4 + 309616 * n**3 + 815324 * n**2 + 30852016 * n - 93987840
    determinant = -2 * p / (5 * s.prod(n - a for a in (16, 12, 10, 8, -3, -4, -6)))
    assert s.factor(matrix.det() - determinant) == 0
    positive_on_halfline(-determinant, n)
    inverse = matrix.inv().applyfunc(s.factor)
    assert (matrix * inverse).applyfunc(s.factor) == s.eye(4)
    born_symbols = s.symbols("a0:4", positive=True)
    jacobian = s.diag(*(2 / a for a in born_symbols)) * matrix
    assert s.factor(jacobian.det() - 16 * determinant / prod(born_symbols)) == 0

    # Graph valences exhaust the original E=4,L=1 one-loop inventory. Species
    # exclusions and sector identifications are the pinned S297 written proof.
    patterns = {
        (v3, v4, v5, v6)
        for v3 in range(5)
        for v4 in range(3)
        for v5 in range(2)
        for v6 in range(2)
        if v3 + 2 * v4 + 3 * v5 + 4 * v6 == 4
    }
    assert patterns == {
        (4, 0, 0, 0),
        (2, 1, 0, 0),
        (0, 2, 0, 0),
        (1, 0, 1, 0),
        (0, 0, 0, 1),
    }

    # Combine the matter and metric virtual poles BEFORE normalizing the rate.
    am, ag, bsoft, trace = s.symbols("A_m A_G Bsoft T0")
    kk = s.Symbol("kappa", positive=True)
    combined_pole = am * bsoft / (8 * s.pi**2 * kk) - 2 * trace * bsoft / (
        16 * s.pi**2 * kk**2
    )
    assert (
        s.factor(
            combined_pole.subs(trace, -kk * ag) - (am + ag) * bsoft / (8 * s.pi**2 * kk)
        )
        == 0
    )
    # Preserve the dimensional Born derivative on BOTH sides. Its finite part
    # cancels in the inclusive conversion; dropping only one copy is incorrect.
    eps, born0, born1, k0, k1, phase1, log_det, log_ref = s.symbols(
        "eps born0 born1 K0 K1 phase1 log_det log_ref"
    )
    numerator = ((born0 + eps * born1) / born0) ** 2 * (
        s.exp(2 * eps * log_det) * (1 + phase1 * eps) * (k0 + k1 * eps)
        - s.exp(2 * eps * log_ref) * k0
    )
    assert numerator.subs(eps, 0) == 0
    finite_conversion = s.diff(numerator, eps).subs(eps, 0)
    assert (
        s.expand(finite_conversion - (k1 + phase1 * k0 + 2 * (log_det - log_ref) * k0))
        == 0
    )
    assert finite_conversion.diff(born1) == 0
    assert sum(F(7, 10) ** j / factorial(j) for j in range(6)) > 2
    assert F(1, 256) < F(1, 192)
    conversion_bound = first_order_conversion_bound(kappa=k_original)
    assert conversion_bound < F(10**13, k_original)
    projected_conversion_bound = 15 * 650 * lam * conversion_bound / 2
    assert projected_conversion_bound < lam / 10**783

    # The fixed local spectator pieces and reference-scale shift occupy the
    # already subtracted axes, not extra independent matching parameters.
    spectator_log, spectator_newton = s.symbols("spectator_log spectator_delta_kappa")
    known_spectator = s.Matrix(
        [
            -spectator_log / (640 * s.pi**2 * kk**2),
            0,
            -12 * spectator_log / (640 * s.pi**2 * kk**2),
            spectator_newton / kk**2,
        ]
    )
    assert (inverse * (matrix * known_spectator) - known_spectator).applyfunc(
        s.factor
    ) == s.zeros(4, 1)

    # Literal covariant contact vertices, retaining the mass-dependent constant.
    ss, tt, mu, alpha, beta, kappa = s.symbols("ss tt mu alpha beta kappa")
    uu = 4 * mu - ss - tt
    grams = {}
    for channel, pairs in zip(
        (ss, tt, uu),
        (((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))),
        strict=True,
    ):
        for i, j in pairs:
            grams[i, j] = grams[j, i] = (channel - 2 * mu) / 2
    literal_y2 = sum(
        grams[i, j] * grams[k, ell] for i, j, k, ell in permutations(range(4))
    )
    sigma2 = ss**2 + tt**2 + uu**2
    assert s.expand(literal_y2 - 2 * (sigma2 - 4 * mu**2)) == 0
    contact = alpha * literal_y2 / (32 * s.pi**2 * kappa**2) + 24 * (
        beta + 4 * alpha
    ) * mu**2 / (384 * s.pi**2 * kappa**2)
    assert (
        s.factor(contact - (alpha * sigma2 + beta * mu**2) / (16 * s.pi**2 * kappa**2))
        == 0
    )
    delta_kappa, newton_shape, g, c_rh = s.symbols("delta_kappa T g c_RH")
    newton = -newton_shape / (kappa + delta_kappa)
    assert s.diff(newton, delta_kappa).subs(delta_kappa, 0) == newton_shape / kappa**2
    heavy = sum((a + 2) / (n - a) for a in (ss, tt, 4 - ss - tt))
    pole_sum = sum(1 / (n - a) for a in (ss, tt, 4 - ss - tt))
    assert s.factor(heavy - (n + 2) * pole_sum + 3) == 0
    assert (
        s.factor(
            2 * g * c_rh / kappa * heavy
            - 2 * g * c_rh * (n + 2) / kappa * pole_sum
            + 6 * g * c_rh / kappa
        )
        == 0
    )

    # The switch factors give exact vanishing jets, not a finite-tube majorant.
    x = s.Symbol("X")
    denominator = x**1024 + (1 - x) ** 1024
    assert denominator.subs(x, 0) == denominator.subs(x, 1) == 1
    assert s.expand((1 - x) ** 1024 - denominator) == -(x**1024)
    # Formal matter-reference targets are positive; no total error bound follows.
    target_distance = F(2, 10**199)
    assert 1 - target_distance > 0

    # A local crossing-symmetric complementary polynomial invisible at all four
    # centers. This is a scope counterexample, not a claimed UV completion.
    sigma_samples = tuple(sum(a * a for a in point) for point in POINTS)
    invisible = s.prod(sigma2.subs(mu, 1) - value for value in sigma_samples)
    assert sigma_samples == (72, 118, 176, 328)
    assert all(invisible.subs({ss: point[0], tt: point[1]}) == 0 for point in POINTS)
    v = s.Symbol("v")
    invisible_forward = invisible.subs({ss: 2 + v, tt: 0})
    assert s.expand(invisible_forward).coeff(v, 2) == -25579520
    independent_forward = 2 * sum(
        prod(8 - value for j, value in enumerate(sigma_samples) if j != i)
        for i in range(4)
    )
    assert independent_forward == -25579520

    coefficient_comparisons = rate_checks = off_point_checks = interval_corners = 0
    for ratio in (F(32), F(64), F(1000), n_original):
        rows = matrix_at(ratio)
        amplitudes = (
            borns_original
            if ratio == n_original
            else tuple(F(20 + i) for i in range(4))
        )
        matter = tuple(
            amplitude * F(i - 2, 10000) for i, amplitude in enumerate(amplitudes)
        )
        hard = tuple(
            amplitude * F((-1) ** i * (i + 1), 7000)
            for i, amplitude in enumerate(amplitudes)
        )
        conversion = tuple(F(2 * i - 3, 9000) for i in range(4))
        coefficients = rate4_coefficients(
            mass_ratio=ratio,
            born=amplitudes,
            other_hard_real=hard,
            rate_conversion=conversion,
        )
        rhs = s.Matrix(
            [
                -rational(value + amplitude * rate / 2)
                for amplitude, value, rate in zip(
                    amplitudes, hard, conversion, strict=True
                )
            ]
        )
        symbolic = matrix.subs(n, rational(ratio)).inv() * rhs
        assert all(
            rational(a) == b for a, b in zip(coefficients, symbolic, strict=True)
        )
        coefficient_comparisons += 4
        slopes = target_slopes(born=amplitudes, matter_hard_real=matter)
        for row, amplitude, light, value, rate, slope in zip(
            rows, amplitudes, matter, hard, conversion, slopes, strict=True
        ):
            assert (
                2 * (light + value + dot(row, coefficients)) / amplitude + rate == slope
            )
            rate_checks += 1

        # A finite hard-reference shift in the four-shape span changes the solved
        # coordinates, not the total hard amplitude. Includes constant/Newton.
        shift = (F(1, 300), F(-2, 700), F(4, 1100), F(-3, 1300))
        shifted_hard = tuple(
            value + dot(row, shift) for value, row in zip(hard, rows, strict=True)
        )
        shifted_coefficients = rate4_coefficients(
            mass_ratio=ratio,
            born=amplitudes,
            other_hard_real=shifted_hard,
            rate_conversion=conversion,
        )
        assert (
            tuple(a + b for a, b in zip(shifted_coefficients, shift, strict=True))
            == coefficients
        )

        # A hard/soft reference redistribution must move BOTH parts of the rate.
        soft_shift = tuple(F(i + 1, 1000) for i in range(4))
        assert (
            rate4_coefficients(
                mass_ratio=ratio,
                born=amplitudes,
                other_hard_real=tuple(
                    value + amplitude * change / 2
                    for value, amplitude, change in zip(
                        hard, amplitudes, soft_shift, strict=True
                    )
                ),
                rate_conversion=tuple(
                    rate - change
                    for rate, change in zip(conversion, soft_shift, strict=True)
                ),
            )
            == coefficients
        )

        # Exact projector: interpolation removes the samples and every basis
        # direction; a second subtraction vanishes. All fixtures are synthetic.
        interpolant = sample_subtraction(mass_ratio=ratio, samples=hard)
        sample_remainders = tuple(
            value - dot(row, interpolant) for value, row in zip(hard, rows, strict=True)
        )
        assert sample_remainders == (0, 0, 0, 0)
        assert sample_subtraction(mass_ratio=ratio, samples=sample_remainders) == (
            0,
            0,
            0,
            0,
        )
        assert (
            sample_subtraction(
                mass_ratio=ratio, samples=tuple(dot(row, shift) for row in rows)
            )
            == shift
        )
        for point in ((14, -5, -5), (9, -2, -3), (18, -7, -7)):
            row = fraction_row(point, ratio)
            assert dot(row, shifted_coefficients) + dot(row, shift) == dot(
                row, coefficients
            )
            off_point_checks += 1

        # Absorptive parts are untouched by the real finite matching coordinates.
        imaginary = tuple(F(i + 1, 5000) for i in range(4))
        for row, value, imag in zip(rows, hard, imaginary, strict=True):
            matched = (
                rational(value)
                + s.I * rational(imag)
                + rational(dot(row, coefficients))
            )
            assert s.im(matched) == rational(imag)

        enclosure = rate4_intervals(
            mass_ratio=ratio,
            born=amplitudes,
            other_hard_real_intervals=tuple(
                (value - amplitude / F(10**9), value + amplitude / F(10**9))
                for value, amplitude in zip(hard, amplitudes, strict=True)
            ),
            rate_conversion_intervals=tuple(
                (value - F(1, 10**6), value + F(1, 10**6)) for value in conversion
            ),
        )
        for value, (lo, hi) in zip(
            coefficients, enclosure["coordinate_intervals"], strict=True
        ):
            assert lo <= value <= hi
        target = [F(2), 2 * (ratio + 2) / (ratio - 2) ** 3, F(0), F(0)]
        corner_coordinates, corner_projections = [], []
        for choices in product((0, 1), repeat=4):
            samples = [
                bounds[choice]
                for bounds, choice in zip(
                    enclosure["sample_residuals"], choices, strict=True
                )
            ]
            corner = solve_fraction(rows, samples)
            corner_coordinates.append(corner)
            corner_projections.append(dot(target, corner))
            interval_corners += 1
        assert enclosure["coordinate_intervals"] == tuple(
            (
                min(c[j] for c in corner_coordinates),
                max(c[j] for c in corner_coordinates),
            )
            for j in range(4)
        )
        assert enclosure["matching_projection"] == (
            min(corner_projections),
            max(corner_projections),
        )
        # This verifies containment without using the looser, decorrelated box
        # as the reported matching projection.
        naive_box = tuple(
            sum(
                weight * bounds[j]
                for weight, bounds in zip(
                    target, enclosure["coordinate_intervals"], strict=True
                )
            )
            for j in range(2)
        )
        assert (
            naive_box[0]
            <= enclosure["matching_projection"][0]
            <= enclosure["matching_projection"][1]
            <= naive_box[1]
        )

    # Conditional fixed-n decoupling is linear in the SUBTRACTED inputs. No
    # numerical parent bound or interchange of quantum limits is asserted.
    inv_k = s.Symbol("inverse_kappa")
    vec1, vec2 = s.Matrix([1, -2, 3, -4]), s.Matrix([-5, 6, -7, 8])
    finite_inverse = matrix.subs(n, 32).inv()
    solution = -finite_inverse * (inv_k * vec1 + inv_k**2 * vec2)
    assert solution.subs(inv_k, 0) == s.zeros(4, 1)
    assert solution.diff(inv_k).subs(inv_k, 0) == -finite_inverse * vec1

    # All four covariant coordinates have the documented amplitude normalization.
    scalar_a, heavy_b, constant_c, newton_d = map(F, (1, 2, 3, 4))
    alpha_over_pi2 = 16 * k_original**2 * scalar_a
    beta_over_pi2 = 16 * k_original**2 * constant_c
    curvature = k_original * heavy_b / (2 * g_original)
    finite_newton = k_original**2 * newton_d
    assert alpha_over_pi2 / (16 * k_original**2) == scalar_a
    assert beta_over_pi2 / (16 * k_original**2) == constant_c
    assert 2 * g_original * curvature / k_original == heavy_b
    assert finite_newton / k_original**2 == newton_d
    assert lam > 0

    valid = {
        "mass_ratio": 32,
        "born": (1, 2, 3, 4),
        "other_hard_real": (1, -2, 3, -4),
        "rate_conversion": (1, 2, 3, 4),
    }
    invalid = [
        dict(valid, mass_ratio=31),
        dict(valid, mass_ratio=32.0),
        dict(valid, born=(1, 2, 3, 0)),
        dict(valid, born=(1, 2, 3, -1)),
        dict(valid, born=(1, 2, 3)),
        dict(valid, other_hard_real=(0.0, 0, 0, 0)),
        dict(valid, rate_conversion=(False, 0, 0, 0)),
    ]
    invalid.extend(
        {key: value for key, value in valid.items() if key != missing}
        for missing in ("other_hard_real", "rate_conversion")
    )
    rejected = 0
    for kwargs in invalid:
        try:
            rate4_coefficients(**kwargs)
        except (TypeError, ValueError):
            rejected += 1
        else:
            raise AssertionError("Invalid or incomplete input accepted")
    interval_valid = {
        "mass_ratio": 32,
        "born": (1, 2, 3, 4),
        "other_hard_real_intervals": ((-1, 2),) * 4,
        "rate_conversion_intervals": ((-2, 3),) * 4,
    }
    rejected_intervals = 0
    for kwargs in (
        dict(interval_valid, other_hard_real_intervals=((2, 1),) * 4),
        dict(interval_valid, rate_conversion_intervals=((0, 1),) * 3),
        {
            key: value
            for key, value in interval_valid.items()
            if key != "rate_conversion_intervals"
        },
    ):
        try:
            rate4_intervals(**kwargs)
        except (TypeError, ValueError):
            rejected_intervals += 1
        else:
            raise AssertionError("Invalid or incomplete interval input accepted")
    assert inputs() == before
    return {
        "candidate": "QG2-H8A420-RATE4, matter-calibrated rate targets, v1",
        "status": "FIRST_ORDER_FOUR_DIRECTION_NORMALIZATION_CONSTRUCTED",
        "protected_input_files": len(before),
        "rate_jacobian": "nonzero for n>=32 and all four Born amplitudes>0",
        "symbolic_fraction_coefficient_comparisons": coefficient_comparisons,
        "synthetic_first_order_rate_equalities": rate_checks,
        "off_point_finite_reference_checks": off_point_checks,
        "synthetic_interval_corners": interval_corners,
        "target_bound_from_retained_matter_result": "abs(R_target-1)<2e-199",
        "first_order_full_born_conversion_bound": "abs(D_rate)<1e-787 at the four centers",
        "first_order_conversion_contribution_to_matching_projection": "<1e-783*lambda",
        "heavy_residue_identity": "H=(n+2)*sum(1/(n-channel))-3",
        "invisible_complementary_contact_b20_per_unit": -25579520,
        "rejected_invalid_or_incomplete_inputs": rejected,
        "rejected_invalid_or_incomplete_interval_inputs": rejected_intervals,
        "not_established": "numerical full matching, omitted-order bounds, complementary curved data, a quantum decoupling theorem, V/G/B or original P8 closure",
    }


if __name__ == "__main__":
    print(json.dumps(audit(), indent=2))
