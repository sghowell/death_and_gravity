"""Conditional RATE4-D8 construction; no physical parent bounds are supplied.

The user approved developing this separately named bounded-EFT candidate.
Every physical input below is required explicitly. Synthetic fixtures and
arithmetic enclosures do not establish its realization, V/G/B or P8 closure.
"""

import hashlib
import json
from fractions import Fraction as F
from functools import cache
from itertools import product

import sympy as s
from p8_match1_audit import POINTS, fraction_row, solve_fraction
from p8_match1_rate_input import REPO, exact, interval
from p8_rate4_candidate import born_vector, dot, four, original_borns, rational
from p8_rate4_parent import inputs as parent_inputs
from p8_rate4_parent import local_case, nonnegative, nuisance_bound, tail_constant
from p8_rate4_remainder import contact, local_words

PARENT_SHA256 = "0840a6d589327a4f4bd340e4fe7d3483f7bca6e4140efb10c7685d549165485a"
TARGET = (0, 2, 0, 32)


def inputs():
    manifest = parent_inputs()
    path = REPO / "scripts/p8_rate4_parent.py"
    assert hashlib.sha256(path.read_bytes()).hexdigest() == PARENT_SHA256
    manifest[str(path.relative_to(REPO))] = PARENT_SHA256
    return manifest


@cache
def matrix():
    return tuple(tuple(local_words(q)[:4]) for q in POINTS)


def four_intervals(values):
    result = tuple(map(interval, values))
    if len(result) != 4:
        raise ValueError("Require four complete intervals in POINTS order")
    return result


def linear_interval(weights, values):
    endpoints = [
        tuple(sorted((weight * lo, weight * hi)))
        for weight, (lo, hi) in zip(weights, values, strict=True)
    ]
    return sum(lo for lo, _ in endpoints), sum(hi for _, hi in endpoints)


def rate4_d8_coefficients(*, born, other_hard_real, rate_conversion):
    """First-order c=(c0,cS,cU,cSS), for fixed COMPLETE remaining inputs.

    other_hard_real includes K+rQ+etaJ+D_N*T+E_tail in one reference,
    excluding the fixed S239 matter term and the four adjustable contacts.
    No central value or error for any missing term is inferred here.
    """
    amplitudes = born_vector(born)
    hard, conversion = four(other_hard_real), four(rate_conversion)
    rhs = tuple(
        -value - amplitude * rate / 2
        for amplitude, value, rate in zip(amplitudes, hard, conversion, strict=True)
    )
    return tuple(solve_fraction(matrix(), rhs))


def rate4_d8_intervals(*, born, other_hard_real_intervals, rate_conversion_intervals):
    """Marginal contact intervals and their DIRECT correlated projection.

    The input boxes must enclose complete residuals. Recovering the direct
    projection from the marginal coordinate boxes loses their correlations.
    """
    amplitudes = born_vector(born)
    hard = four_intervals(other_hard_real_intervals)
    conversion = four_intervals(rate_conversion_intervals)
    residuals = tuple(
        (-hi - a * dhi / 2, -lo - a * dlo / 2)
        for a, (lo, hi), (dlo, dhi) in zip(amplitudes, hard, conversion, strict=True)
    )
    inverse = tuple(
        solve_fraction(
            list(zip(*matrix(), strict=True)), [int(i == j) for j in range(4)]
        )
        for i in range(4)
    )
    return {
        "coordinate_intervals": tuple(
            linear_interval(row, residuals) for row in inverse
        ),
        "local_projection_interval": linear_interval(local_case(4)[1], residuals),
    }


def rate_order_coefficients(*, born, known_rate_coefficient, target_rate_coefficient):
    """One FORMAL recursive step, not a higher-order calculation or bound.

    known_rate_coefficient must contain the COMPLETE finite coefficient at
    the desired order with this order's four adjustable contacts omitted.
    It includes all lower-order insertions, absorptive squares and radiation.
    Existence/finiteness of that input is an independent hypothesis.
    """
    amplitudes = born_vector(born)
    known, target = four(known_rate_coefficient), four(target_rate_coefficient)
    rhs = tuple(
        amplitude * (goal - value) / 2
        for amplitude, goal, value in zip(amplitudes, target, known, strict=True)
    )
    return tuple(solve_fraction(matrix(), rhs))


def tail_norm_from_bidisk(*, radius, analytic_supremum):
    """Sufficient Cauchy witness for N64,D8, IF its analytic premises hold.

    Require holomorphy on a neighbourhood of |S|<=R^2, |U|<=R^3 and the
    supplied supremum there. This is not an inference from real samples,
    a cutoff, a mass gap or the full nonanalytic light/gravity amplitude.
    """
    radius, bound = exact(radius), nonnegative(analytic_supremum)
    if radius <= 64:
        raise ValueError("An absolutely summable Cauchy bound requires R>64")
    t = 64 / radius
    return bound * (1 / ((1 - t**2) * (1 - t**3)) - 1 - t**2 - t**3 - t**4)


def matched_coefficient_interval(
    *,
    mass_ratio,
    born,
    reference_coefficient_interval,
    known_hard_L_interval,
    known_hard_sample_intervals,
    rate_conversion_intervals,
    remaining_rate_errors,
    higher_L_error,
    tail_norm64,
    residue_abs,
    mass_insertion_abs,
    newton_abs,
):
    """Conditional enclosure of the total coefficient, NOT a V/G verdict.

    The reference is tree plus S239 matter in the SAME admissible functional.
    Known-hard samples and L must cover the COMPLETE non-matter inventory.
    remaining_rate_errors covers all normalized rate/target/conversion and
    omitted-order errors not already included in the conversion intervals.
    higher_L_error is separate: sample bounds do not bound that functional.
    Physical provenance and existence/admissibility of L are not certified.
    """
    amplitudes = born_vector(born)
    reference = interval(reference_coefficient_interval)
    known_L = interval(known_hard_L_interval)
    known = four_intervals(known_hard_sample_intervals)
    conversion = four_intervals(rate_conversion_intervals)
    errors = tuple(map(nonnegative, four(remaining_rate_errors)))
    weights = local_case(4)[1]
    residuals = tuple(
        (klo + a * dlo / 2, khi + a * dhi / 2)
        for a, (klo, khi), (dlo, dhi) in zip(amplitudes, known, conversion, strict=True)
    )
    projected = linear_interval(weights, residuals)
    remaining = (
        sum(
            abs(w) * a * e / 2
            for w, a, e in zip(weights, amplitudes, errors, strict=True)
        )
        + nonnegative(higher_L_error)
        + tail_constant(4) * nonnegative(tail_norm64)
        + nuisance_bound(
            degree_limit=4,
            mass_ratio=mass_ratio,
            residue_abs=residue_abs,
            mass_insertion_abs=mass_insertion_abs,
            newton_abs=newton_abs,
        )
    )
    return (
        reference[0] + known_L[0] - projected[1] - remaining,
        reference[1] + known_L[1] - projected[0] + remaining,
    )


def audit():
    before = inputs()
    symbolic = s.Matrix(matrix())
    assert symbolic.det() == 10011648
    inverse_comparisons = 0
    for axis in range(4):
        independent = solve_fraction(
            list(zip(*matrix(), strict=True)), [int(j == axis) for j in range(4)]
        )
        for value, entry in zip(independent, symbolic.inv().row(axis), strict=True):
            assert rational(value) == entry
            inverse_comparisons += 1
    n0, g, kappa, lam, original = original_borns()
    assert 16 * symbolic.det() / rational(s.prod(original)) > 0
    weights = local_case(4)[1]

    # Nonzero synthetic complete inputs; not claimed to be physical evaluations.
    born = tuple(map(F, (11, 13, 17, 19)))
    hard = (F(1, 7), F(-2, 9), F(3, 11), F(-5, 13))
    conversion = (F(1, 103), F(-1, 107), F(2, 109), F(-3, 113))
    coefficients = rate4_d8_coefficients(
        born=born, other_hard_real=hard, rate_conversion=conversion
    )
    for row, a, value, rate in zip(matrix(), born, hard, conversion, strict=True):
        assert 2 * (value + dot(row, coefficients)) / a + rate == 0
    shift = (F(2, 3), F(-1, 7), F(1, 13), F(-3, 17))
    shifted = rate4_d8_coefficients(
        born=born,
        other_hard_real=tuple(
            h + dot(row, shift) for row, h in zip(matrix(), hard, strict=True)
        ),
        rate_conversion=conversion,
    )
    assert shifted == tuple(a - b for a, b in zip(coefficients, shift, strict=True))
    for q in POINTS + ((10, -2, -4), (11, -1, -6)):
        row = local_words(q)[:4]
        assert dot(row, coefficients) == dot(row, shifted) + dot(row, shift)
    # Hard/soft redistribution changes both inputs, not just one side.
    redistributed = rate4_d8_coefficients(
        born=born,
        other_hard_real=tuple(
            h + a * z / 2 for h, a, z in zip(hard, born, shift, strict=True)
        ),
        rate_conversion=tuple(d - z for d, z in zip(conversion, shift, strict=True)),
    )
    assert redistributed == coefficients

    hard_boxes = tuple((h - F(1, 100), h + F(1, 100)) for h in hard)
    rate_boxes = tuple((d - F(1, 1000), d + F(1, 1000)) for d in conversion)
    boxes = rate4_d8_intervals(
        born=born,
        other_hard_real_intervals=hard_boxes,
        rate_conversion_intervals=rate_boxes,
    )
    interval_corners = 0
    for corners in product((0, 1), repeat=8):
        actual = rate4_d8_coefficients(
            born=born,
            other_hard_real=tuple(
                box[c] for box, c in zip(hard_boxes, corners[:4], strict=True)
            ),
            rate_conversion=tuple(
                box[c] for box, c in zip(rate_boxes, corners[4:], strict=True)
            ),
        )
        for value, (lo, hi) in zip(actual, boxes["coordinate_intervals"], strict=True):
            assert lo <= value <= hi
        lo, hi = boxes["local_projection_interval"]
        assert lo <= dot(TARGET, actual) <= hi
        interval_corners += 1

    # A four-order synthetic rate retains imaginary squares and lower insertions.
    formal_equalities = 0
    absorptive_square_controls = 0
    for fixture in (1, 2):
        real, imaginary = [born], [(F(0),) * 4]
        for order in range(1, 5):
            kh = tuple(F(fixture * (i + 1), order + i + 5) for i in range(4))
            imag = tuple(F(fixture * (i + 2), order + i + 9) for i in range(4))
            extra = tuple(F((-1) ** i * fixture, 100 + order + i) for i in range(4))
            target = tuple(F(i - 2, 17) if order == 1 else F(0) for i in range(4))
            cross = tuple(
                sum(
                    real[j][i] * real[order - j][i]
                    + imaginary[j][i] * imaginary[order - j][i]
                    for j in range(1, order)
                )
                / born[i] ** 2
                for i in range(4)
            )
            known = tuple(2 * kh[i] / born[i] + cross[i] + extra[i] for i in range(4))
            fitted = rate_order_coefficients(
                born=born, known_rate_coefficient=known, target_rate_coefficient=target
            )
            real.append(tuple(kh[i] + dot(matrix()[i], fitted) for i in range(4)))
            imaginary.append(imag)
            for i in range(4):
                actual = (
                    sum(
                        real[j][i] * real[order - j][i]
                        + imaginary[j][i] * imaginary[order - j][i]
                        for j in range(order + 1)
                    )
                    / born[i] ** 2
                    + extra[i]
                )
                assert actual == target[i]
                formal_equalities += 1
                if order == 2:
                    assert imaginary[1][i] ** 2 / born[i] ** 2 > 0
                    absorptive_square_controls += 1

    t = s.Symbol("t")
    generating = 1 / ((1 - t**2) * (1 - t**3)) - 1 - t**2 - t**3 - t**4
    assert (
        s.factor(
            generating - t**5 * (1 + 2 * t - t**3 - t**4) / ((1 - t**2) * (1 - t**3))
        )
        == 0
    )
    assert tail_norm_from_bidisk(radius=128, analytic_supremum=1) == F(29, 336)
    assert tail_norm_from_bidisk(radius=128, analytic_supremum=23000 * lam) < 2000 * lam
    # Exact infinite rational analytic fixtures, with poles outside the witness bidisk.
    infinite_fixtures = 0
    for outer in (192, 256, 512):
        c, radius = 7 * lam, F(128)
        inner_ratio = F(64, outer)
        actual_norm = c * (
            1 / ((1 - inner_ratio**2) * (1 - inner_ratio**3))
            - 1
            - inner_ratio**2
            - inner_ratio**3
            - inner_ratio**4
        )
        boundary_bound = c / ((1 - (radius / outer) ** 2) * (1 - (radius / outer) ** 3))
        assert actual_norm <= tail_norm_from_bidisk(
            radius=radius, analytic_supremum=boundary_bound
        )
        infinite_fixtures += 1

    # The two known blind contacts are controlled ONLY IF the D8 tail cap holds.
    S, U, v = s.symbols("S U v")
    discr = S**3 / 2 - 20 * S**2 - 36 * S * U + 256 * S - 27 * U**2 + 320 * U - 1024
    a, b = s.symbols("a b")
    c = 4 - a - b
    assert (
        s.expand(
            discr.subs({S: a**2 + b**2 + c**2, U: a * b * c})
            - (a - b) ** 2 * (b - c) ** 2 * (c - a) ** 2
        )
        == 0
    )
    discr_norm = F(1, 2) * 64**6 + 36 * 64**5 + 27 * 64**6
    assert discr_norm == 1796 * 64**5
    for q in POINTS:
        assert contact(q) == 0
        assert discr.subs({S: sum(a * a for a in q), U: s.prod(q)}) == 0
    assert s.expand(discr.subs({S: 8 + 2 * v**2, U: 0})).coeff(v, 2) == 64
    r5_limit = F(902256 * 2000, 3259 * 64**5)
    discr_limit = F(64 * 2000, discr_norm)
    assert r5_limit == F(7048875, 13669236736) < F(1, 1000)
    assert discr_limit == F(125, 1883242496) < F(1, 10**7)
    assert F(8, 902256) * 3259 * 64**5 > 2000

    # Conditional total-coefficient enclosure, with nonzero terms of every type.
    mass, rcap, jcap, dcap, ecoef = F(64), F(1, 9), F(1, 11), F(1, 13), F(1, 17 * 64**5)
    kwargs = {
        "mass_ratio": mass,
        "born": born,
        "reference_coefficient_interval": (F(4) - F(1, 101), F(4) + F(1, 101)),
        "known_hard_L_interval": (F(2) - F(1, 103), F(2) + F(1, 103)),
        "known_hard_sample_intervals": hard_boxes,
        "rate_conversion_intervals": rate_boxes,
        "remaining_rate_errors": (F(1, 107),) * 4,
        "higher_L_error": F(1, 109),
        "tail_norm64": ecoef * 64**5,
        "residue_abs": rcap,
        "mass_insertion_abs": jcap,
        "newton_abs": dcap,
    }
    enclosure = matched_coefficient_interval(**kwargs)
    matched_corners = 0
    for signs in product((-1, 1), repeat=8):
        r, eta, dn, ec = (
            bound * sign
            for bound, sign in zip((rcap, jcap, dcap, ecoef), signs[:4], strict=True)
        )
        sampled = []
        for i, q in enumerate(POINTS):
            rem = r * sum(1 / (mass - a) for a in q) + eta * sum(
                1 / (mass - a) ** 2 for a in q
            )
            rem += dn * fraction_row(q, mass)[3] + ec * sum(a * a for a in q) * s.prod(
                q
            )
            sampled.append(
                hard[i]
                + F(signs[0], 100)
                + rem
                + born[i]
                * (conversion[i] + F(signs[1], 1000) + F(signs[4 + i], 107))
                / 2
            )
        exact_value = (
            6
            + F(signs[0], 101)
            + F(signs[1], 103)
            + F(signs[2], 109)
            + 2 * r / (mass - 2) ** 3
            + 6 * eta / (mass - 2) ** 4
            - dot(weights, sampled)
        )
        assert enclosure[0] <= exact_value <= enclosure[1]
        matched_corners += 1
    # The proposed allocation is reproduced with original Borns and no physical fill-in.
    allocation = dict(
        kwargs,
        mass_ratio=n0,
        born=original,
        reference_coefficient_interval=(4 * lam, 4 * lam),
        known_hard_L_interval=(-lam / 10, lam / 10),
        known_hard_sample_intervals=((-lam / 100, lam / 100),) * 4,
        rate_conversion_intervals=((-F(1, 20000), F(1, 20000)),) * 4,
        remaining_rate_errors=(F(1, 20000),) * 4,
        higher_L_error=lam / 4,
        tail_norm64=2000 * lam,
        residue_abs=g**2 / 4,
        mass_insertion_abs=g**2 * n0 / 4,
        newton_abs=F(1, 4 * kappa),
    )
    alo, ahi = matched_coefficient_interval(**allocation)
    assert F(7, 2) * lam < alo < 4 * lam < ahi < F(9, 2) * lam

    bad = [
        {key: value for key, value in kwargs.items() if key != missing}
        for missing in kwargs
    ] + [
        dict(kwargs, born=(1, 2, 3, 0)),
        dict(kwargs, remaining_rate_errors=(0, 0, 0, -1)),
        dict(kwargs, reference_coefficient_interval=(2, 1)),
        dict(kwargs, known_hard_L_interval=(0.0, 1)),
        dict(kwargs, known_hard_sample_intervals=((0, 1),) * 3),
        dict(kwargs, higher_L_error=True),
        dict(kwargs, tail_norm64=None),
        dict(kwargs, mass_ratio=16),
        dict(kwargs, residue_abs=-1),
        dict(kwargs, mass_insertion_abs=-1),
        dict(kwargs, newton_abs=-1),
    ]
    rejected = 0
    for packet in bad:
        try:
            matched_coefficient_interval(**packet)
        except (ValueError, TypeError):
            rejected += 1
        else:
            raise AssertionError("Malformed or incomplete coefficient packet accepted")
    for packet in (
        {"radius": 64, "analytic_supremum": 1},
        {"radius": 32, "analytic_supremum": 1},
        {"radius": 128.0, "analytic_supremum": 1},
        {"radius": 128, "analytic_supremum": -1},
        {"radius": 128},
    ):
        try:
            tail_norm_from_bidisk(**packet)
        except (ValueError, TypeError):
            rejected += 1
        else:
            raise AssertionError("Invalid Cauchy witness accepted")
    assert inputs() == before
    return {
        "candidate": "QG2-H8A420-RATE4-D8 v1, parent-bounded conditional candidate",
        "status": "CONDITIONAL_CONSTRUCTION_NOT_PHYSICAL_MATCHING",
        "development_authorized": True,
        "parent_bounds_established": False,
        "protected_input_files": len(before),
        "matching_basis": ["1", "S", "U", "S^2"],
        "jacobian": "16*10011648/prod(a_i)>0",
        "inverse_Fraction_SymPy_comparisons": inverse_comparisons,
        "complete_input_interval_corners": interval_corners,
        "synthetic_formal_order_rate_equalities": formal_equalities,
        "nonzero_absorptive_square_controls": absorptive_square_controls,
        "infinite_analytic_tail_fixtures": infinite_fixtures,
        "bidisk_R128_sufficient_norm_factor": "29/336",
        "bidisk_R128_sufficient_supremum": "23000*lambda; NOT established for this parent",
        "conditional_R5_shift_over_lambda": str(r5_limit),
        "conditional_angular_discriminant_shift_over_lambda": str(discr_limit),
        "complete_projection_fixture_corners": matched_corners,
        "rejected_invalid_or_incomplete_inputs": rejected,
        "not_established": "physical parent or tail/pole/source bounds, complete hard values, higher-order errors, admissible V/G functional, curved/state control or P8 closure",
    }


if __name__ == "__main__":
    print(json.dumps(audit(), indent=2))
