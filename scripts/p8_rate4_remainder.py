"""Read-only RATE4 remainder identifiability and conditional-budget audit.

This is not a full hard-amplitude calculation or a positivity certificate.
An explicit real, local first-order EFT comparison leaves the four rates
unchanged but changes their proposed analytic remainder projection. Its
UV realization, higher-order matching and curved stability are NOT asserted.
"""

import hashlib
import json
from fractions import Fraction as F
from itertools import permutations, product
from math import prod

import sympy as s
from p8_match1_audit import POINTS, positive_on_halfline, solve_fraction
from p8_match1_rate_input import REPO, exact, interval
from p8_rate4_candidate import dot, matrix_at, original_borns, rational
from p8_rate4_candidate import inputs as candidate_inputs

CANDIDATE_SHA256 = "acba381c608e85c30d40185f537523030cf7afbf5b33b495159163f7e362c0de"
LOCAL_COEFFICIENTS = (21316224, -152808, 2463528, -18645, 3259)


def inputs():
    manifest = candidate_inputs()
    path = REPO / "scripts/p8_rate4_candidate.py"
    assert hashlib.sha256(path.read_bytes()).hexdigest() == CANDIDATE_SHA256
    manifest[str(path.relative_to(REPO))] = CANDIDATE_SHA256
    return manifest


def local_words(point):
    sigma, cubic = sum(a * a for a in point), prod(point)
    return (1, sigma, cubic, sigma**2, sigma * cubic)


def contact(point):
    """Ten-derivative-or-lower local comparison in original mass-one units."""
    return dot(LOCAL_COEFFICIENTS, local_words(point))


def symbolic_matching(n):
    rows = []
    for point in POINTS:
        rows.append(
            [
                sum(a * a for a in point),
                sum((s.Integer(a) + 2) / (n - a) for a in point),
                1,
                sum(
                    s.Rational(2 - 2 * a - point[(i + 1) % 3] * point[(i + 2) % 3], a)
                    for i, a in enumerate(point)
                ),
            ]
        )
    matrix = s.Matrix(rows)
    target = s.Matrix([[2, 2 * (n + 2) / (n - 2) ** 3, 0, 0]])
    weights = (target * matrix.inv()).applyfunc(s.factor)
    assert (weights * matrix - target).applyfunc(s.factor) == s.zeros(1, 4)
    return matrix, weights


def gapped_spectator_bound(*, heavy_mass_squared, proca_mass_squared, kappa):
    """Bound ONLY the S305 H/Proca nonlocal insertion on |a|<=16.

    Its fixed local S/1 and Newton pieces are removed by the RATE4
    projector. Analyticity of this gapped insertion on |v|<=1 then bounds
    its subtracted v^2 coefficient by 16 times this amplitude bound.
    The massless M1 contribution and other hard sectors are excluded.
    """
    n, m, coupling = map(exact, (heavy_mass_squared, proca_mass_squared, kappa))
    if n < 10**12 or m < 10**6 or coupling <= 0:
        raise ValueError("Require n>=1e12, m>=1e6 and kappa>0")
    spin2 = 16 * (F(1, 210) / (4 * n - 16) + F(3, 14) / (4 * m - 16))
    trace = 16 * (F(68, 35) / (4 * n - 16) + F(36, 35) / (4 * m - 16))
    assert spin2 < F(1, 10**6) and trace < F(1, 10**5)
    # Three channels, |N2|<=344, |N0|<=108 and pi^2>9.
    bound = 3 * (344 * spin2 / 144 + 108 * trace / 6912) / coupling**2
    assert bound < F(1, 10**5) / coupling**2
    return bound


def conditional_remainder_interval(
    *,
    known_projection,
    complementary_projection,
    evaluation_abs,
    radiation_abs,
    omitted_abs,
):
    """Combine supplied COMPLETE projection intervals and absolute errors.

    known_projection must already enclose L F_known - w.Re F_known(q_i)
    for an independently justified L. complementary_projection must include
    every remaining finite direction. No input defaults to zero, and this
    arithmetic does not authenticate physical completeness or provenance.
    """
    known, complementary = (
        interval(known_projection),
        interval(complementary_projection),
    )
    errors = tuple(map(exact, (evaluation_abs, radiation_abs, omitted_abs)))
    if any(value < 0 for value in errors):
        raise ValueError("Require nonnegative complete absolute error bounds")
    return known[0] + complementary[0] - sum(errors), known[1] + complementary[1] + sum(
        errors
    )


def audit():
    before = inputs()
    n_original, g, kappa, lam, borns = original_borns()
    n = s.Symbol("n", positive=True)
    matrix, weights = symbolic_matching(n)

    # Heavy residue is an exactly removed shape. A mass displacement is not.
    residue_samples = s.Matrix([sum(1 / (n - a) for a in q) for q in POINTS])
    assert (
        s.Matrix([matrix[i, 1] for i in range(4)])
        - (n + 2) * residue_samples
        + 3 * s.ones(4, 1)
    ).applyfunc(s.factor) == s.zeros(4, 1)
    assert s.factor(2 / (n - 2) ** 3 - (weights * residue_samples)[0]) == 0
    mass_samples = s.Matrix([sum(1 / (n - a) ** 2 for a in q) for q in POINTS])
    mass_projection = s.factor(6 / (n - 2) ** 4 - (weights * mass_samples)[0])
    assert s.factor(mass_projection + (weights.diff(n) * residue_samples)[0]) == 0
    cubic_projection = s.factor(-(weights * s.Matrix([prod(q) for q in POINTS]))[0])
    for expression in (
        mass_projection,
        10**6 / n**6 - mass_projection,
        -cubic_projection,
        10**5 / n + cubic_projection,
    ):
        positive_on_halfline(expression, n)
    assert s.limit(n**6 * mass_projection, n, s.oo) == s.Rational(8591792, 3659)
    assert s.limit(n * cubic_projection, n, s.oo) == -s.Rational(8591792, 10977)
    # Only a CONDITIONAL first-order displacement bound, not an imposed cap.
    capped_mass_error = g**2 * (n_original / 4) * F(10**6) / n_original**6
    assert capped_mass_error < lam / 10**388

    comparisons = 0
    for ratio in (F(32), F(64), F(1000), n_original):
        numeric = matrix_at(ratio)
        target = (F(2), 2 * (ratio + 2) / (ratio - 2) ** 3, F(0), F(0))
        independent = solve_fraction(list(zip(*numeric, strict=True)), target)
        for value, symbolic in zip(independent, weights, strict=True):
            assert rational(value) == symbolic.subs(n, rational(ratio))
            comparisons += 1
        mass = 6 / (ratio - 2) ** 4 - dot(
            independent, [sum(1 / (ratio - a) ** 2 for a in q) for q in POINTS]
        )
        cubic = -dot(independent, [prod(q) for q in POINTS])
        assert rational(mass) == mass_projection.subs(n, rational(ratio))
        assert rational(cubic) == cubic_projection.subs(n, rational(ratio))
        comparisons += 2

    # Symmetric scalar contacts with s+t+u=4: the ring is Q[S,U].
    # All weighted-degree<=4 words are 1,S,U,S^2; <=5 adds S*U.
    evaluation = s.Matrix([local_words(q) for q in POINTS])
    assert evaluation[:, :4].det() == 10011648
    assert evaluation.rank() == 4
    assert evaluation * s.Matrix(LOCAL_COEFFICIENTS) == s.zeros(4, 1)
    independent_coefficients = solve_fraction(
        [local_words(q)[:4] for q in POINTS],
        [-3259 * local_words(q)[4] for q in POINTS],
    )
    assert tuple(independent_coefficients) == LOCAL_COEFFICIENTS[:4]

    ss, tt, uu, mu, v = s.symbols("ss tt uu mu v")
    sigma, cubic = ss**2 + tt**2 + uu**2, ss * tt * uu
    polynomial = (
        3259 * sigma * cubic
        - 18645 * mu * sigma**2
        + 2463528 * mu**2 * cubic
        - 152808 * mu**3 * sigma
        + 21316224 * mu**5
    )
    for permutation in permutations((ss, tt, uu)):
        assert (
            s.expand(
                polynomial.xreplace(dict(zip((ss, tt, uu), permutation))) - polynomial
            )
            == 0
        )
    for point in POINTS:
        assert (
            s.expand(polynomial.subs(dict(zip((ss, tt, uu), [mu * a for a in point]))))
            == 0
        )
    forward = s.expand(polynomial.subs({ss: 2 * mu + v, tt: 0, uu: 2 * mu - v}))
    assert forward == 18900480 * mu**5 - 902256 * mu**3 * v**2 - 74580 * mu * v**4
    # Flat contact normalization from the literal 24 assignments to slots.
    pairs = {(0, 1): ss, (2, 3): ss, (0, 2): tt, (1, 3): tt, (0, 3): uu, (1, 2): uu}
    assignments = []
    for labels in permutations(range(4)):
        channels = [pairs[tuple(sorted((labels[0], labels[j])))] for j in (1, 2, 3)]
        assignments.append(polynomial.xreplace(dict(zip((ss, tt, uu), channels))))
    assert s.expand(sum(assignments) / 24 - polynomial) == 0

    rate_equalities = 0
    for sign in (-1, 1):
        zeta = sign * 8 * lam / 902256
        assert -902256 * zeta == -sign * 8 * lam
        for point, born in zip(POINTS, borns, strict=True):
            assert 2 * zeta * contact(point) / born == 0
            rate_equalities += 1
    assert 4 * lam - 8 * lam + 4 * lam / 10**203 < 0
    assert 4 * lam + 8 * lam > 0
    # This bounds amplitudes absolutely, NOT relative corrections or UV validity.
    contact_box_bound = (
        3259 * 768 * 4096 + 18645 * 768**2 + 2463528 * 4096 + 152808 * 768 + 21316224
    )
    assert contact_box_bound == 31478479488
    assert 8 * F(contact_box_bound, 902256) < 3 * 10**5
    assert contact((10, -2, -4)) == -37140096

    # A separate direction invisible to ANY number of t=u normalizations.
    discriminant = (ss - tt) ** 2 * (tt - uu) ** 2 * (uu - ss) ** 2
    assert discriminant.subs(uu, tt) == 0
    assert s.expand(discriminant.subs({ss: 2 * mu + v, tt: 0, uu: 2 * mu - v})) == (
        64 * mu**4 * v**2 - 32 * mu**2 * v**4 + 4 * v**6
    )
    assert discriminant.subs({ss: 10, tt: -2, uu: -4}) == 112896

    gap_bound = gapped_spectator_bound(
        heavy_mass_squared=n_original, proca_mass_squared=10**6, kappa=kappa
    )
    assert gap_bound <= gapped_spectator_bound(
        heavy_mass_squared=10**12, proca_mass_squared=10**6, kappa=kappa
    )
    assert 16 * gap_bound < F(1, 1000) / kappa**2 == lam / 10**1003
    # Independently sum the exact local S305 source tensors before P removes them.
    on_shell_channels = (ss, tt, 4 - ss - tt)
    local_sum = 0
    for i, a in enumerate(on_shell_channels):
        b, c = on_shell_channels[(i + 1) % 3], on_shell_channels[(i + 2) % 3]
        source0 = (a + 2) ** 2 / 3
        source2 = 2 - 2 * a - b * c + source0 / 2
        local_sum += source2 / 960 + source0 / 384
    assert s.expand(local_sum - (sum(a * a for a in on_shell_channels) + 12) / 640) == 0
    # S302 at delta=1 controls the sample term ONLY, not its forward germ.
    assert F(15 * 10**9, kappa**2) < lam / 10**989
    # S305 massless M1 at these centers: |log(tau)|<2, hence |M1|<6/kappa^2.
    assert F(90, kappa**2) < lam / 10**998

    supplied = {
        "known_projection": (F(-1, 8), F(1, 4)),
        "complementary_projection": (F(-1, 16), F(1, 32)),
        "evaluation_abs": F(1, 128),
        "radiation_abs": F(1, 256),
        "omitted_abs": F(1, 64),
    }
    enclosure = conditional_remainder_interval(**supplied)
    corners = 0
    for known, complement, errors in product(
        supplied["known_projection"],
        supplied["complementary_projection"],
        product((-1, 1), repeat=3),
    ):
        value = (
            known
            + complement
            + dot(
                errors,
                [
                    supplied[key]
                    for key in ("evaluation_abs", "radiation_abs", "omitted_abs")
                ],
            )
        )
        assert enclosure[0] <= value <= enclosure[1]
        corners += 1
    rejected = 0
    malformed = [
        {
            key: value
            for key, value in supplied.items()
            if key != "complementary_projection"
        },
        dict(supplied, complementary_projection=None),
        dict(supplied, complementary_projection=(1, -1)),
        dict(supplied, omitted_abs=-1),
        dict(supplied, radiation_abs=0.0),
        dict(supplied, evaluation_abs=True),
        dict(supplied, known_projection=(0,)),
    ]
    for packet in malformed:
        try:
            conditional_remainder_interval(**packet)
        except (TypeError, ValueError):
            rejected += 1
        else:
            raise AssertionError("Malformed or incomplete budget accepted")
    for packet in (
        {"heavy_mass_squared": 32, "proca_mass_squared": 10**6, "kappa": kappa},
        {"heavy_mass_squared": n_original, "proca_mass_squared": 0, "kappa": kappa},
        {"heavy_mass_squared": n_original, "proca_mass_squared": 10**6, "kappa": 0},
    ):
        try:
            gapped_spectator_bound(**packet)
        except ValueError:
            rejected += 1
        else:
            raise AssertionError("Out-of-domain gapped bound accepted")
    assert inputs() == before
    return {
        "milestone": "RATE4.REMAINDER",
        "outcome": "CONDITIONAL_ONLY_COMPLEMENTARY_PROJECTION_NOT_BOUNDED",
        "protected_input_files": len(before),
        "symbolic_fraction_comparisons": comparisons,
        "halfline_sensitivity_bounds": 4,
        "contact_only_basis_rank_through_eight_derivatives": 4,
        "contact_only_basis_dimension_through_ten_derivatives": 5,
        "contact_coefficients_in_1_S_U_S2_SU_order": LOCAL_COEFFICIENTS,
        "contact_b20_per_unit_at_mu1": -902256,
        "contact_only_lambda_budget_requires": "abs(zeta)<=lambda/902256",
        "literal_contact_assignments": len(assignments),
        "first_order_rate_null_equalities": rate_equalities,
        "gapped_spectator_subtracted_projection_bound": "<1e-1003*lambda",
        "conditional_first_order_mass_shift_bound": "<1e-388*lambda IF abs(delta_n)<=n/4",
        "metric_sample_projection_bound_only": "<1e-989*lambda",
        "massless_M1_sample_projection_bound_only": "<1e-998*lambda",
        "synthetic_complete_budget_corners": corners,
        "rejected_invalid_or_incomplete_inputs": rejected,
        "not_established": "full hard bound, admissible full V/G functional, physical complementary matching, higher-order bounds, a UV realization, curved/state control or P8 closure",
    }


if __name__ == "__main__":
    print(json.dumps(audit(), indent=2))
