"""Read-only RATE4-D8.WITNESS source and finite-order identifiability audit.

The controls are formal EFT comparisons, not unitary UV counterexamples or
members asserted to satisfy D8's extra bounds. No parent input is filled in.
"""

import hashlib
import json
from fractions import Fraction as F

import sympy as s
from p8_match1_audit import POINTS, fraction_row, solve_fraction
from p8_match1_rate_input import REPO, exact
from p8_rate4_candidate import dot, original_borns, rational
from p8_rate4_d8 import inputs as d8_inputs
from p8_rate4_d8 import matrix
from p8_rate4_parent import local_case, positive_after
from p8_rate4_remainder import contact, local_words

D8_SHA256 = "1b21c1bd43bb68a36d30671217df9748d13862a1fb941a10b59e4b1b723710f5"
SCOPE_PATHS = tuple(
    f"problems/P8/s6/continuation/s6_{number}/{name}"
    for number, name in (
        (238, "FORMULATION.md"),
        (239, "FORMULATION.md"),
        (240, "FORMULATION.md"),
        (294, "notes/matching.md"),
        (297, "FORMULATION.md"),
        (297, "notes/source.md"),
        (297, "notes/matching.md"),
        (297, "notes/forward.md"),
        (302, "notes/matching.md"),
        (336, "notes/matching.md"),
    )
)
R5_TERMS = {
    (1, 1): 3259,
    (2, 0): -18645,
    (0, 1): 2463528,
    (1, 0): -152808,
    (0, 0): 21316224,
}


def inputs():
    manifest = d8_inputs()
    path = REPO / "scripts/p8_rate4_d8.py"
    assert hashlib.sha256(path.read_bytes()).hexdigest() == D8_SHA256
    manifest[str(path.relative_to(REPO))] = D8_SHA256
    # These are actual retained sources, not a claim that their text proves new bounds.
    assert all(path in manifest for path in SCOPE_PATHS)
    return manifest


def blind_tail(*, preserved_degree, scale):
    """A finite polynomial beyond ANY requested finite S/U weighted jet.

    scale is an explicit formal comparison parameter, not a physical value.
    k>=3 puts every word in the D8 tail; 2k>preserved_degree preserves the jet.
    """
    if type(preserved_degree) is not int or preserved_degree < 0:
        raise ValueError("Require a nonnegative integer weighted degree")
    scale = exact(scale)
    k = max(3, preserved_degree // 2 + 1)
    return k, {(i + k, j): scale * c for (i, j), c in R5_TERMS.items()}


def norm64(terms):
    return sum(abs(c) * 64 ** (2 * i + 3 * j) for (i, j), c in terms.items())


def evaluate(terms, S, U):
    return sum(c * S**i * U**j for (i, j), c in terms.items())


def forward_coefficient(terms):
    return sum(2 * i * 8 ** (i - 1) * c for (i, j), c in terms.items() if j == 0)


def pole_shapes(point, mass_ratio):
    n = exact(mass_ratio)
    if n < 32:
        raise ValueError("Require n>=32")
    return (
        sum(1 / (n - a) for a in point),
        sum(1 / (n - a) ** 2 for a in point),
        fraction_row(point, n)[3],
    )


def nuisance_compensation(*, mass_ratio, direction):
    """Local fit that removes the four samples of one supplied Q/J/T direction."""
    if direction not in ("Q", "J", "T"):
        raise ValueError("Declare Q, J or T")
    index = ("Q", "J", "T").index(direction)
    values = tuple(pole_shapes(q, mass_ratio)[index] for q in POINTS)
    return tuple(solve_fraction(matrix(), tuple(-value for value in values)))


def audit():
    before = inputs()
    n0, g, kappa, lam, _ = original_borns()
    S, U, v = s.symbols("S U v")
    base = 3259 * S * U - 18645 * S**2 + 2463528 * U - 152808 * S + 21316224
    assert (
        s.Poly(base, S, U).terms()
        == s.Poly(sum(c * S**i * U**j for (i, j), c in R5_TERMS.items()), S, U).terms()
    )
    base_norm = norm64(R5_TERMS)
    assert base_norm == 4458582098560
    assert (
        s.expand(base.subs({S: 8 + 2 * v * v, U: 0}))
        == 18900480 - 902256 * v**2 - 74580 * v**4
    )

    forward_checks = sample_checks = normalized_tail_controls = 0
    supremum_controls = 0
    for degree in (4, 8, 16, 32, 128):
        k, terms = blind_tail(preserved_degree=degree, scale=1)
        assert all(2 * i + 3 * j > degree and 2 * i + 3 * j >= 5 for i, j in terms)
        assert norm64(terms) == 64 ** (2 * k) * base_norm
        expected = 8**k * (4725120 * k - 902256)
        assert expected > 0
        literal = s.expand(
            (8 + 2 * v * v) ** k * base.subs({S: 8 + 2 * v * v, U: 0})
        ).coeff(v, 2)
        assert forward_coefficient(terms) == literal == expected
        forward_checks += 1
        for q in POINTS:
            sample_S, sample_U = local_words(q)[1:3]
            assert evaluate(terms, sample_S, sample_U) == sample_S**k * contact(q) == 0
            sample_checks += 1
        # Equal first-order data and exact finite jets, but arbitrary projected shifts.
        for sign in (-1, 1):
            _, scaled = blind_tail(
                preserved_degree=degree, scale=sign * 8 * lam / expected
            )
            assert forward_coefficient(scaled) == sign * 8 * lam
            assert norm64(scaled) > 2000 * lam
            forward_checks += 1
        # Analyticity alone does not fix a SUPREMUM: this polynomial is entire.
        point_value = evaluate(terms, F(128**2), F(0))
        assert point_value != 0
        _, enlarged = blind_tail(
            preserved_degree=degree, scale=46000 * lam / abs(point_value)
        )
        assert abs(evaluate(enlarged, F(128**2), F(0))) == 46000 * lam > 23000 * lam
        supremum_controls += 1
        # Conversely, failure of a strong tail norm need not mean a large b20.
        _, norm_control = blind_tail(
            preserved_degree=degree, scale=4000 * lam / norm64(terms)
        )
        assert norm64(norm_control) == 4000 * lam
        assert abs(forward_coefficient(norm_control)) < lam / 1000
        normalized_tail_controls += 1

    # Symbolic Q/J/T compensations and nonzero subtracted forward sensitivities.
    n = s.Symbol("n", positive=True)
    _, weights = local_case(4)
    symbolic_weights = tuple(map(rational, weights))
    sampled_Q = [sum(1 / (n - a) for a in q) for q in POINTS]
    sampled_J = [sum(1 / (n - a) ** 2 for a in q) for q in POINTS]
    KQ = s.factor(2 / (n - 2) ** 3 - dot(symbolic_weights, sampled_Q))
    KJ = s.factor(6 / (n - 2) ** 4 - dot(symbolic_weights, sampled_J))
    KT = -dot(weights, [fraction_row(q, F(32))[3] for q in POINTS])
    assert s.factor(KJ + s.diff(KQ, n)) == 0
    positive_after(-KQ, n, 32)
    positive_after(-KJ, n, 32)
    assert KT == F(536987, 2346480) > 0
    assert s.limit(n**6 * KQ, n, s.oo) == -s.Rational(2255640, 3259)
    assert s.limit(n**7 * KJ, n, s.oo) == -s.Rational(13533840, 3259)
    null_equalities = projection_comparisons = 0
    for mass in (F(32), F(64), F(1024), n0):
        for index, direction in enumerate(("Q", "J", "T")):
            local = nuisance_compensation(mass_ratio=mass, direction=direction)
            for row, q in zip(matrix(), POINTS, strict=True):
                assert dot(row, local) + pole_shapes(q, mass)[index] == 0
                null_equalities += 1
            functional = (2 / (mass - 2) ** 3, 6 / (mass - 2) ** 4, F(0))[index]
            actual = functional + dot((0, 2, 0, 32), local)
            symbolic_value = (
                KQ.subs(n, rational(mass)),
                KJ.subs(n, rational(mass)),
                rational(KT),
            )[index]
            assert rational(actual) == symbolic_value
            projection_comparisons += 1
    # Moving a nuisance coordinate does not require an analytic-tail addition.
    cap_controls = 0
    for cap, direction in zip(
        (g * g / 4, g * g * n0 / 4, F(1, 4 * kappa)), ("Q", "J", "T"), strict=True
    ):
        local = nuisance_compensation(mass_ratio=n0, direction=direction)
        for sign in (-1, 1):
            scale = sign * 2 * cap
            for row, q in zip(matrix(), POINTS, strict=True):
                i = ("Q", "J", "T").index(direction)
                assert scale * (dot(row, local) + pole_shapes(q, n0)[i]) == 0
            cap_controls += 1
        assert (
            4 * cap > 2 * cap
        )  # the two shifts cannot both fit an interval of width 2*cap

    # Laurent coordinates show independence modulo the fitted polynomial.
    z, t = s.symbols("z t")
    heavy_channels = (z, s.Integer(-2), 6 - z)
    Q = sum(1 / (n - a) for a in heavy_channels)
    J = sum(1 / (n - a) ** 2 for a in heavy_channels)
    assert s.limit((z - n) * Q, z, n) == -1
    assert s.limit((z - n) ** 2 * J, z, n) == 1
    assert s.limit(s.diff(s.cancel((z - n) ** 2 * J), z), z, n) == 0
    channels = (s.Integer(8), t, -4 - t)
    T = sum(
        (2 - 2 * a - channels[(i + 1) % 3] * channels[(i + 2) % 3]) / a
        for i, a in enumerate(channels)
    )
    assert s.limit(t * T, t, 0) == 34

    # Fixing the total flat residue does not fix c_RH separately.
    H = sum((a + 2) / (n - a) for a in heavy_channels)
    assert s.factor(H - (n + 2) * Q + 3) == 0
    delta_B = s.Symbol("delta_B", real=True)
    assert s.factor(delta_B * H - (n + 2) * delta_B * Q + 3 * delta_B) == 0

    rejected = 0
    malformed = (
        {"preserved_degree": -1, "scale": 1},
        {"preserved_degree": True, "scale": 1},
        {"preserved_degree": 4.0, "scale": 1},
        {"preserved_degree": 4, "scale": None},
        {"preserved_degree": 4, "scale": 1.0},
        {"preserved_degree": 4},
    )
    for packet in malformed:
        try:
            blind_tail(**packet)
        except (ValueError, TypeError):
            rejected += 1
        else:
            raise AssertionError("Invalid formal comparison accepted")
    for packet in (
        {"mass_ratio": 16, "direction": "Q"},
        {"mass_ratio": 32.0, "direction": "J"},
        {"mass_ratio": 32, "direction": "unknown"},
        {"mass_ratio": 32},
    ):
        try:
            nuisance_compensation(**packet)
        except (ValueError, TypeError):
            rejected += 1
        else:
            raise AssertionError("Invalid nuisance direction accepted")
    assert inputs() == before
    return {
        "milestone": "RATE4-D8.WITNESS",
        "outcome": "RETAINED_SOURCE_ROUTE_INSUFFICIENT_NOT_A_PHYSICAL_EXCLUSION",
        "physical_witness_established": False,
        "protected_input_files": len(before),
        "direct_scope_sources": len(SCOPE_PATHS),
        "arbitrary_prefix_obstruction": "S^k R5, k>=3, 2k>preserved_degree; a written all-k argument",
        "tested_prefix_degrees": [4, 8, 16, 32, 128],
        "forward_coefficient_cross_checks": forward_checks,
        "tail_sample_null_equalities": sample_checks,
        "analytic_domain_without_supremum_controls": supremum_controls,
        "large_norm_small_projection_controls": normalized_tail_controls,
        "nuisance_rate_null_equalities": null_equalities,
        "nuisance_Fraction_SymPy_projection_comparisons": projection_comparisons,
        "nuisance_halfline_sign_proofs": 2,
        "independent_cap_shift_controls": cap_controls,
        "laurent_coordinate_identities": 4,
        "curvature_residue_compensation_identity": "delta_B H-(n+2)delta_B Q+3delta_B=0",
        "invalid_inputs_rejected": rejected,
        "scope": "Finite-order local EFT comparisons and source ownership only. No unitary UV counterexample, violation by a verified D8 member, all-order indistinguishability, supplied parent bound or P8 closure.",
    }


if __name__ == "__main__":
    print(json.dumps(audit(), indent=2))
