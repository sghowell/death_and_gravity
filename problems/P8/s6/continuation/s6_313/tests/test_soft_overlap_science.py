"""General overlap identities, integration budgets, exact amplitudes and scope."""

import pytest
import sympy as s
from p8_vacuum_affine_complete_two_graviton_tree import trees as e
from p8_vacuum_affine_massive_scalar_graviton_bremsstrahlung import tree as matter
from p8_vacuum_affine_two_real_soft_overlap import (
    analytic,
    audit,
    external,
    source,
)
from p8_vacuum_affine_two_real_soft_overlap import (
    subtraction as sub,
)

ROWS, GATES = audit.residuals(), audit.gates()


@pytest.mark.parametrize("name", tuple(ROWS))
def test_all_exact_residuals(name):
    value = ROWS[name]
    assert all(
        v == 0 for v in (list(value) if isinstance(value, s.MatrixBase) else [value])
    )
    assert not s.sympify(value).has(s.Float)


@pytest.mark.parametrize("name", tuple(GATES))
def test_all_written_proof_gates(name):
    assert GATES[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[r[0] for r in audit.bad_cases()]
)
def test_scope_rejections(name, call, args):
    assert name
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize(
    "first,second",
    (
        (s.Rational(1, 16), s.Rational(1, 16)),
        (s.Rational(1, 10**400), s.Rational(1, 20)),
        (s.Rational(1, 10**40), s.Rational(1, 10**60)),
    ),
)
def test_entropy_bound_is_symmetric_and_homogeneous(first, second):
    entropy = sub.entropy(first, second)
    assert s.expand_log(entropy - sub.entropy(second, first), force=True) == 0
    assert (
        s.expand_log(sub.entropy(first / 2, second / 2) - entropy / 2, force=True) == 0
    )
    assert (
        s.expand_log(
            sub.remainder_upper(first / 2, second / 2)
            - 2 * sub.remainder_upper(first, second),
            force=True,
        )
        == 0
    )
    assert not entropy.has(s.Float)
    assert sub.remainder_upper(
        first, second
    ) == analytic.REMAINDER_COEFFICIENT * entropy / (first * second)


@pytest.mark.parametrize(
    "invalid",
    (
        True,
        False,
        1.0,
        s.Float(1),
        "1",
        None,
        s.oo,
        s.I,
        s.nan,
        s.Symbol("missing"),
        0,
        -1,
        s.Rational(1, 4),
    ),
)
def test_energy_and_resolution_domains_reject_inexact_or_outside_inputs(invalid):
    for call in (sub.entropy, sub.remainder_upper):
        with pytest.raises((TypeError, ValueError)):
            call(invalid, s.Rational(1, 64))
    for call in (
        sub.require_resolution,
        sub.measure_upper,
        sub.leading_reference_relative_upper,
    ):
        with pytest.raises((TypeError, ValueError)):
            call(invalid)


@pytest.mark.parametrize(
    "first,second",
    (
        (s.Rational(1, 8), s.Rational(1, 8)),
        (s.Rational(1, 10), s.Rational(1, 10)),
        (s.Rational(1, 8), s.Rational(1, 10**40)),
    ),
)
def test_joint_total_energy_domain_is_enforced(first, second):
    with pytest.raises(ValueError):
        sub.remainder_upper(first, second)


@pytest.mark.parametrize(
    "resolution",
    (
        s.Rational(1, 8),
        s.Rational(1, 1000),
        s.Rational(1, 10**400),
    ),
)
def test_signed_measure_and_reference_bounds_stay_exact(resolution):
    result = sub.measure_upper(resolution)
    relative = sub.leading_reference_relative_upper(resolution)
    assert (
        result
        == 2 * s.Rational(1, 10**725) * resolution
        + s.Rational(1, 10**652) * resolution**2
    )
    assert relative == 2 * result / s.sqrt(resolution)
    assert result > 0 and not result.has(s.Float)
    assert relative > 0 and not relative.has(s.Float)


@pytest.mark.parametrize(
    "pair,expected",
    (
        ((0, 0), 2642),
        ((0, 1), 432),
        ((1, 0), 352),
        ((1, 1), 1764),
    ),
)
def test_every_complex_pair_angular_coefficient_budget(pair, expected):
    checks, budgets = analytic.angular_coefficients()
    assert all(value == 0 for value in checks.values())
    assert budgets[pair] == expected < 10**4


@pytest.mark.parametrize("index", range(6))
def test_entire_original_amplitude_samples_are_exact_and_not_reduced_leading_products(
    index,
):
    _, _, records = sub.original_samples()
    family, denominator, a, b, F2, B2 = records[index]
    assert family in ("simultaneous", "hierarchical")
    assert denominator in (100, 1000, 10000)
    assert 0 < a + b <= s.Rational(1, 8)
    assert not F2.has(s.Float) and not B2.has(s.Float)
    assert F2 != B2


@pytest.mark.parametrize("kind", ("C", "H", "GR"))
def test_every_double_external_sector_needs_its_nonleading_correction(kind):
    checks, gates, corrections = external.independent_tree()
    assert all(value == 0 for value in checks.values())
    assert gates[kind + "_leading_only_grouping_has_nonzero_correction"]
    assert corrections[kind] != 0


def test_nonopposite_marked_null_leg_changes_the_physical_TT_current():
    energy = s.Rational(5, 4)
    weight = s.Rational(1, 32)
    marked = e.imm([weight, 3 * weight / 5, 0, 4 * weight / 5])
    soft = e.imm([1, 0, 0, 1])
    outgoing = e.imm([0, s.Rational(3, 5), s.Rational(4, 5)])
    points, _, _ = source.recoil.momenta(energy, [marked], outgoing)
    full = matter.soft_current((*points, marked), soft)
    massive = matter.soft_current(points, soft)
    assert (full * e.ETA * soft).applyfunc(s.simplify) == e.VECTOR_ZERO
    assert (massive * e.ETA * soft + marked).applyfunc(s.simplify) == e.VECTOR_ZERO
    difference = s.simplify(e.old.pair(sub.EPS, full - massive))
    assert difference == (marked.T * sub.EPS * marked)[0] / e.old.dot(marked, soft)
    assert difference != 0


def test_general_ordering_controls_are_nonzero():
    checks, controls = external.general_identities()
    assert all(
        v == 0
        for value in checks.values()
        for v in (list(value) if isinstance(value, s.MatrixBase) else [value])
    )
    assert all(controls.values())


def test_Cauchy_coefficient_uses_the_complex_radius_twice():
    assert analytic.CAUCHY_COEFFICIENT == analytic.LINEAR_SLOPE / analytic.RADIUS**2
    assert analytic.CAUCHY_COEFFICIENT < analytic.REMAINDER_COEFFICIENT
    assert analytic.LINEAR_SLOPE < s.Rational(1, 10**350)


def test_entropy_integral_and_cutoff_limit_are_exact():
    a, b = s.symbols("a b", positive=True)
    I = (a + b) * s.log(a + b) - a * s.log(a) - b * s.log(b)
    assert s.simplify(s.diff(I, a, b) - 1 / (a + b)) == 0
    assert s.limit(I, a, 0, dir="+") == 0
    assert s.limit(I, b, 0, dir="+") == 0
    assert sub.entropy_quadratures() == (True, True)


def test_measure_difference_retains_interference():
    br, bi, rr, ri = s.symbols("br bi rr ri", real=True)
    B = br + s.I * bi
    R = rr + s.I * ri
    difference = s.expand((B + R) * s.conjugate(B + R) - B * s.conjugate(B))
    assert s.expand(difference - 2 * s.re(B * s.conjugate(R)) - R * s.conjugate(R)) == 0
    assert s.expand(difference - R * s.conjugate(R)) != 0


def test_original_source_and_unclosed_frontiers_preserved():
    assert (
        source.KAPPA == 10**800 and source.HEAVY_MASS2 == s.Rational(10**200, 512) + 2
    )
    assert source.CUBIC == s.Rational(1, 8192)
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 169
    assert len(audit.qualifications()) == 6
    assert audit.matching()[:-1] == audit.previous.matching()
    assert (
        "integrated counterterm/virtual matching"
        in audit.observable()["not_established"]
    )
