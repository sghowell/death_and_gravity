"""Independent cutoff integrals, rare-event normalization and original marks."""

from itertools import product

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_affine_quantitative_soft_cutoff import (
    audit,
    increment,
    source,
)
from p8_vacuum_affine_quantitative_soft_cutoff import (
    conditioning as c,
)

ROWS, GATES = audit.residuals(), audit.gates()


@pytest.mark.parametrize("name", tuple(ROWS))
def test_every_exact_residual(name):
    value = ROWS[name]
    assert all(
        v == 0 for v in (list(value) if isinstance(value, s.MatrixBase) else [value])
    )
    assert not s.sympify(value).has(s.Float)


@pytest.mark.parametrize("name", tuple(GATES))
def test_every_analytic_proof_gate(name):
    assert GATES[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[r[0] for r in audit.bad_cases()]
)
def test_rejected_inputs_and_scope(name, call, args):
    assert name
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize(
    "bad",
    (
        None,
        True,
        False,
        "1",
        1.0,
        s.Float(1),
        [],
        {},
        s.Matrix([1]),
        s.I,
        s.oo,
        s.nan,
        s.Symbol("a"),
    ),
)
def test_exact_apis_reject_inexact_or_unresolved_values(bad):
    x, eta = s.Rational(1, 8), s.Rational(1, 64)
    for call, args in (
        (increment.positive_increment_upper, (bad,)),
        (c.parameters, (bad, x, eta)),
        (c.parameters, (1, bad, eta)),
        (c.parameters, (1, x, bad)),
        (c.original_physical_mean_cutoff_upper, (x, bad)),
    ):
        with pytest.raises((TypeError, ValueError)):
            call(*args)


@pytest.mark.parametrize("bad", (-1, -s.Rational(1, 100), s.Rational(1, 7), 2))
def test_positive_increment_domain(bad):
    with pytest.raises(ValueError):
        increment.positive_increment_upper(bad)


@pytest.mark.parametrize("bad", (-1, s.Rational(11, 10), 2))
def test_relative_normalization_index_domain(bad):
    with pytest.raises(ValueError):
        c.parameters(bad, s.Rational(1, 8), s.Rational(1, 64))


@pytest.mark.parametrize(
    "x,eta",
    (
        (0, s.Rational(1, 64)),
        (s.Rational(1, 7), s.Rational(1, 64)),
        (s.Rational(1, 8), 0),
        (s.Rational(1, 8), -1),
        (s.Rational(1, 8), s.Rational(1, 7)),
    ),
)
def test_ordered_cutoff_domain(x, eta):
    with pytest.raises(ValueError):
        c.parameters(1, x, eta)


@pytest.mark.parametrize(
    "a,x,h",
    tuple(
        product(
            (0, s.Rational(1, 10**800), s.Rational(1, 5), s.S.One),
            (s.Rational(1, 10**60), s.Rational(1, 128), s.Rational(1, 8)),
            (s.Rational(1, 4), s.Rational(3, 4), s.S.One),
        )
    ),
)
def test_exact_tail_and_quantitative_cutoff_including_rare_index(a, x, h):
    eta = h * x
    mean, second = c.conditional_lost_energy_moments(a, x, eta)
    assert mean == a * x * (1 - (1 - h) ** (a + 1)) / (a + 1)
    assert c.conditioning_mass_loss_upper(a, x, eta) == a * h
    assert c.normalized_mean_cutoff_upper(a, x, eta) == 87334 * a * eta
    assert (
        c.normalized_second_moment_cutoff_upper(a, x, eta) == 2453486889 * a * x * eta
    )
    assert not mean.has(s.Float) and not second.has(s.Float)
    if h == 1:
        assert s.factor(second - a * x * x / (a + 2)) == 0
    if not a:
        assert mean == second == 0


@pytest.mark.parametrize(
    "x,h",
    tuple(
        product(
            (s.Rational(1, 10**60), s.Rational(1, 128), s.Rational(1, 8)),
            (s.Rational(1, 4), s.Rational(3, 4), s.S.One),
        )
    ),
)
def test_original_physical_kappa_powers(x, h):
    eta = x * h
    assert (
        c.original_physical_mean_cutoff_upper(x, eta) == 7 * eta / s.Integer(10) ** 1996
    )
    assert (
        c.original_physical_second_moment_cutoff_upper(x, eta)
        == 2 * x * eta / s.Integer(10) ** 3191
    )
    assert increment.positive_increment_upper(eta) == 4 * eta / s.Integer(10) ** 1196
    assert source.KAPPA == s.Integer(10) ** 800
    assert source.HEAVY_MASS2 == s.Integer(10) ** 200 / 512 + 2


@pytest.mark.parametrize(
    "a,h",
    tuple(
        product(
            (s.Rational(1, 5), s.Rational(1, 2), s.S.One),
            (s.Rational(1, 10), s.Rational(1, 2), s.Rational(3, 4), s.S.One),
        )
    ),
)
def test_independent_conditioned_tail_integrals(a, h):
    with mp.workdps(100):
        aa, hh = mp.mpf(str(a.evalf(110))), mp.mpf(str(h.evalf(110)))
        x = mp.mpf(1) / 8
        mean, second = c.conditional_lost_energy_moments(a, s.Rational(1, 8), h / 8)
        actual_mean = aa * x * mp.quad(lambda u: (1 - u) ** aa, [0, hh])
        diagonal = aa * x * x * mp.quad(lambda u: u * (1 - u) ** aa, [0, hh])

        def inner(u):
            upper = min(hh, 1 - u)
            if upper <= 0:
                return mp.mpf(0)
            return ((1 - u) ** (aa + 1) - (1 - u - upper) ** (aa + 1)) / (aa + 1)

        intervals = sorted({mp.mpf(0), min(hh, 1 - hh), hh})
        pair = aa * aa * x * x * mp.quad(inner, intervals)
        assert abs(actual_mean - mp.mpf(str(mean.evalf(110)))) < mp.mpf("1e-75")
        assert abs(diagonal + pair - mp.mpf(str(second.evalf(110)))) < mp.mpf("1e-75")
        assert 0 < actual_mean <= aa * x * hh
        assert diagonal + pair <= aa * (x * hh) ** 2 / 2 + (aa * x * hh) ** 2


@pytest.mark.parametrize(
    "a,y",
    tuple(
        product(
            (s.Rational(1, 10), s.Rational(1, 2), s.S.One),
            (s.S.One, s.Rational(3, 2), s.Integer(2), s.Rational(5, 2), s.Integer(3)),
        )
    ),
)
def test_independent_finite_Poisson_sectors_and_relative_normalization(a, y):
    # Here eta=x/y and y<=3: all nonzero count sectors fit in n<=2.
    # The full Bose factors and the finite zero-count atom are retained.
    with mp.workdps(100):
        aa, yy = mp.mpf(str(a.evalf(110))), mp.mpf(str(y.evalf(110)))
        pair = (
            mp.quad(lambda u: mp.log(yy - u) / u, [1, yy - 1]) if yy > 2 else mp.mpf(0)
        )
        finite_scaled = 1 + aa * mp.log(yy) + aa * aa * pair / 2
        limiting_scaled = mp.exp(-mp.euler * aa) / mp.gamma(1 + aa) * yy**aa
        ratio = limiting_scaled / finite_scaled
        assert 1 - aa / yy <= ratio <= 1
        # Extremely small x is not expanded or used as a loose error denominator.
        x = mp.mpf("1e-10000")
        eta = x / yy
        finite = eta**aa * finite_scaled
        limiting = mp.exp(-mp.euler * aa) / mp.gamma(1 + aa) * x**aa
        assert abs(limiting / finite - ratio) < mp.mpf("1e-90")
        if yy > 2:
            derivative_pair = mp.quad(lambda u: 1 / (u * (yy - u)), [1, yy - 1])
            assert abs(yy * derivative_pair - 2 * mp.log(yy - 1)) < mp.mpf("1e-90")


@pytest.mark.parametrize("h", (s.Rational(1, 100), s.Rational(3, 4), s.S.One))
def test_rare_index_tail_scale_without_rounding_it_to_zero(h):
    a = s.Rational(1, 10**800)
    mean, second = c.conditional_lost_energy_moments(a, s.Rational(1, 8), h / 8)
    with mp.workdps(900):
        aa = mp.mpf(1) / mp.mpf(10) ** 800
        hh = mp.mpf(str(h.evalf(910)))
        if hh == 1:
            ratio = mp.mpf(1) / (1 + aa)
        else:
            ratio = -mp.expm1((1 + aa) * mp.log1p(-hh)) / (1 + aa)
        assert abs(mp.mpf(str((8 * mean / a).evalf(910))) - ratio) < mp.mpf("1e-850")
        assert mp.mpf(str(second.evalf(910))) > 0


@pytest.mark.parametrize("row", range(16))
def test_original_complete_coefficient_positive_tail_calibrations(row):
    name, pol, t, slope = increment.original_calibration()[row]
    assert name in ("born", "one", "two", "four") and pol in ("plus", "complex")
    assert t in (s.Rational(1, 128), s.Rational(1, 10**20))
    assert mp.mpf(slope) < 40000
    assert isinstance(slope, str)


def test_conditioning_does_not_preserve_independent_counts():
    # Two independent pre-cut Bernoulli indicators become anticorrelated
    # under the total-count cut. The argument uses independence only before conditioning.
    states = tuple(product((0, 1), repeat=2))
    selected = tuple(row for row in states if sum(row) <= 1)
    efirst = s.Rational(sum(row[0] for row in selected), len(selected))
    esecond = s.Rational(sum(row[1] for row in selected), len(selected))
    ejoint = s.Rational(sum(row[0] * row[1] for row in selected), len(selected))
    assert ejoint - efirst * esecond == -s.Rational(1, 9)


def test_second_moment_keeps_diagonal_and_ordered_pair():
    x = s.Rational(1, 8)
    a = s.Rational(1, 5)
    mean, second = c.conditional_lost_energy_moments(a, x, x)
    diagonal = a * x * x / ((a + 1) * (a + 2))
    pair = a * a * x * x / ((a + 1) * (a + 2))
    assert second == diagonal + pair
    assert second > mean * mean and pair > 0


def test_zero_positive_increment_is_exact():
    assert increment.positive_increment_upper(0) == 0


def test_scoped_reference_preserves_all_original_frontiers():
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 186
    assert audit.matching()[:-1] == audit.previous.matching()
    assert len(audit.qualifications()) == 6
    assert "Arbitrary signed" in audit.observable()["not_established"]
    assert "dimensional/evanescent" in audit.observable()["not_established"]
