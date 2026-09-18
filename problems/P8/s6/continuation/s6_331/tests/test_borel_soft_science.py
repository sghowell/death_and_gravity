"""Independent original-current, cloud-logarithm and cutoff controls."""

from itertools import product

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_affine_borel_soft_conversion import (
    audit,
    current,
    cutoff,
    moments,
    source,
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
    "name,call,args", audit.bad_cases(), ids=[row[0] for row in audit.bad_cases()]
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
        (current.index_increment_upper, (bad,)),
        (current.conversion_increment_upper, (bad,)),
        (moments.parameters, (bad, x)),
        (moments.parameters, (1, bad)),
        (cutoff.index_mean_cutoff_upper, (1, x, bad)),
        (cutoff.conversion_mean_cutoff_upper, (bad, x, eta)),
    ):
        with pytest.raises((TypeError, ValueError)):
            call(*args)


@pytest.mark.parametrize("bad", (-1, -s.Rational(1, 100), s.Rational(1, 7), 2))
def test_positive_energy_domain(bad):
    for call in (
        current.index_increment_upper,
        current.conversion_increment_upper,
        current.regulated_increment_upper,
    ):
        with pytest.raises(ValueError):
            call(bad)


@pytest.mark.parametrize("bad", (-1, s.Rational(11, 10), 2))
def test_leading_index_domain(bad):
    with pytest.raises(ValueError):
        moments.parameters(bad, s.Rational(1, 8))


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
        cutoff.index_mean_cutoff_upper(1, x, eta)
    with pytest.raises(ValueError):
        cutoff.conversion_mean_cutoff_upper(1, x, eta)


@pytest.mark.parametrize(
    "a,x",
    tuple(
        product(
            (s.S.Zero, s.Rational(1, 10**800), s.Rational(1, 5), s.S.One),
            (s.Rational(1, 10**60), s.Rational(1, 128), s.Rational(1, 8)),
        )
    ),
)
def test_exact_conditioned_logarithmic_moments(a, x):
    m1, mlog, mrem = moments.weighted_log_moments(a, x)
    assert m1 == a * x / (a + 1)
    assert s.simplify(mlog - m1 * (1 - s.log(x) + 1 / (a + 1))) == 0
    assert (
        s.simplify(mrem - m1 * (-s.log(x) + s.polygamma(0, a + 2) + s.EulerGamma)) == 0
    )
    assert moments.conditional_index_change_upper(a, x) == 1400 * m1 / source.KAPPA
    assert (
        moments.conditional_conversion_change_upper(a, x) == 10000 * mlog / source.KAPPA
    )
    assert (
        moments.conditional_soft_transfer_tv_upper(a, x)
        == (10000 * mlog + 1400 * mrem) / source.KAPPA
    )
    assert (
        moments.conditional_regulator_dominator_upper(a, x)
        == (15000 * mlog + 1400 * mrem) / source.KAPPA
    )
    assert all(not v.has(s.Float) for v in (m1, mlog, mrem))
    if a == 0:
        assert m1 == mlog == mrem == 0


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
        moments.original_soft_transfer_tv_upper(x)
        == 19 * x * (1 - s.log(x)) / s.Integer(10) ** 1597
    )
    assert (
        cutoff.original_index_mean_cutoff_upper(x, eta)
        == 4 * eta / s.Integer(10) ** 1597
    )
    assert (
        cutoff.original_conversion_mean_cutoff_upper(x, eta)
        == 34 * eta * (1 - s.log(eta)) / s.Integer(10) ** 1597
    )
    assert source.KAPPA == s.Integer(10) ** 800
    assert source.HEAVY_MASS2 == s.Integer(10) ** 200 / 512 + 2


@pytest.mark.parametrize(
    "t", (s.S.Zero, s.Rational(1, 10**60), s.Rational(1, 128), s.Rational(1, 8))
)
def test_positive_increment_moduli_and_zero(t):
    assert current.index_increment_upper(t) == 1700 * t / source.KAPPA
    if t:
        assert (
            current.conversion_increment_upper(t)
            == 11000 * t * (1 - s.log(t)) / source.KAPPA
        )
        assert (
            current.regulated_increment_upper(t)
            == 18000 * t * (1 - s.log(t)) / source.KAPPA
        )
    else:
        assert (
            current.conversion_increment_upper(t)
            == current.regulated_increment_upper(t)
            == 0
        )


@pytest.mark.parametrize("row", range(36))
def test_original_complete_fixed_ball_calibrations(row):
    name, step, ni, rad, t = current.calibration()[row]
    assert name in ("born", "two", "split") and ni in (0, 1)
    assert rad in (s.S.Zero, s.Rational(3, 5), s.S.One)
    assert step in (s.Rational(1, 1000), s.Rational(1, 10**20))
    assert 0 < t < s.Rational(1, 8)


@pytest.mark.parametrize(
    "a,x",
    tuple(
        product(
            (s.Rational(1, 10), s.Rational(1, 2), s.S.One),
            (s.Rational(1, 128), s.Rational(1, 8)),
        )
    ),
)
def test_independent_conditioned_logarithm_integrals(a, x):
    with mp.workdps(100):
        aa, xx = mp.mpf(str(a.evalf(110))), mp.mpf(str(x.evalf(110)))
        m1, mlog, mrem = moments.weighted_log_moments(a, x)
        actual1 = aa * xx * mp.quad(lambda u: u**aa, [0, 1])
        actuallog = (
            aa
            * xx
            * mp.quad(lambda u: u**aa * (1 - mp.log(xx * u)) if u else 0, [0, 1])
        )
        actualrem = (
            aa
            * xx
            * mp.quad(lambda u: u**aa * (-mp.log(xx * (1 - u))) if u < 1 else 0, [0, 1])
        )
        for actual, exact in ((actual1, m1), (actuallog, mlog), (actualrem, mrem)):
            assert abs(actual - mp.mpf(str(exact.evalf(110)))) < mp.mpf("1e-75")
        assert mp.digamma(aa + 2) + mp.euler <= mp.mpf(3) / 2
        assert 10000 * actuallog + 1400 * actualrem < 23000 * aa * xx * (1 - mp.log(xx))
        assert 15000 * actuallog + 1400 * actualrem < 33000 * aa * xx * (1 - mp.log(xx))


@pytest.mark.parametrize(
    "x",
    (s.Rational(1, 128), s.Rational(1, 10**10000)),
    ids=("ordinary", "extremely_small"),
)
def test_rare_index_and_extremely_small_cut_not_expanded(x):
    a = s.Rational(1, 10**800)
    m1, mlog, mrem = moments.weighted_log_moments(a, x)
    with mp.workdps(900):
        aa = mp.mpf(1) / mp.mpf(10) ** 800
        xx = mp.mpf(str(x.evalf(910)))
        exact_log_ratio = 1 - mp.log(xx) + 1 / (aa + 1)
        exact_remaining_ratio = -mp.log(xx) + mp.digamma(aa + 2) + mp.euler
        assert m1 > 0
        assert abs(mp.mpf(str((mlog / m1).evalf(910))) - exact_log_ratio) < mp.mpf(
            "1e-850"
        )
        assert abs(
            mp.mpf(str((mrem / m1).evalf(910))) - exact_remaining_ratio
        ) < mp.mpf("1e-850")


@pytest.mark.parametrize(
    "A,B,t",
    tuple(
        (A, B, t)
        for A, B in (
            (s.Integer(120000), s.Integer(64000)),
            (s.Integer(80000), s.Integer(42000)),
        )
        for t in (s.Rational(1, 128), s.Rational(1, 1024), s.Rational(1, 10**30))
    ),
)
def test_independent_split_radial_majorant_integral(A, B, t):
    # Integrate after the fourth-root substitution; the two branches meet at z0.
    with mp.workdps(100):
        aa, bb, tt = mp.mpf(int(A)), mp.mpf(int(B)), mp.mpf(str(t.evalf(110)))
        z0 = aa * tt / bb
        first = mp.quad(lambda z: 4 * bb, [0, z0])
        second = mp.quad(lambda z: 4 * aa * tt / z, [z0, 1])
        expected = 4 * aa * tt * (1 + mp.log(bb / (aa * tt)))
        assert abs(first + second - expected) < mp.mpf("1e-65")


@pytest.mark.parametrize(
    "a,x,h",
    tuple(
        product(
            (s.S.Zero, s.Rational(1, 10**800), s.Rational(1, 5), s.S.One),
            (s.Rational(1, 128), s.Rational(1, 8)),
            (s.Rational(1, 4), s.S.One),
        )
    ),
)
def test_cutoff_bounds_keep_the_rare_index(a, x, h):
    eta = h * x
    assert cutoff.index_mean_cutoff_upper(a, x, eta) == 4500 * a * eta / source.KAPPA
    assert (
        cutoff.conversion_mean_cutoff_upper(a, x, eta)
        == 42000 * a * eta * (1 - s.log(eta)) / source.KAPPA
    )
    if not a:
        assert cutoff.index_mean_cutoff_upper(a, x, eta) == 0
        assert cutoff.conversion_mean_cutoff_upper(a, x, eta) == 0


def test_tail_log_subadditivity_keeps_distinguished_emissions():
    t = s.Rational(1, 128)
    count = 4
    f = lambda r: r * (1 - s.log(r))
    difference = s.expand_log(count * f(t) - f(count * t), force=True)
    assert s.simplify(difference - count * t * s.log(count)) == 0
    assert difference > 0


def test_remaining_energy_log_is_not_the_full_resolution_log():
    a, x = s.Rational(1, 5), s.Rational(1, 8)
    m1, _, mrem = moments.weighted_log_moments(a, x)
    extra = s.simplify(mrem + m1 * s.log(x))
    assert mp.mpf(str(extra.evalf(80))) > 0


def test_scoped_reference_preserves_all_original_frontiers():
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 187
    assert audit.matching()[:-1] == audit.previous.matching()
    assert len(audit.qualifications()) == 6
    assert "outer hard dimensional/evanescent" in audit.observable()["not_established"]
    assert "remaining-energy logarithm" in audit.observable()["not_established"]
