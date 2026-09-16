"""Independent whole-soft sums, Laplace inversions and order-of-limit tests."""

import itertools
import math

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_affine_calorimetric_soft_resummation import audit, gamma, limits, poisson

ROWS, GATES = audit.residuals(), audit.gates()


@pytest.mark.parametrize("name", tuple(ROWS))
def test_every_exact_identity(name):
    value = ROWS[name]
    assert all(
        x == 0 for x in (list(value) if isinstance(value, s.MatrixBase) else [value])
    )
    assert not s.sympify(value).has(s.Float)


@pytest.mark.parametrize("name", tuple(GATES))
def test_every_written_proof_gate(name):
    assert GATES[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[r[0] for r in audit.bad_cases()]
)
def test_every_scope_guard(name, call, args):
    assert name
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize("count", range(1, 8))
def test_independent_eikonal_permutation_sum(count):
    weights = tuple(s.Rational(2 * i + 1, i + 1) for i in range(count))
    result = 0
    for order in itertools.permutations(weights):
        partial, term = 0, s.S.One
        for value in order:
            partial += value
            term /= partial
        result += term
    assert result == 1 / math.prod(weights)


def poisson_expectation(rate, epsilon):
    mode = int(mp.floor(rate))
    center = mp.exp(-rate + mode * mp.log(rate) - mp.loggamma(mode + 1))
    result = center / mp.gamma(1 + 2 * epsilon * mode)
    weight = center
    for n in range(mode, 0, -1):
        weight *= n / rate
        result += weight / mp.gamma(1 + 2 * epsilon * (n - 1))
    weight, n = center, mode
    while True:
        n += 1
        weight *= rate / n
        result += weight / mp.gamma(1 + 2 * epsilon * n)
        if n > rate + 20 * mp.sqrt(rate + 1) and weight < mp.mpf("1e-75"):
            break
        assert n < 100000
    return result


def finite_sum(a, delta, epsilon, ratio):
    rate = (a + 2 * epsilon * delta) * mp.gamma(2 * epsilon) * ratio ** (2 * epsilon)
    return mp.exp(rate - a / (2 * epsilon)) * poisson_expectation(rate, epsilon)


@pytest.mark.parametrize(
    "a,delta,ratio",
    (
        (".2", ".03", ".4"),
        (".7", "-.02", ".125"),
        ("1.2", ".1", ".8"),
        (".01", "0", "1e-12"),
    ),
)
def test_finite_regulator_whole_sum_converges(a, delta, ratio):
    with mp.workdps(80):
        a, delta, ratio = map(mp.mpf, (a, delta, ratio))
        expected = mp.exp(delta - mp.euler * a) * ratio**a / mp.gamma(1 + a)
        errors = [
            abs(finite_sum(a, delta, e, ratio) / expected - 1)
            for e in map(mp.mpf, (".05", ".01", ".002"))
        ]
        assert errors[2] < errors[1] < errors[0]
        assert errors[2] < mp.mpf(".02")


@pytest.mark.parametrize("a", (".0001", ".2", "1", "5", "20"))
def test_whole_Gamma_bound(a):
    with mp.workdps(80):
        a = mp.mpf(a)
        logF = -mp.euler * a - mp.loggamma(1 + a)
        assert -(mp.pi**2) * a * a / 12 <= logF <= 0
        assert 0 <= -mp.expm1(logF) <= mp.pi**2 * a * a / 12


@pytest.mark.parametrize("a,ratio", ((".3", ".2"), (".7", ".75")))
def test_independent_finite_reference_Laplace_inversion(a, ratio):
    with mp.workdps(45):
        a, ratio = mp.mpf(a), mp.mpf(ratio)
        transform = lambda z: mp.exp(-a * (mp.euler + mp.log(z) + mp.e1(z))) / z
        actual = mp.invertlaplace(transform, ratio, method="dehoog")
        expected = mp.exp(-mp.euler * a) * ratio**a / mp.gamma(1 + a)
        assert mp.almosteq(actual, expected, rel_eps=mp.mpf("1e-25"))


def test_above_reference_energy_power_law_is_not_the_probability():
    with mp.workdps(45):
        a, ratio = mp.mpf(".7"), mp.mpf("1.4")
        transform = lambda z: mp.exp(-a * (mp.euler + mp.log(z) + mp.e1(z))) / z
        actual = mp.invertlaplace(transform, ratio, method="dehoog")
        naive = mp.exp(-mp.euler * a) * ratio**a / mp.gamma(1 + a)
        assert 0 < actual < naive
        assert abs(actual - naive) > mp.mpf(".01")


def test_two_real_total_energy_loss_wedge():
    with mp.workdps(80):
        value = mp.quad(lambda x: -mp.log1p(-x) / x, [0, mp.mpf(".5"), 1])
        assert mp.almosteq(value, mp.pi**2 / 6, rel_eps=mp.mpf("1e-65"))


@pytest.mark.parametrize("rho", (".001", ".1", "1"))
def test_joint_regulator_scaling_error_and_finite_sum(rho):
    with mp.workdps(80):
        rho, chi, A, D = map(mp.mpf, (rho, ".7", ".5", ".1"))
        exponent = A * mp.expm1(-2 * rho * chi) / (2 * rho)
        assert 0 <= exponent + A * chi <= A * rho * chi**2
        errors = []
        for k in (32, 128, 512):
            actual = finite_sum(A / k, D / k, rho / k, mp.exp(-k * chi))
            errors.append(abs(actual / mp.exp(exponent) - 1))
        assert errors[2] < errors[1] < errors[0]
        assert errors[2] < mp.mpf(".001")


@pytest.mark.parametrize(
    "args",
    ((0, 1), (1, s.Rational(1, 8)), (s.Rational(1, 10), s.Rational(1, 10**1000))),
)
def test_valid_exact_detector_domains(args):
    assert limits.require_domain(*args) == args


@pytest.mark.parametrize(
    "args",
    (
        (-1, 1),
        (1, 0),
        (1, 2),
        (1.0, 1),
        (True, 1),
        (1, s.Float(".5")),
        (s.Symbol("a"), 1),
        (1, s.oo),
    ),
)
def test_invalid_detector_domains(args):
    with pytest.raises((ValueError, TypeError)):
        limits.require_domain(*args)


@pytest.mark.parametrize("count", (-1, True, 1.0, s.Rational(1, 2), s.Symbol("N")))
def test_invalid_emission_counts(count):
    with pytest.raises((ValueError, TypeError)):
        poisson.simplex_term(count, 1)


def test_original_scope_and_unmatched_remainder_are_retained():
    assert len(audit.frontier()) == 9
    assert len(audit.matching()) == 155
    assert audit.matching()[:-1] == audit.previous.matching()
    assert len(audit.qualifications()) == 6
    assert "Uniform in multiplicity" in audit.observable()["not_established"]
    assert gamma.uniform_error_bound(10**800) == s.Rational(2, 10**1596)
