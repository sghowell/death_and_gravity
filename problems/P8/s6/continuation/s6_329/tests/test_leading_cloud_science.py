"""Exact rare-index moments and independent leading-reference controls."""

from itertools import product

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_affine_continuous_logarithmic_coefficient import kernels
from p8_vacuum_affine_leading_cloud_logarithmic_coefficient import (
    audit,
    cloud,
    moments,
    source,
)
from p8_vacuum_affine_radiative_logarithmic_soft_coefficient import radiative

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
    for call, args in (
        (cloud.index, (bad,)),
        (cloud.infrared_cutoff, (bad,)),
        (cloud.calorimetric_probability, (bad, s.Rational(1, 8))),
        (moments.conditional_energy_moment, (1, s.Rational(1, 8), bad)),
        (moments.original_physical_mean_upper, (bad,)),
    ):
        with pytest.raises((TypeError, ValueError)):
            call(*args)


@pytest.mark.parametrize("bad", (-1, -s.Rational(1, 10)))
def test_nonnegative_index_domain(bad):
    with pytest.raises(ValueError):
        cloud.index(bad)


@pytest.mark.parametrize("bad", (0, -1, 2))
def test_sharp_infrared_cutoff_domain(bad):
    with pytest.raises(ValueError):
        cloud.infrared_cutoff(bad)


@pytest.mark.parametrize("bad", (0, -1, s.Rational(1, 7), 2))
def test_calorimetric_resolution_domain(bad):
    with pytest.raises(ValueError):
        cloud.calorimetric_probability(1, bad)
    with pytest.raises(ValueError):
        moments.original_physical_second_moment_upper(bad)


@pytest.mark.parametrize("bad", (0, -1))
def test_positive_moment_power_domain(bad):
    with pytest.raises(ValueError):
        moments.conditional_energy_moment(1, s.Rational(1, 8), bad)


@pytest.mark.parametrize(
    "a,x,p",
    tuple(
        product(
            (0, s.Rational(1, 10**800), s.Rational(1, 5), s.Integer(2)),
            (s.Rational(1, 10**60), s.Rational(1, 128), s.Rational(1, 8)),
            (s.Rational(1, 2), s.S.One, s.Integer(2), s.Integer(3)),
        )
    ),
)
def test_exact_conditioned_moments_including_original_rare_scale(a, x, p):
    assert moments.conditional_energy_moment(a, x, p) == a * x**p / (a + p)
    assert moments.normalized_mean_change_upper(a, x) == 23667 * a * x / (a + 1)
    assert moments.normalized_second_moment_upper(a, x) == 23667**2 * a * x * x / (
        a + 2
    )
    if a:
        assert moments.conditional_energy_moment(a, x, p) > 0
    else:
        assert cloud.calorimetric_probability(a, x) == 1


@pytest.mark.parametrize(
    "x", (s.Rational(1, 10**60), s.Rational(1, 128), s.Rational(1, 8))
)
def test_original_physical_kappa_powers(x):
    assert moments.original_physical_mean_upper(x) == 19 * x / s.Integer(10) ** 1997
    assert (
        moments.original_physical_second_moment_upper(x)
        == 23 * x * x / s.Integer(10) ** 3193
    )
    assert source.KAPPA == s.Integer(10) ** 800
    assert source.HEAVY_MASS2 == s.Integer(10) ** 200 / 512 + 2


@pytest.mark.parametrize(
    "a,eta",
    tuple(
        product(
            (0, s.Rational(1, 10**800), s.Rational(1, 5), s.Integer(2)),
            (s.S.One, s.Rational(1, 2), s.Rational(1, 1024), s.Rational(1, 10**30)),
        )
    ),
)
def test_poisson_tail_energy_and_count_are_distinct(a, eta):
    mean, second = cloud.lost_energy_moments(a, eta)
    assert mean == a * eta
    assert second - mean * mean == a * eta * eta / 2
    assert cloud.finite_cutoff_count_mean(a, eta) == -a * s.log(eta)
    if a:
        assert second > mean * mean
    else:
        assert mean == second == 0


@pytest.mark.parametrize(
    "a,p",
    tuple(
        product(
            (s.Rational(1, 5), s.S.One, s.Integer(2)),
            (s.Rational(1, 2), s.S.One, s.Integer(2)),
        )
    ),
)
def test_independent_inverse_conditional_CDF_quadrature(a, p):
    with mp.workdps(100):
        aa, pp = mp.mpf(str(a.evalf(110))), mp.mpf(str(p.evalf(110)))
        x = mp.mpf(1) / 8
        actual = mp.quad(lambda u: (x * u ** (1 / aa)) ** pp, [0, 1])
        expected = aa * x**pp / (aa + pp)
        assert abs(actual - expected) < mp.mpf("1e-70")


@pytest.mark.parametrize(
    "z,eta",
    tuple(
        product(
            (s.Rational(1, 2), s.S.One, s.Integer(16)),
            (s.Rational(1, 2), s.Rational(1, 1024), s.Rational(1, 10**30)),
        )
    ),
)
def test_independent_missing_low_energy_Laplace_exponent_bound(z, eta):
    with mp.workdps(100):
        zz, ee = mp.mpf(str(z.evalf(110))), mp.mpf(str(eta.evalf(110)))
        v = zz * ee
        actual = mp.euler + mp.log(v) + mp.e1(v)
        integrated = mp.quad(lambda t: -mp.expm1(-v * t) / t if t else v, [0, 1])
        assert abs(actual - integrated) < mp.mpf("1e-70")
        assert 0 < actual < v


@pytest.mark.parametrize(
    "name,pol",
    tuple(product(("born", "one", "two", "split", "four"), ("plus", "complex"))),
)
def test_mark_is_the_original_complete_complex_coefficient(name, pol):
    E, states = radiative.calibration_states()
    u = s.Matrix([0, s.Rational(4, 5), s.Rational(3, 5)])
    n, p, c = kernels.frame_x(0, "y")
    A = p if pol == "plus" else (p + s.I * c) / s.sqrt(2)
    value = kernels.coefficient(E, states[name], u, n, A)
    born = kernels.coefficient(E, [], u, n, A)
    R = sum(q[0] for q in states[name])
    if R:
        assert abs(s.N(value - born, 90)) < 23667 * R
    else:
        assert value == born


def test_nested_total_cut_event_not_individual_particle_cuts():
    energies = (
        s.Rational(1, 16),
        s.Rational(1, 32),
        s.Rational(1, 64),
        s.Rational(1, 128),
    )
    cut = s.Rational(1, 16)
    assert all(w <= cut for w in energies) and sum(energies) > cut
    indicators = [sum(energies[:n]) <= cut for n in range(1, len(energies) + 1)]
    assert indicators == [True, False, False, False]


def test_conditioned_second_moment_is_not_the_square_of_its_mean():
    a, x = s.symbols("a x", positive=True)
    variance = s.factor(a * x * x / (a + 2) - (a * x / (a + 1)) ** 2)
    assert variance == a * x * x / ((a + 1) ** 2 * (a + 2))
    assert variance > 0


def test_scoped_reference_preserves_all_original_frontiers():
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 185
    assert audit.matching()[:-1] == audit.previous.matching()
    assert len(audit.qualifications()) == 6
    assert "interacting quantum state" in audit.observable()["not_established"]
    assert "dimensional/evanescent conversion" in audit.observable()["not_established"]
