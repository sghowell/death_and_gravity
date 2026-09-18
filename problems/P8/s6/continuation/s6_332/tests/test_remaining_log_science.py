"""Independent logarithmic integration, density controls and exact API tests."""

from itertools import product

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_affine_borel_soft_conversion import current, moments
from p8_vacuum_affine_radiative_soft_state_transfer import continuity
from p8_vacuum_affine_remaining_log_cutoff import audit, density, endpoint, source, tail

ROWS, GATES = audit.residuals(), audit.gates()


@pytest.mark.parametrize("name", tuple(ROWS))
def test_every_exact_residual(name):
    value = ROWS[name]
    assert all(
        v == 0 for v in (list(value) if isinstance(value, s.MatrixBase) else [value])
    )
    assert not s.sympify(value).has(s.Float)


@pytest.mark.parametrize("name", tuple(GATES))
def test_every_written_proof_gate(name):
    assert GATES[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[row[0] for row in audit.bad_cases()]
)
def test_rejected_inputs_and_scope(name, call, args):
    assert name
    with pytest.raises((TypeError, ValueError)):
        call(*args)


CALLS = (
    density.conversion_density_cutoff_upper,
    density.common_space_cutoff_upper,
    tail.tail_remaining_log_upper,
    tail.logratio_mark_upper,
    tail.same_event_log_upper,
    endpoint.singular_layer_integral,
    endpoint.log_event_change_upper,
)


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
def test_exact_apis_reject_inexact_or_unresolved_inputs(bad):
    for call in CALLS:
        for args in (
            (bad, s.Rational(1, 8), s.Rational(1, 64)),
            (1, bad, s.Rational(1, 64)),
            (1, s.Rational(1, 8), bad),
        ):
            with pytest.raises((TypeError, ValueError)):
                call(*args)


@pytest.mark.parametrize(
    "a,x,eta",
    (
        (-1, s.Rational(1, 8), s.Rational(1, 64)),
        (s.Rational(11, 10), s.Rational(1, 8), s.Rational(1, 64)),
        (1, 0, s.Rational(1, 64)),
        (1, s.Rational(1, 7), s.Rational(1, 64)),
        (1, s.Rational(1, 8), 0),
        (1, s.Rational(1, 8), -1),
        (1, s.Rational(1, 8), s.Rational(1, 7)),
    ),
)
def test_ordered_exact_physical_domain(a, x, eta):
    for call in CALLS:
        with pytest.raises(ValueError):
            call(a, x, eta)


@pytest.mark.parametrize(
    "a,x,h",
    tuple(
        product(
            (s.S.Zero, s.Rational(1, 10**800), s.Rational(1, 5), s.S.One),
            (s.Rational(1, 10**60), s.Rational(1, 128), s.Rational(1, 8)),
            (s.Rational(1, 100), s.Rational(1, 4), s.Rational(3, 4), s.S.One),
        )
    ),
)
def test_exact_bound_budgets_and_zero_index(a, x, h):
    eta = h * x
    L = 1 - s.log(eta)
    assert (
        density.conversion_density_cutoff_upper(a, x, eta)
        == 42000 * a * eta * L / source.KAPPA
    )
    assert (
        density.common_space_cutoff_upper(a, x, eta)
        == 63000 * a * eta * L**2 / source.KAPPA
    )
    assert tail.tail_remaining_log_upper(a, x, eta) == a * eta * (2 - s.log(x))
    assert tail.logratio_mark_upper(a, x, eta) == 11200 * a * eta * L / source.KAPPA
    assert tail.same_event_log_upper(a, x, eta) == 14600 * a * eta * L / source.KAPPA
    assert (
        endpoint.log_event_change_upper(a, x, eta)
        == 6300 * a * eta * L**2 / source.KAPPA
    )
    for call in CALLS:
        assert not call(a, x, eta).has(s.Float)
        if a == 0:
            assert call(a, x, eta) == 0


@pytest.mark.parametrize(
    "t", (s.Rational(1, 1000), s.Rational(1, 10**40), s.Rational(1, 8))
)
def test_actual_frozen_physical_moduli(t):
    assert current.index_increment_upper(t) == 1700 * t / source.KAPPA
    assert (
        s.simplify(
            current.conversion_increment_upper(t)
            - 11000 * t * (1 - s.log(t)) / source.KAPPA
        )
        == 0
    )
    assert continuity.index_change_upper(t) == 1400 * t / source.KAPPA
    assert (
        s.simplify(
            continuity.finite_conversion_change_upper(t)
            - 10000 * t * (1 - s.log(t)) / source.KAPPA
        )
        == 0
    )


@pytest.mark.parametrize(
    "x,h",
    tuple(
        product(
            (s.Rational(1, 10**60), s.Rational(1, 128), s.Rational(1, 8)),
            (s.Rational(1, 4), s.S.One),
        )
    ),
)
def test_original_physical_power_and_bound(x, h):
    eta = h * x
    assert (
        density.original_common_space_cutoff_upper(x, eta)
        == 51 * eta * (1 - s.log(eta)) ** 2 / s.Integer(10) ** 1597
    )
    assert source.KAPPA == s.Integer(10) ** 800
    assert source.HEAVY_MASS2 == s.Integer(10) ** 200 / 512 + 2


@pytest.mark.parametrize(
    "a,h",
    tuple(
        product(
            (s.Rational(1, 10), s.Rational(1, 2), s.S.One),
            (s.Rational(1, 100), s.Rational(1, 4), s.Rational(3, 4), s.S.One),
        )
    ),
)
def test_independent_tail_log_integral_and_full_tail_moment(a, h):
    with mp.workdps(100):
        aa, hh = mp.mpf(str(a.evalf(110))), mp.mpf(str(h.evalf(110)))
        xx = mp.mpf(1) / 8
        H = mp.digamma(aa + 1) + mp.euler
        actual = (
            aa
            * xx
            * mp.quad(
                lambda z: (1 - z) ** aa * (-mp.log(xx * (1 - z)) + H) if z < 1 else 0,
                [0, hh],
            )
        )
        assert 0 < actual <= aa * xx * hh * (-mp.log(xx) + 2)
        if h == 1:
            _, _, exact = moments.weighted_log_moments(a, s.Rational(1, 8))
            assert abs(actual - mp.mpf(str(exact.evalf(110)))) < mp.mpf("1e-75")


@pytest.mark.parametrize(
    "astr,estr",
    tuple(product(("0.1", "0.5", "1", "1e-800"), ("0.125", "0.001", "1e-30"))),
)
def test_independent_excluded_layer_quadrature(astr, estr):
    with mp.workdps(110):
        aa, ee = mp.mpf(astr), mp.mpf(estr)
        b = mp.mpf(1) / 16
        c = aa * ee
        ell = -mp.log(c)
        if c <= b:
            # Divide by c before comparing; tiny absolute errors cannot pass vacuously.
            scaled = mp.quad(lambda z: ell - mp.log(z) if z else 0, [0, 1])
            scaled += mp.quad(lambda z: z, [-mp.log(b), ell])
        else:
            scaled = b / c * mp.quad(lambda z: -mp.log(b * z) if z else 0, [0, 1])
        exact = endpoint.singular_layer_integral(
            s.Rational(astr), s.Rational(1, 8), s.Rational(estr)
        ) / (s.Rational(astr) * s.Rational(estr))
        assert abs(scaled - mp.mpf(str(exact.evalf(120)))) < mp.mpf("1e-90")
        assert 0 < scaled <= 1 + ell + ell * ell / 2
        assert aa * scaled <= 2 * (1 - mp.log(ee)) ** 2


@pytest.mark.parametrize(
    "astr,ystr", tuple(product(("0.1", "0.5", "1"), ("1", "1.5", "2", "2.5", "3")))
)
def test_independent_finite_Poisson_model_mark(astr, ystr):
    # Bounded-class scalar model only; not the original current or an L1 proof.
    with mp.workdps(100):
        aa, yy = mp.mpf(astr), mp.mpf(ystr)
        xx = mp.mpf("1e-100")
        ee = xx / yy
        i2 = mp.quad(lambda u: mp.log(yy - u) / u, [1, yy - 1]) if yy > 2 else mp.mpf(0)
        norm = 1 + aa * mp.log(yy) + aa * aa * i2 / 2

        def weighted(u):
            if u == yy:
                return mp.mpf(0)
            weight = 1 + aa * mp.log(u - 1) if u > 2 else 1
            return mp.log(ee * (yy - u)) * weight

        cuts = [1, yy] if yy <= 2 else [1, 2, yy]
        finite = aa * ee * mp.quad(weighted, cuts) / norm if yy > 1 else mp.mpf(0)
        full = -aa * xx / (aa + 1) * (-mp.log(xx) + mp.digamma(aa + 2) + mp.euler)
        assert abs(finite - full) <= 63000 * aa * ee * (1 - mp.log(ee)) ** 2


@pytest.mark.parametrize(
    "a,x",
    tuple(
        product(
            (s.Rational(1, 10**800), s.S.One),
            (s.Rational(1, 128), s.Rational(1, 10**10000)),
        )
    ),
)
def test_rare_index_extreme_resolution_without_cut_probability_division(a, x):
    eta = x / 4
    exact = density.common_space_cutoff_upper(a, x, eta) * source.KAPPA / (a * eta)
    with mp.workdps(100):
        xx = mp.mpf(str(x.evalf(110)))
        target = 63000 * (1 - mp.log(xx / 4)) ** 2
        assert abs(mp.mpf(str(exact.evalf(110))) / target - 1) < mp.mpf("1e-90")


@pytest.mark.parametrize("p", (s.Rational(1, 4), s.Rational(1, 2), s.Rational(3, 4)))
def test_both_density_normalization_terms_are_necessary(p):
    # Two-atom probability model with A_eta=Omega and constant mark=1.
    actual = p * abs(1 - 1 / p) + (1 - p)
    r = 1 - p
    assert actual == 2 * r and actual > r
    assert p * (1 - 1 / p) + (1 - p) == 0  # Means alone miss the whole L1 discrepancy.


def test_log_sum_subadditivity_has_strict_cross_term():
    y, w, v = s.Rational(1, 64), s.Rational(1, 128), s.Rational(1, 256)
    assert (1 + w / y) * (1 + v / y) > 1 + (w + v) / y


def test_configuration_law_TV_is_not_the_established_observable():
    assert "common underlying Poisson space" in audit.observable()["domain"]
    assert "emission-configuration laws" in audit.observable()["not_established"]
    with pytest.raises(ValueError):
        audit.require_observable("finite_configuration_laws_converge_in_TV")


def test_complete_historical_frontier_preserved():
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 188
    assert audit.matching()[:-1] == audit.previous.matching()
    assert len(audit.qualifications()) == 6
    assert "outer hard dimensional/evanescent" in audit.observable()["not_established"]
