"""Independent signed-series, marked-energy and original-source checks."""

from functools import cache
from itertools import product

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_affine_single_residual_soft_dressing import (
    audit,
    bounds,
    dressing,
    series,
    source,
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
def test_written_proof_gates(name):
    assert GATES[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[r[0] for r in audit.bad_cases()]
)
def test_original_scope_rejections(name, call, args):
    assert name
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize("total", range(1, 5))
def test_independent_state_dependent_marked_Bose_sum(total):
    states = (s.Rational(1, 9), s.Rational(1, 7), s.Rational(1, 5))
    cut = s.Rational(3, 5)
    residual = lambda w: 1 - 7 * w
    soft = lambda mark, w: 1 + mark + 2 * w
    left = s.S.Zero
    for row in product(states, repeat=total):
        if sum(row) <= cut:
            left += sum(
                residual(row[j])
                * s.prod(soft(row[j], row[i]) for i in range(total) if i != j)
                for j in range(total)
            ) / s.factorial(total)
    right = sum(
        residual(mark) * s.prod(soft(mark, w) for w in row) / s.factorial(total - 1)
        for mark in states
        for row in product(states, repeat=total - 1)
        if mark + sum(row) <= cut
    )
    assert left == right


@pytest.mark.parametrize("power", (1, 2, 3))
@pytest.mark.parametrize("count", (0, 1, 2))
def test_literal_marked_energy_integral(power, count):
    with mp.workdps(55):
        x = mp.mpf(1) / 7
        e = mp.mpf(1) / 8
        ae = mp.mpf(3) / 11
        got = mp.quad(
            lambda u: (
                x**power
                * u ** (power - 1)
                * (ae * mp.gamma(2 * e) * (x * (1 - u)) ** (2 * e)) ** count
                / (mp.factorial(count) * mp.gamma(1 + 2 * e * count))
            ),
            [0, 1],
        )
        target = series.marked_power_term(
            count, power, s.Rational(3, 11), s.Rational(1, 8), s.Rational(1, 7)
        )
        assert abs(got / mp.mpf(str(s.N(target, 50))) - 1) < mp.mpf("1e-40")


def finite_sum(a, d, e, x, p):
    lam = (a + 2 * e * d) * mp.gamma(2 * e) * x ** (2 * e)
    limit = int(abs(lam) + 25 * mp.sqrt(abs(lam) + 1) + 180)
    weight = mp.mpf(1)
    terms = []
    for n in range(limit + 1):
        terms.append(weight / mp.gamma(p + 1 + 2 * e * n))
        weight *= lam / (n + 1)
    assert abs(terms[-1]) < max(mp.mpf(1), max(map(abs, terms))) * mp.mpf("1e-60")
    return mp.gamma(p) * x**p * mp.exp(-a / (2 * e)) * mp.fsum(terms)


@pytest.mark.parametrize("physical_index,delta", (("0.2", "0.07"), ("0", "-0.3")))
@pytest.mark.parametrize("power", ("0.5", "1", "2"))
def test_complete_dimensional_series_limit_including_signed_zero_index(
    physical_index, delta, power
):
    with mp.workdps(75):
        a, d, p = map(mp.mpf, (physical_index, delta, power))
        x = mp.mpf(1) / 7
        target = (
            mp.gamma(p) * mp.exp(d - mp.euler * a) * x ** (p + a) / mp.gamma(p + 1 + a)
        )
        errors = [
            abs(finite_sum(a, d, mp.mpf(e), x, p) / target - 1)
            for e in ("0.002", "0.0005", "0.000125")
        ]
        assert errors[1] < errors[0] / 3 and errors[2] < errors[1] / 3
        assert errors[-1] < mp.mpf("0.002")


@pytest.mark.parametrize("count", (1, 3, 5))
def test_negative_finite_dimensional_index_is_not_rejected_as_probability(count):
    term = series.finite_term(
        count, -s.Rational(1, 50), s.Rational(1, 8), s.Rational(1, 4)
    )
    assert term.is_negative is True


@pytest.mark.parametrize("power", range(1, 7))
def test_individual_energy_cut_is_not_total_energy_cut(power):
    a = s.Rational(1, 5)
    ratio = s.factorial(power) / s.prod(a + j for j in range(1, power + 1))
    assert 0 < ratio < 1
    assert ratio != 1


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
def test_invalid_soft_regulator_rejected(invalid):
    with pytest.raises((TypeError, ValueError)):
        series.require_epsilon(invalid)


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
        2,
    ),
)
def test_invalid_remaining_energy_rejected(invalid):
    with pytest.raises((TypeError, ValueError)):
        series.require_remaining(invalid)


@pytest.mark.parametrize("count", (True, False, 1.0, s.Float(1), -1, s.Rational(1, 2)))
def test_invalid_emission_count_rejected(count):
    with pytest.raises((TypeError, ValueError)):
        series.finite_term(count, s.Rational(1, 10), s.Rational(1, 8), s.Rational(1, 8))


@pytest.mark.parametrize("power", (0, -1, True, 1.0))
def test_nonintegrable_or_inexact_marked_power_rejected(power):
    with pytest.raises((TypeError, ValueError)):
        series.marked_power_term(
            1, power, s.Rational(1, 10), s.Rational(1, 8), s.Rational(1, 8)
        )


@pytest.mark.parametrize(
    "resolution", (s.Rational(1, 8), s.Rational(1, 10**20), s.exp(-10 * source.KAPPA))
)
def test_unexpanded_original_remaining_energy_weight_at_extreme_threshold(resolution):
    a = s.Rational(1, 2) / source.KAPPA
    weight = dressing.original_weight(a, 0, a, 0, resolution / 2, resolution)
    assert s.simplify(weight - 2 ** (-a)) == 0
    assert weight.is_positive is True


@pytest.mark.parametrize(
    "mutation",
    (
        "negative_index",
        "wrong_index_gap",
        "bad_conversion",
        "zero_remaining",
        "too_large_resolution",
    ),
)
def test_original_weight_domain_is_enforced(mutation):
    k = source.KAPPA
    args = [
        s.Rational(1, 2) / k,
        0,
        s.Rational(1, 2) / k,
        0,
        s.Rational(1, 64),
        s.Rational(1, 8),
    ]
    if mutation == "negative_index":
        args[0] = -1 / k
    if mutation == "wrong_index_gap":
        args[0] = 2 / k
    if mutation == "bad_conversion":
        args[1] = 2001 / k
    if mutation == "zero_remaining":
        args[4] = args[5]
    if mutation == "too_large_resolution":
        args[5] = s.Rational(1, 4)
    with pytest.raises((TypeError, ValueError)):
        dressing.original_weight(*args)


@cache
def physical(row):
    _ps, q, p0 = source.vertices.sample(row)
    e = -p0[0][0]
    u = p0[2][1:4, 0] / s.sqrt(e * e - 1)
    return dressing.physical_kernel(e, q[0], q[1:4, 0] / q[0], u)


@pytest.mark.parametrize("row", range(3))
def test_complete_original47_tree_residual_is_actually_signed(row):
    packet = physical(row)
    assert packet["whole_Born"] > 0
    assert packet["kappa_times_energy_weighted_signed_density"] < 0
    assert 0 < packet["rho"] < 1
    assert not packet["kappa_times_energy_weighted_signed_density"].has(s.Float)


@pytest.mark.parametrize("row", range(3))
def test_original_physical_radiative_index_and_recoil(row):
    packet = physical(row)
    q = packet["ray"]
    assert all(s.simplify(v) == 0 for v in sum(packet["points"], q.copy()))
    assert source.vertices.dot(q, q) == 0
    assert s.N(packet["whole_radiative_K0"], 30) > 0
    assert (
        abs(s.N(packet["whole_radiative_K0"] - packet["whole_elastic_K0"], 30))
        < 1440 * q[0]
    )


@pytest.mark.parametrize(
    "resolution", (s.Rational(1, 8), s.Rational(1, 10**10), s.Rational(1, 10**400))
)
@pytest.mark.parametrize(
    "transfer", (-s.Rational(1, 10), -s.Rational(1, 10**204), -s.Rational(1, 10**400))
)
def test_uniform_dressed_bound_without_angular_or_resolution_window(
    resolution, transfer
):
    sharp = bounds.relative_bound(9, transfer, resolution)
    assert sharp < bounds.coarse_relative_bound(9, transfer, resolution)
    assert sharp < bounds.vanishing_relative_bound(9, transfer, resolution)
    assert sharp < 1
    assert bounds.coarse_relative_bound(9, transfer, resolution) == s.Rational(
        2, 10**768
    )


def test_wrong_virtual_state_has_nonzero_IR_pole():
    e = source.EP
    assert dressing.data()["whole_wrong_virtual_negative_control"] != 0
    assert e in dressing.data()["whole_wrong_virtual_negative_control"].free_symbols


def test_all_original_scope_and_rejected_history_preserved():
    assert audit.matching()[:-1] == audit.previous.matching()
    assert audit.frontier() == audit.previous.frontier()
    assert audit.qualifications() == audit.previous.qualifications()
    assert audit.require_parameters(audit.parameters()) == audit.parameters()
    assert "full all-N nonleading" in audit.observable()["not_established"]
