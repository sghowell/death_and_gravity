"""Independent full-channel, Bose, soft-limit and explicit-domain checks."""

from functools import cache

import pytest
import sympy as s
from p8_vacuum_affine_minimal_gravity_radiation import audit, bounds, source, ward
from p8_vacuum_affine_minimal_gravity_radiation import vertices as v

ROWS, GATES = audit.residuals(), audit.gates()


@pytest.mark.parametrize("name", tuple(ROWS))
def test_all_exact_residuals(name):
    value = ROWS[name]
    assert all(
        x == 0 for x in (list(value) if isinstance(value, s.MatrixBase) else [value])
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
    with pytest.raises((ValueError, TypeError)):
        call(*args)


def substitutions(ps, k, xi):
    return {
        ward.mu: v.dot(ps[0], ps[0]),
        ward.a: v.dot(ps[0], ps[1]),
        ward.b: v.dot(ps[0], ps[2]),
        ward.u1: v.dot(ps[0], k),
        ward.u2: v.dot(ps[1], k),
        ward.u3: v.dot(ps[2], k),
        ward.x1: v.dot(ps[0], xi),
        ward.x2: v.dot(ps[1], xi),
        ward.x3: v.dot(ps[2], xi),
        ward.xk: v.dot(k, xi),
        s.Symbol("xx"): v.dot(xi, xi),
    }


@pytest.mark.parametrize("row", range(3))
@pytest.mark.parametrize("partition", range(3))
def test_independent_invariant_and_component_full_seven_graphs(row, partition):
    ps, k, _ = v.sample(row)
    left, right = v.PARTS[partition]
    ordered = tuple(ps[i] for i in (*left, *right))
    xi = s.Matrix([s.Rational(2, 3), 1, -s.Rational(3, 5), s.Rational(4, 7)])
    eps = (v.ETA * xi) * (v.ETA * xi).T
    sub = substitutions(ordered, k, xi)
    got = ward.rank_one_channel().xreplace(sub)
    expected = sum(v.channel(ps, k, eps, left, right))
    assert s.factor(got - expected) == 0


@pytest.mark.parametrize("row", range(3))
@pytest.mark.parametrize("permutation", ((1, 0, 2, 3), (2, 1, 0, 3), (0, 1, 3, 2)))
def test_full_tree_Bose_symmetry(row, permutation):
    ps, k, _ = v.sample(row)
    eps = s.diag(0, 1, -1, 0)
    assert (
        s.factor(
            v.amplitude(tuple(ps[i] for i in permutation), k, eps)
            - v.amplitude(ps, k, eps)
        )
        == 0
    )


@pytest.mark.parametrize("row", range(3))
def test_component_Born_matches_original_input(row):
    _, _, ps = v.sample(row)
    energy = v.dot(ps[0] + ps[1], ps[0] + ps[1])
    transfer = v.dot(ps[0] + ps[2], ps[0] + ps[2])
    cosine = 1 + 2 * transfer / (energy - 4)
    assert s.factor(v.born(ps) - source.forward.gravity_born(energy, cosine, 1)) == 0


@cache
def soft_row(denominator, polarization):
    h = s.Rational(1, denominator)
    alpha = 2 - h
    E = s.Rational(5, 4)
    ep = (alpha + 1 / alpha) / 2
    omega = E - ep**2 / E
    ps, k, p0 = v.recoil(E, omega, (0, 0, 1), (s.Rational(3, 5), s.Rational(4, 5), 0))
    eps = s.diag(0, 1, -1, 0) / s.sqrt(2)
    if polarization == "cross":
        eps = s.zeros(4)
        eps[1, 2] = eps[2, 1] = 1 / s.sqrt(2)
    full = v.amplitude(ps, k, eps)
    soft = sum((p.T * eps * p)[0] / v.dot(p, k) for p in p0)
    return (
        omega,
        s.factor(full - v.born(p0) * soft),
        s.factor(v.born(p0) * soft),
        v.born(p0),
    )


@pytest.mark.parametrize("polarization", ("plus", "cross"))
def test_exact_physical_soft_limit_and_finite_remainder(polarization):
    rows = [soft_row(n, polarization) for n in (1000, 10000, 100000)]
    for omega, remainder, leading, born in rows:
        assert abs(remainder / born) < bounds.REMAINDER
        assert omega < s.Rational(1, 192)
        assert abs(remainder / leading) < 100 * omega
    errors = [abs(rows[i + 1][1] - rows[i][1]) for i in range(2)]
    assert errors[1] < errors[0] / 5


@pytest.mark.parametrize(
    "resolution",
    (
        True,
        False,
        1.0,
        s.Float(1),
        "1",
        None,
        0,
        -1,
        s.oo,
        s.I,
        s.Symbol("x"),
        s.Rational(1, 100),
    ),
)
def test_invalid_or_excessive_resolution_rejected(resolution):
    with pytest.raises((TypeError, ValueError)):
        bounds.require_domain(9, -1, resolution)


@pytest.mark.parametrize(
    "tau", (s.Rational(1, 10), s.Rational(1, 10**204), s.Rational(1, 10**400))
)
def test_every_nonzero_transfer_has_admissible_resolution(tau):
    domain = bounds.require_domain(9, -tau, tau / 192)
    assert domain[-2:] == (tau, tau / 192)
    assert bounds.real_rate_bound(9, -tau, tau / 192, source.KAPPA) < s.Rational(
        1, 10**789
    )
    with pytest.raises(ValueError):
        bounds.require_domain(9, -tau, tau / 191)


def test_wrong_cubic_sign_has_nonzero_general_Ward_remainder():
    assert ward.data()["whole_wrong_cubic_sign_negative_control"] != 0


def test_no_original_status_or_parameter_change():
    assert audit.matching()[:-1] == audit.previous.matching()
    assert audit.frontier() == audit.previous.frontier()
    assert audit.qualifications() == audit.previous.qualifications()
    assert audit.require_parameters(audit.parameters()) == audit.parameters()


@pytest.mark.parametrize(
    "h", (s.Rational(1, 2), s.Rational(1, 10**6), s.Rational(1, 10**20))
)
@pytest.mark.parametrize(
    "alpha", (s.Rational(19, 10), 2 - s.Rational(1, 10**8), 2 - s.Rational(1, 10**30))
)
def test_exact_physical_gap_and_both_uniform_energy_regions(h, alpha):
    E = s.Rational(5, 4)
    ep = (alpha + 1 / alpha) / 2
    omega = E - ep**2 / E
    n = (s.Rational(3, 5), s.Rational(4, 5), 0)
    u = (2 * h / (1 + h * h), 0, (1 - h * h) / (1 + h * h))
    ps, k, p0 = v.recoil(E, omega, n, u)
    taus = (-v.dot(p0[0] + p0[2], p0[0] + p0[2]), -v.dot(p0[0] + p0[3], p0[0] + p0[3]))
    delta = min(s.S.One, *taus)
    for pair in ((0, 2), (1, 3), (0, 3), (1, 2)):
        D = -v.dot(ps[pair[0]] + ps[pair[1]], ps[pair[0]] + ps[pair[1]])
        Dj = -v.dot(p0[pair[0]] + p0[pair[1]], p0[pair[0]] + p0[pair[1]])
        assert D > (Dj + omega**2) / 30000
    e = s.Matrix([0, s.Rational(4, 5), -s.Rational(3, 5), 0])
    f = s.Matrix([0, 0, 0, 1])
    eps = (e * e.T - f * f.T) / 2
    full = v.amplitude(ps, k, eps)
    born = v.born(p0)
    soft = sum((p.T * eps * p)[0] / v.dot(p, k) for p in p0)
    assert abs(full / born) < 10**16 * delta / (omega * (delta + omega**2))
    assert abs(soft) < 160 * s.sqrt(delta) / omega
    if omega <= s.sqrt(delta) / 192:
        assert abs(full / born - soft) < 10**8 / s.sqrt(delta)


@pytest.mark.parametrize(
    "tau", (s.Rational(1, 10), s.Rational(1, 10**204), s.Rational(1, 10**400))
)
def test_fixed_resolution_not_shrunk_at_either_forward_endpoint(tau):
    x = s.Rational(1, 8)
    for transfer in (-tau, -5 + tau):
        assert bounds.require_fixed_resolution_domain(9, transfer, x)[-1] == x
        assert bounds.fixed_resolution_rate_bound(
            9, transfer, x, source.KAPPA
        ) < s.Rational(1, 10**768)
        with pytest.raises(ValueError):
            bounds.require_domain(9, transfer, x)


@pytest.mark.parametrize(
    "resolution",
    (
        True,
        False,
        1.0,
        s.Float(1),
        "1",
        None,
        0,
        -1,
        s.oo,
        s.I,
        s.Symbol("x"),
        s.Rational(1, 4),
    ),
)
def test_invalid_fixed_resolution_rejected(resolution):
    with pytest.raises((TypeError, ValueError)):
        bounds.require_fixed_resolution_domain(9, -1, resolution)


def test_uniform_bound_keeps_both_energy_budgets():
    packet = bounds.uniform_data()
    assert (
        packet["whole_combined_coefficient"]
        == packet["whole_low_band_coefficient"] + packet["whole_high_band_coefficient"]
    )
    assert packet["whole_combined_coefficient"] < 10**32
    assert all(value > 0 for value in packet["whole_uniform_positive_margins"].values())
