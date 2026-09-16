"""Independent full-state polarization, massless-limit and stability checks."""

from functools import cache

import mpmath as mp
import numpy as np
import pytest
import sympy as s
from p8_vacuum_affine_radiative_state_soft_index import audit, index, recoil, stability

ROWS, GATES = audit.residuals(), audit.gates()
SIGNS = (-1, -1, 1, 1)
DIRECTIONS = (
    (1, 0, 0),
    (-1, 0, 0),
    (0, 1, 0),
    (0, 0, 1),
    (s.Rational(3, 5), 0, s.Rational(4, 5)),
    (0, -s.Rational(4, 5), s.Rational(3, 5)),
    (-s.Rational(4, 5), s.Rational(3, 5), 0),
)
CASES = tuple((energy, count) for energy in ("5/4", "2") for count in (1, 2, 3, 5, 20))
U = s.Matrix([0, s.Rational(4, 5), s.Rational(3, 5)])


@pytest.mark.parametrize("name", tuple(ROWS))
def test_every_exact_identity(name):
    value = ROWS[name]
    assert all(
        x == 0 for x in (list(value) if isinstance(value, s.MatrixBase) else [value])
    )
    assert not s.sympify(value).has(s.Float)


@pytest.mark.parametrize("name", tuple(GATES))
def test_every_written_gate(name):
    assert GATES[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[r[0] for r in audit.bad_cases()]
)
def test_scope_rejections(name, call, args):
    assert name
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@cache
def make_state(energy, count, total="1/8"):
    E, R = s.Rational(energy), s.Rational(total)
    rays = tuple(
        s.Matrix(
            [R / count, *((R / count) * s.Matrix(DIRECTIONS[i % len(DIRECTIONS)]))]
        )
        for i in range(count)
    )
    return recoil.momenta(E, rays, U)


def mpvector(p):
    return [mp.mpf(str(s.N(x, 80))) for x in p]


def mdot(p, q):
    return p[0] * q[0] - sum(p[i] * q[i] for i in range(1, 4))


def direct_TT_index(points, nhat, nodes=100):
    nhat = np.array(nhat, dtype=float).reshape(3)
    e1 = np.array([1.0, 0.0, 0.0])
    if abs(np.dot(e1, nhat)) > 0.8:
        e1 = np.array([0.0, 1.0, 0.0])
    e1 -= np.dot(e1, nhat) * nhat
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(nhat, e1)
    z, weights = np.polynomial.legendre.leggauss(nodes)
    phi = 2 * np.pi * (np.arange(2 * nodes) + 0.5) / (2 * nodes)
    zz, pp = np.meshgrid(z, phi, indexing="ij")
    ss = np.sqrt(1 - zz * zz)
    ns = (
        ss[..., None] * np.cos(pp)[..., None] * e1
        + ss[..., None] * np.sin(pp)[..., None] * e2
        + zz[..., None] * nhat
    )
    t1 = (
        zz[..., None] * np.cos(pp)[..., None] * e1
        + zz[..., None] * np.sin(pp)[..., None] * e2
        - ss[..., None] * nhat
    )
    t2 = -np.sin(pp)[..., None] * e1 + np.cos(pp)[..., None] * e2
    jxx = np.zeros_like(zz)
    jxy = np.zeros_like(zz)
    jyy = np.zeros_like(zz)
    for p in np.array(points, dtype=float):
        den = p[0] - np.einsum("...i,i->...", ns, p[1:])
        xx = np.einsum("...i,i->...", t1, p[1:])
        yy = np.einsum("...i,i->...", t2, p[1:])
        jxx += xx * xx / den
        jxy += xx * yy / den
        jyy += yy * yy / den
    norm = (jxx - jyy) ** 2 / 2 + 2 * jxy * jxy
    return np.dot(weights, np.mean(norm, axis=1)) / 2


@pytest.mark.parametrize("energy,count", CASES)
def test_whole_state_index_against_independent_TT_integral(energy, count):
    points, rays, born = make_state(energy, count)
    with mp.workdps(80):
        actual = mp.mpf(
            str(s.N(index.kernel(points, rays).subs(index.source.MU, 1), 80))
        )
        base = mp.mpf(str(s.N(index.kernel(born, []).subs(index.source.MU, 1), 80)))
        assert actual >= 0
        assert abs(actual - base) < mp.mpf(1416) / 8 + mp.mpf(3) / 64
        direct = direct_TT_index(
            [mpvector(p) for p in (*points, *rays)], [0, 0, 1], nodes=200
        )
        assert abs(direct - float(actual)) < mp.mpf("1e-5")


@pytest.mark.parametrize("count", (1, 2, 5))
def test_full_current_Ward_and_massive_subset_failure(count):
    points, rays, _ = make_state("5/4", count)
    soft = s.Matrix([1, 0, -1, 0])

    def current(items):
        return sum((p * p.T / recoil.dot(p, soft) for p in items), s.zeros(4))

    all_current = current((*points, *rays))
    assert (all_current * recoil.ETA * soft).applyfunc(s.simplify) == s.zeros(4, 1)
    Q = sum(rays, s.zeros(4, 1))
    assert (current(points) * recoil.ETA * soft + Q).applyfunc(s.simplify) == s.zeros(
        4, 1
    )
    assert Q[0] > 0


@pytest.mark.parametrize("energy,count", (("5/4", 2), ("2", 3)))
def test_massless_pair_regulators_cancel(energy, count):
    points, rays, _ = make_state(energy, count)
    with mp.workdps(80):
        target = mp.mpf(
            str(s.N(index.kernel(points, rays).subs(index.source.MU, 1), 80))
        )
        ps = [[SIGNS[i] * v for v in mpvector(p)] for i, p in enumerate(points)] + [
            mpvector(q) for q in rays
        ]
        signs = [*SIGNS, *([1] * count)]
        errors = []
        for reg in (mp.mpf("1e-8"), mp.mpf("1e-16"), mp.mpf("1e-24")):
            masses = [mp.mpf(1)] * 4 + [reg**2 * q[0] ** 2 for q in ps[4:]]
            result = sum(masses) / 2
            for i, p in enumerate(ps):
                for j in range(i + 1, len(ps)):
                    d = mdot(p, ps[j])
                    uv = masses[i] * masses[j]
                    result += (
                        2
                        * signs[i]
                        * signs[j]
                        * (d * d - uv / 2)
                        * mp.acosh(d / mp.sqrt(uv))
                        / mp.sqrt(d * d - uv)
                    )
            errors.append(abs(result - target))
        assert errors[2] < errors[1] < errors[0]
        assert errors[2] < mp.mpf("1e-40")


@pytest.mark.parametrize("energy,count", (("5/4", 3), ("2", 5)))
def test_exact_rational_boost_reference_cancellation(energy, count):
    points, rays, _ = make_state(energy, count)
    boost = s.Matrix(
        [
            [s.Rational(5, 4), -s.Rational(3, 4), 0, 0],
            [-s.Rational(3, 4), s.Rational(5, 4), 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1],
        ]
    )
    original = index.kernel(points, rays).subs(index.source.MU, 1)
    moved = index.kernel(
        tuple(boost * p for p in points), tuple(boost * q for q in rays)
    ).subs(index.source.MU, 1)
    assert abs(s.N(original - moved, 70)) < s.Rational(1, 10**55)


@pytest.mark.parametrize("count", (1, 2, 5))
def test_arbitrary_collinear_partition(count):
    points, rays, _ = make_state("5/4", count)
    split = (*rays[:-1], *(rays[-1] / (count + 1) for _ in range(count + 1)))
    change = (index.kernel(points, split) - index.kernel(points, rays)).subs(
        index.source.MU, 1
    )
    assert abs(s.N(change, 70)) < s.Rational(1, 10**55)


@pytest.mark.parametrize("total", ("1/8", "1/1000", "1/1000000"))
def test_known_power_and_Gamma_factor_after_pairing(total):
    points, rays, born = make_state("5/4", 3, total)
    with mp.workdps(80):
        k = mp.mpf(1000)
        a = mp.mpf(
            str(s.N(index.kernel(points, rays).subs(index.source.MU, 1), 80))
        ) / (4 * mp.pi**2 * k)
        a0 = mp.mpf(str(s.N(index.kernel(born, []).subs(index.source.MU, 1), 80))) / (
            4 * mp.pi**2 * k
        )
        x = mp.mpf(str(s.N(s.Rational(total), 80)))
        logratio = (
            (a - a0) * mp.log(x)
            - mp.euler * (a - a0)
            - mp.loggamma(1 + a)
            + mp.loggamma(1 + a0)
        )
        assert abs(logratio) < 12 / k + 60 / k**2
        assert abs(mp.expm1(logratio)) < 26 / k


BAD_STATES = (
    (True, [], (1, 0, 0)),
    (1.5, [], (1, 0, 0)),
    (1, [], (1, 0, 0)),
    (3, [], (1, 0, 0)),
    (s.Rational(5, 4), [(0, 0, 0, 0)], (1, 0, 0)),
    (s.Rational(5, 4), [(-1, -1, 0, 0)], (1, 0, 0)),
    (s.Rational(5, 4), [(s.Rational(1, 10), 0, 0, 0)], (1, 0, 0)),
    (s.Rational(5, 4), [(s.Rational(1, 4), s.Rational(1, 4), 0, 0)], (1, 0, 0)),
    (s.Rational(5, 4), [], (2, 0, 0)),
    (s.Symbol("energy"), [], (1, 0, 0)),
)


@pytest.mark.parametrize("args", BAD_STATES)
def test_exact_recoil_domain_guards(args):
    with pytest.raises((TypeError, ValueError)):
        recoil.momenta(*args)


def test_zero_radiation_and_original_scope():
    points, rays, born = recoil.momenta(s.Rational(5, 4), [], U)
    assert points == born and rays == ()
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 156
    assert audit.matching()[:-1] == audit.previous.matching()
    assert len(audit.qualifications()) == 6
    assert stability.original_factor_ratio_error() == s.Rational(26, 10**800)
    assert (
        "full d finite angular conversion"
        in audit.observable()["not_established"].lower()
    )


@pytest.mark.parametrize("direction", (-1, 1))
def test_exact_zero_angle_zero_radiation_endpoint(direction):
    points, rays, _ = recoil.momenta(s.Rational(5, 4), [], (0, 0, direction))
    assert s.simplify(index.kernel(points, rays)) == 0
