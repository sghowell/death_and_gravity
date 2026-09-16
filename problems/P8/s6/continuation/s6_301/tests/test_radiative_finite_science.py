"""Independent finite-dimensional, angular, collinear and scope checks."""

from functools import cache

import mpmath as mp
import numpy as np
import pytest
import sympy as s
from p8_vacuum_affine_radiative_angular_finite import (
    angular,
    audit,
    bounds,
    conversion,
    source,
)
from scipy.special import roots_jacobi

recoil, index = source.recoil, source.index
ROWS, GATES = audit.residuals(), audit.gates()
U = s.Matrix([0, s.Rational(4, 5), s.Rational(3, 5)])
DIRECTIONS = [
    (1, 0, 0),
    (-1, 0, 0),
    (0, 1, 0),
    (0, 0, 1),
    (s.Rational(3, 5), 0, s.Rational(4, 5)),
]


def make(E, count, total=None):
    total = s.Rational(1, 8) if total is None else total
    rays = [
        s.Matrix([total / count, *((total / count) * s.Matrix(DIRECTIONS[i % 5]))])
        for i in range(count)
    ]
    points, rays, born = recoil.momenta(s.Rational(E), rays, U)
    return np.array([list(p) for p in (*points, *rays)], float), points, rays, born


def grid(nodes):
    z, w = np.polynomial.legendre.leggauss(nodes)
    phi = 2 * np.pi * (np.arange(2 * nodes) + 0.5) / (2 * nodes)
    zz, pp = np.meshgrid(z, phi, indexing="ij")
    n = np.stack(
        (np.sqrt(1 - zz * zz) * np.cos(pp), np.sqrt(1 - zz * zz) * np.sin(pp), zz),
        axis=-1,
    ).reshape(-1, 3)
    weight = np.repeat(w / (4 * nodes), 2 * nodes)
    return n, weight


def current_moments(ps, n, t):
    radial = np.sqrt(t)
    nhat = np.column_stack((radial * n, np.full(len(n), np.sqrt(1 - t))))
    current = np.zeros((len(n), 4, 4))
    for p in ps:
        pp = np.array([*p[1:], 0.0])
        projected = pp - nhat * np.einsum("bi,i->b", nhat, pp)[:, None]
        den = p[0] - radial * np.einsum("bi,i->b", n, p[1:])
        current += np.einsum("bi,bj->bij", projected, projected) / den[:, None, None]
    trace = np.trace(current, axis1=1, axis2=2)
    trace2 = trace * trace
    return np.einsum("bij,bij->b", current, current) - trace2 / 2, trace2


def evaluate(ps, angular=52, radial=56, eps=(0.02, 0.01, 0.005)):
    n, weight = grid(angular)
    f0, tr0 = current_moments(ps, n, 1.0)
    k0 = float(weight @ f0)
    t0 = float(weight @ tr0)
    nodes, w = np.polynomial.legendre.leggauss(radial)
    u = (nodes + 1) / 2
    ww = w / 2
    integral = 0.0
    for x, wt in zip(u, ww):
        t = (1 - x * x) ** 2
        f, _ = current_moments(ps, n, t)
        integral += wt * 4 * (1 - x * x) ** 2 * (weight @ (f - f0)) / (x * (2 - x * x))
    k1 = t0 / 2 + integral
    values = []
    for e in eps:
        # Independent normalized beta quadrature; do not form a tiny difference
        # using a subtraction inside the radial quadrature.
        z, wj = roots_jacobi(radial, e - 1, 0.5)
        weights = wj / np.sum(wj)
        ke = 0.0
        for x, wt in zip((z + 1) / 2, weights):
            f, tr = current_moments(ps, n, x)
            ke += wt * (weight @ (f + e * tr / (2 * (1 + e))))
        values.append((e, float(ke), float((ke - k0) / e)))
    return k0, k1, t0, values


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
def test_original_scope_rejections(name, call, args):
    assert name
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@cache
def evaluated(energy, count):
    ps, points, rays, born = make(energy, count)
    return evaluate(ps), points, rays, born


CASES = (("5/4", 0), ("2", 0), ("5/4", 1), ("5/4", 3), ("2", 5), ("2", 20))


@pytest.mark.parametrize("energy,count", CASES)
def test_full_radiative_dimensional_limit_and_finite_derivative(energy, count):
    (k0, k1, trace, values), points, rays, _ = evaluated(energy, count)
    exact = float(s.N(index.kernel(points, rays), 35))
    assert abs(k0 - exact) < 2e-4
    assert 0 <= k0 < 265 and abs(k1) < 130000 and 0 <= trace < 265
    errors = [abs(row[2] - k1) for row in values]
    assert errors[2] < errors[1] < errors[0]
    for e, ke, _ in values:
        assert abs(ke - k0 - e * k1) < 1400000 * e * e


def massive_pair_first(E, c):
    def ratio(b):
        return (
            1 + b / 3 + b * b / 5
            if abs(b) < mp.mpf("1e-30")
            else mp.atanh(mp.sqrt(b)) / mp.sqrt(b)
        )

    value = 4 * (mp.mpf(".5") + 1 - ratio(1 - 1 / E**2))
    for sign, d in zip(
        (1, -1, -1), (2 * E**2 - 1, 1 + (E**2 - 1) * (1 - c), 1 + (E**2 - 1) * (1 + c))
    ):
        j0 = mp.acosh(d) / mp.sqrt(d * d - 1)

        def first(x, d=d):
            mx = 1 + 2 * x * (1 - x) * (d - 1)
            return 2 * (1 - ratio(1 - mx / E**2)) / mx

        value += (
            4
            * sign
            * (j0 / 2 + (d * d - mp.mpf(".5")) * mp.quad(first, [0, mp.mpf(".5"), 1]))
        )
    return value


@pytest.mark.parametrize("energy", ("5/4", "2"))
def test_independent_massive_Feynman_parameter_derivative(energy):
    (_, k1, _, _), _, _, _ = evaluated(energy, 0)
    with mp.workdps(40):
        E = mp.mpf(str(s.N(s.Rational(energy), 50)))
        expected = massive_pair_first(E, mp.mpf(".6"))
    assert abs(k1 - float(expected)) < 2e-6


@pytest.mark.parametrize("tau", (1.0, 0.1, 0.001, 1e-6, 1e-10))
def test_independent_massless_nuclear_norm_Holder_calibration(tau):
    z, w = np.polynomial.legendre.leggauss(180)
    n = np.column_stack((np.sqrt(1 - z * z), np.zeros(len(z)), z))
    nt = np.column_stack((np.sqrt(1 - tau) * n, np.full(len(n), np.sqrt(tau))))
    n1 = np.column_stack((n, np.zeros(len(n))))
    p = np.array([0.0, 0.0, 1.0, 0.0])
    vt = p - nt * (nt @ p)[:, None]
    v1 = p - n1 * (n1 @ p)[:, None]
    At = np.einsum("bi,bj->bij", vt, vt) / (1 - np.sqrt(1 - tau) * z)[:, None, None]
    A1 = np.einsum("bi,bj->bij", v1, v1) / (1 - z)[:, None, None]
    result = np.dot(w, np.linalg.svd(At - A1, compute_uv=False).sum(axis=1)) / 2
    assert result <= 14 * tau**0.25


@pytest.mark.parametrize("radial", (0, s.Rational(3, 5), 1))
def test_exact_collinear_split_preserves_entire_current(radial):
    q = s.Matrix([s.Rational(1, 8), s.Rational(1, 8), 0, 0])
    first = angular.transverse_current(s.Rational(5, 4), [q], U, (0, 1, 0), radial)
    split = angular.transverse_current(
        s.Rational(5, 4), [q / 7] * 7, U, (0, 1, 0), radial
    )
    assert (first - split).applyfunc(s.simplify) == s.zeros(4)


@pytest.mark.parametrize("outgoing", ((0, 0, 1), (0, 0, -1)))
def test_exact_forward_backward_zero_radiation(outgoing):
    for radial in (0, s.Rational(3, 5), 1):
        current = angular.transverse_current(
            s.Rational(5, 4), [], outgoing, (0, 1, 0), radial
        )
        assert current == s.zeros(4)


def test_bounded_null_collinear_boundary_representative():
    q = s.Matrix([s.Rational(1, 8), s.Rational(1, 8), 0, 0])
    current = angular.transverse_current(s.Rational(5, 4), [q], U, (1, 0, 0), 1)
    value, trace = angular.contractions(current)
    assert not current.has(s.nan, s.zoo, s.oo)
    assert 0 <= float(value) <= 265 and 0 <= float(trace) <= 265


BAD = (
    (angular.require_direction, ((1, 0),)),
    (angular.require_direction, ((1.0, 0, 0),)),
    (angular.require_direction, ((1, 1, 0),)),
    (angular.require_radial, (1.0,)),
    (angular.require_radial, (-1,)),
    (angular.require_radial, (2,)),
    (bounds.require_epsilon, (0.1,)),
    (bounds.require_epsilon, (-1,)),
    (bounds.require_epsilon, (s.Rational(1, 4),)),
    (conversion.require_resolution, (0,)),
    (conversion.require_resolution, (-1,)),
    (conversion.require_resolution, (s.Rational(1, 4),)),
    (conversion.require_resolution, (0.1,)),
    (bounds.logarithmic_holder_moment, (True,)),
    (bounds.logarithmic_holder_moment, (-1,)),
    (bounds.logarithmic_holder_moment, (1.0,)),
)


@pytest.mark.parametrize("call,args", BAD)
def test_domain_guards(call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize("power", range(5))
def test_logarithmic_Holder_integral_independent_quadrature(power):
    with mp.workdps(35):
        # tau=exp(-u) removes the endpoint singularity.
        result = mp.quad(lambda u: mp.exp(-u / 4) * u**power, [0, 4, 16, mp.inf])
        expected = int(bounds.logarithmic_holder_moment(power))
        assert abs(result - expected) < mp.mpf("1e-25")


@pytest.mark.parametrize("energy,count", (("5/4", 1), ("2", 5)))
def test_finite_conversion_and_fixed_resolution_error(energy, count):
    (k0, k1, _, values), _, _, _ = evaluated(energy, count)
    with mp.workdps(35):
        phase_first = mp.euler - 2 - mp.log(mp.pi)
        delta = (k1 + phase_first * k0) / (8 * mp.pi**2)
        assert abs(delta) < 2000
        x = mp.mpf(1) / 8
        for epsilon, ke, _ in values:
            e = mp.mpf(str(epsilon))
            phase = (
                (4 * mp.pi) ** (-e)
                * mp.gamma(mp.mpf("1.5"))
                / mp.gamma(mp.mpf("1.5") + e)
            )
            paired = x ** (2 * e) * (phase * ke - k0) / (8 * mp.pi**2 * e)
            assert abs(paired - delta) < e * (28000 + 4000 * abs(mp.log(x)))


def test_original_bounds_and_full_scope_preserved():
    assert conversion.original_conversion_bound() == s.Rational(2000, 10**800)
    assert conversion.original_reference_ratio_bound() == s.Rational(4250, 10**800)
    assert conversion.original_unexpanded_factor_remainder() == s.Rational(
        5000000, 10**1600
    )
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 157
    assert audit.matching()[:-1] == audit.previous.matching()
    assert len(audit.qualifications()) == 6
    assert "hard evanescent" in audit.observable()["not_established"].lower()
